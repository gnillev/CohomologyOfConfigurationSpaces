import pprint
load("cohomology.sage")
load("algebra.sage")


G = 2; # Genus of complex curve
N = 3; # Number of points in configuration space.
L = 1; # Complex dimension; No impact at the moment.

### Set genus for all instances of CohomologyBasisElement in the following
CohomologyBasisElement._genus = G


# cohomXn = CohomologyObject([1,2*G,1],"h")**N
# confAlgebraRn = ConfAlgebraObject(N);

# print cohomXn;
# print confAlgebraRn;


class dgBasisTuple(object):
    """docstring for dgBasisTuple"""
    def __init__(self, element, coord = 1):
        self.element = element
        self.coord = coord
        self.hpart = CohomologyBasisElement(element[0])
        self.gpart = AlgebraElement([AlgebraBasisElement(element[1])])

        i = self.hpart.degree
        j = self.gpart.degree
        self.multiIndex = (i,j)

    def __str__(self):
        return str(self.coord) + "*" + str(self.element)

    def __repr__(self):
        return str(self)

    def __mul__(a,b):
        # We multiply some ([(d_1,k_1),...,(d_i,k_i)],[(a_1, b_1), ..., (a_j, b_j)]) with ([(d'_1,k'_1),...,(d'_i',k'_i')],[(a'_1, b'_1), ..., (a'_j', b'_j')]).
        # First commute algebra generators from a with cohomology generators from b.
        # Calculate signed coordinate
        sign = (-1)**(a.multiIndex[1] * b.multiIndex[0])

        # Multiply cohomology generators
        hpart = a.hpart * b.hpart
        # Multiply algebra generators
        
        gpart = a.gpart * b.gpart

        # Get base coefficient
        coord = sign * a.coord * b.coord * hpart.coord 
        
        
        # Replace a[j]*G_ij with a[i]*G_ij
        helm = hpart.element
        elements = []
        if coord != 0:
            for gfactor in gpart.elements:
                if gfactor.coord != 0:
                    out_coord = coord * gfactor.coord
                    if gfactor.element != 1:
                        for g in gfactor:
                            (i,j) = (g[0], g[1])
                            if helm[j] != (0,0):
                                hsub = CohomologyBasisElement([helm[i]]) * CohomologyBasisElement([helm[j]])

                                # Factor in sign from "switching" places
                                sub_sign = 1
                                for k in range(i+1,j):
                                    sub_sign *= (-1)**(helm[j][0] * helm[k][0])

                                helm[i] = hsub.element[0]
                                helm[j] = (0,0)
                                out_coord *= sub_sign * hsub.coord
                    if out_coord != 0:
                        elements.append(dgBasisTuple((helm, gfactor.element), out_coord))

        
        return dgElement(elements)

    def __rmul__(self, other):
        return dgBasisTuple(self.element, self.coord * other)


class dgElement(object):
    """docstring for dgElement"""
    def __init__(self, elements):
        self.elements = elements;
        
        if elements == []:
            self.multiIndex = (0,0)
        else:
            if all(x.multiIndex == elements[0].multiIndex for x in elements):
                self.multiIndex = elements[0].multiIndex
            else:
                raise ValueError("Basiselements must be of the same bi-degree.")

    def __add__(a,b):
        if a.multiIndex == b.multiIndex:
            return dgElement(a.elements + b.elements)
        else:
            raise ValueError("To add two dgElements, they must reside in the same vector space.")


    def __mul__(a,b):
        if a == 0 or b == 0:
            return 0;

        out_elements = []
        for aelm in a.elements:
            for belm in b.elements:
                out_elements += aelm * belm
        
        return dgElement(out_elements)

    def __rmul__(self, other):
        return dgElement(map(lambda elm : other * elm, self.elements))

    def __str__(self):
        return str(self.elements);

    def __repr__(self):
        return str(self)

    def __getitem__(self, item):
        return self.elements[item];



class CohomConfSpaceComplexCurve(object):
    """docstring for ConfSpaceCohomComplexCurve"""
    def __init__(self, npoints, genus):
        self.npoints = npoints
        self.genus = genus
        self.complex_dimension = 1

        self.cohomObject = ComplexCurveCohomPower(genus,npoints)
        self.algebraObject = ConfAlgebraObject(npoints)

        self.E2 = [];
        self.E2Vector = [];

        for i,h in enumerate(self.cohomObject):
            self.E2.append([]);
            self.E2Vector.append([]);
            for j,g in enumerate(self.algebraObject):
                self.E2[i].append([]);
                for x in h:
                    for y in g:
                        if self.relGH1(x,y):
                            self.E2[i][j].append((x,y));

                self.E2Vector[i].append(FreeModule(QQ,len(self.E2[i][j])));

    def GetModule(self,i,j):
        if i < 0 or j < 0 or i >= len(self.E2) or j >= len(self.E2[i]):
            return FreeModule(QQ,0)
        else:
            return self.E2Vector[i][j]

    def GetTensor(self,i,j):
        if i < 0 or j < 0 or i >= len(self.E2) or j >= len(self.E2[i]):
            return []
        else:
            return self.E2[i][j]

    def GetVectorFromElement(self, elems):
        (i,j) = elems.multiIndex
        out_vector = 0
        for elem in elems:
            indx = self.E2[i][j].index(elem.element)
            out_vector += elem.coord * self.E2Vector[i][j].basis()[indx] 
        return out_vector

    def GetElementFromVector(self, vect, i, j):
        out_elements = []
        for indx, coord in enumerate(self.E2Vector[i][j].coordinates(vect)):
            if coord != 0:
                out_elements.append(dgBasisTuple(self.E2[i][j][indx],coord))
        return dgElement(out_elements)

    def relGH1(self, x, y):
        if y != 1:
            for gelm in y:
                for j in xrange(0,N):
                    if x[j] != (0,0) and gelm[1] == j:
                        return False
        return True     

    def d_on_vectorbasis(self, x, i, j):
        if (i,j) <= (0,0):
            return 0;
        elif j == 0:
            return 0;
        elif (i,j) == (0,1):
            """ Send G_{a,b} to class of diagonal in (0+2*L, 1+1-2*L) = (2,0) """
            cohom_indx = self.E2Vector[i][j].basis().index(x)
            cohom_elm = self.E2[i][j][cohom_indx]

            diagonal_sum = self.cohomObject.Diagonal(cohom_elm[1][0][0],cohom_elm[1][0][1])

            out_vector = 0
            for elm in diagonal_sum:
                out_vector += self.cohomObject.GetVector(elm)
            return out_vector
        else:
            inElement = self.GetElementFromVector(x,i,j);
            i_out = i+2
            j_out = j-1
            hpart = inElement[0].hpart
            gpart = inElement[0].gpart
            
            
            g0Index = self.algebraObject[1].index([gpart[0][0]]);
            
            if j == 1:
                grIndex = 0;
            else:
                grIndex = self.algebraObject[j-1].index(gpart[0][1:]);
            
            
            hIndex = self.cohomObject[i].index(hpart.element)
            hElement = dgElement([dgBasisTuple(self.E2[i][0][hIndex])])
            
            g0Element = dgElement([dgBasisTuple(self.E2[0][1][g0Index])])
            g0Vector = self.GetVectorFromElement(g0Element)
            grElement = dgElement([dgBasisTuple(self.E2[0][j-1][grIndex])])
            grVector = self.GetVectorFromElement(grElement)

            #print inElement, g0Element, grElement
            
            d1Element = self.GetElementFromVector(self.d_on_vectorbasis(g0Vector,0,1),2,0) * grElement
            d1Vector = self.GetVectorFromElement(d1Element)
            
            d2Element = g0Element * self.GetElementFromVector(self.d_on_vectorbasis(grVector,0,j-1),2,j-2) 
            d2Vector = self.GetVectorFromElement(d2Element)

            sign = (-1)**(hElement.multiIndex[0])
            

            outElement =  (sign * hElement) * self.GetElementFromVector(d1Vector - d2Vector, 2, j-1) ;

            outVector = self.GetVectorFromElement(outElement)
            
            return outVector

    # Must be able to evaluate on combinations of basis-elements
    def d_on_vector(self, x, i, j):
        if (i,j) <= (0,0):
            return 0;
        elif j == 0:
            return 0;
        elif (i,j) == (0,1):
            """ Send G_{a,b} to class of diagonal in (0+2*L, 1+1-2*L) = (2,0) """
            
            inElement = self.GetElementFromVector(x,i,j)
            
            vectorCombination = []
            for basis_elm in inElement:
                
                diagonal_sum = self.cohomObject.Diagonal(basis_elm.gpart[0][0][0],basis_elm.gpart[0][0][1])

                out_vector = 0
                for elm in diagonal_sum:
                    out_vector += self.cohomObject.GetVector(elm)
                vectorCombination.append(out_vector)

            return reduce(lambda x,y : x+y, vectorCombination)
        else:
            inElement = self.GetElementFromVector(x,i,j);
            i_out = i+2
            j_out = j-1

            vectorCombination = []
            #Loop through the linear combination
            for basis_elm in inElement:

                #Fetch the cohomObject and algebraObject from the basistuple
                hpart = basis_elm.hpart
                gpart = basis_elm.gpart
                
                #print basis_elm, hpart, gpart
                
                #Fetch the index of the first G_ab-element in the product
                g0Index = self.algebraObject[1].index([gpart[0][0]]);
                #Fetch the corresponding element embedded in E2[0][1]
                g0Element = dgElement([dgBasisTuple(self.E2[0][1][g0Index])])
                #Convert to vector
                g0Vector = self.GetVectorFromElement(g0Element)


                # If j>1 there are at least one more G_ab factor
                if j == 1:
                    grIndex = 0;
                else:
                    grIndex = self.algebraObject[j-1].index(gpart[0][1:]);
                
                #Fetch the corresponding element embedded in E2[0][j-1] (products of G_ij with j-1 factors)
                grElement = dgElement([dgBasisTuple(self.E2[0][j-1][grIndex])])
                #Convert to vector
                grVector = self.GetVectorFromElement(grElement)

                #Fetch the index of the hpart
                hIndex = self.cohomObject[i].index(hpart.element)
                #Fetch the element embedded in E2[i][0] and attach coefficient of basis_elm
                hElement = dgElement([dgBasisTuple(self.E2[i][0][hIndex],basis_elm.coord)])
            
                #print inElement, hElement, g0Element, grElement
                
                #Calculate d(G_0) * (G_1 * ... * G_j-1)
                d1Element = self.GetElementFromVector(self.d_on_vector(g0Vector,0,1),2,0) * grElement
                d1Vector = self.GetVectorFromElement(d1Element)
                
                #Calculate G_0 * d(G_1 * ... * G_j-1)
                d2Element = g0Element * self.GetElementFromVector(self.d_on_vector(grVector,0,j-1),2,j-2) 
                d2Vector = self.GetVectorFromElement(d2Element)

                # Sign comes from d(H * G) = d(H) * G + (-1)**(deg(H)) H * d(G) = (-1)**(deg(H)) H * d(G),
                # since d(h) = 0 for all h in H^*(X^n)
                sign = (-1)**(hElement.multiIndex[0])
            
                # Calculate (-1)**(deg(H)) H * d(G), where 
                # d(G) = d(G_0) * (G_1 * ... * G_j-1) + (-1)^(deg(G_0)) G_0 * d(G_1 * ... * G_j-1)
                #      = d1Element - d2Element
                outElement =  (sign * hElement) * self.GetElementFromVector(d1Vector - d2Vector, 2, j-1) ;
                outVector = self.GetVectorFromElement(outElement)

                vectorCombination.append(outVector)
            
            #Sum vectors and return
            return reduce(lambda x,y : x+y, vectorCombination)

    def d(self,i,j):
        if i+2 > 2*N or j-1 < 0:
            return linear_transformation(self.E2Vector[i][j], VectorSpace(QQ,0), lambda elm : self.d_on_vector(elm,i,j) );
        if i < 0 or j > N-1:
            return linear_transformation(VectorSpace(QQ,0), self.E2Vector[i+2][j-1], lambda elm : self.d_on_vector(elm,i,j));
        else:
            return linear_transformation(self.E2Vector[i][j], self.E2Vector[i+2][j-1], lambda elm : self.d_on_vector(elm,i,j));


        
testObject = CohomConfSpaceComplexCurve(N,G);
# testBasisElm = dgBasisTuple(([(1, 1), (0, 0), (0, 0)], [(0, 2)]),5);
# testBasisElm2 = dgBasisTuple(([(0, 0), (0, 0), (1, 0)], [(0, 1)]),-3);
# testBasisElms = dgElement([testBasisElm, testBasisElm2]);
# (i,j) = testBasisElm.multiIndex
# testVect = testObject.GetVectorFromElement(testBasisElms)
# print testObject.GetElementFromVector(testVect,i,j)
#testBasisElm3 = dgBasisTuple(([(0, 0), (0, 0), (0, 0)], [(0, 1), (0, 2)]),2)
#print testObject.GetVectorFromElement(dgElement([testBasisElm3]))
#print testObject.E2;
#print testObject.d_on_vectorbasis(vector((0,)),0,0)
#print testObject.GetTensor(2,1);

# print "In: " + str(testObject.GetElementFromVector(vector((1,0)),0,2))
# vect = testObject.d_on_vector(vector((1,0)),0,2)
# elem = testObject.GetElementFromVector(vect,2,1)
# print "Out1: " + str(elem)

# vect2 = testObject.d_on_vector(vect,2,1)
# elem2 = testObject.GetElementFromVector(vect2,4,0)
# print "Out2: " + str(elem2)

# sum_vect = 0
# for indx, coord in enumerate(testObject.E2Vector[2][1].coordinates(vect)):
#     if coord != 0:
#         print coord, testObject.E2[2][1][indx]
#         vect2 = testObject.d_on_vector(coord*testObject.E2Vector[2][1].basis()[indx],2,1)
#         elem2 = testObject.GetElementFromVector(vect2,4,0)
#         sum_vect += vect2
#         print elem2
# print sum_vect
# print testObject.GetElementFromVector(sum_vect,4,0)
# print testObject.E2[0][1][testObject.E2Vector[0][1].basis().index(vector((1,0,0)))]

# elem1 = dgElement([dgBasisTuple(([(1,0),(0,0),(0,0)],1))])
# elem2 = dgElement([dgBasisTuple(([(0,0),(1,0),(0,0)],1))])
# elem3 = dgElement([dgBasisTuple(([(0,0),(0,0),(1,0)],1))])
# print elem2 * elem1
# elem4 = dgElement([dgBasisTuple(([(0,0),(1,0),(1,0)],1))])
# elem5 = dgElement([dgBasisTuple(([(0,0),(0,0),(0,0)],[(0,2)]))])
# print elem4 * elem5

# #(i,j) = (2,1)
# #print testObject.d(i,j).kernel();
# #print testObject.d(i-2,j+1).image();
# #print testObject.d(i,j).kernel();
# #print testObject.d(i,j).kernel().quotient(testObject.d(i-2,j+1).image());

# testd2 = testObject.d_on_vectorbasis(testObject.d_on_vectorbasis(vector((1,0)),0,2),2,1)
# # testd2 = testObject.d_on_vectorbasis(vector((1,0)),0,2)
# print testObject.GetElementFromVector(testd2,4,0)

for i in xrange(2*N+1):
    for j in xrange(round((N/2)**2)+1):
        print i,j;
        
        # print testObject.d(i,j).kernel();
        # print testObject.d(i-2,j+1).image();
        print testObject.d(i,j).kernel().quotient(testObject.d(i-2,j+1).image()).dimension();
        # print testObject.d(i,j).kernel().dimension() - testObject.d(i-2,j+1).image().dimension();