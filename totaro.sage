import pprint
load("cohomology.sage")
load("algebra.sage")


G = 2; # Genus of complex curve
N = 2; # Number of points in configuration space.
L = 1; # Complex dimension; No impact at the moment.

### Set genus for all instances of CohomologyBasisElement in the following
CohomologyBasisElement._genus = G


# cohomXn = CohomologyObject([1,2*G,1],"h")**N
# confAlgebraRn = ConfAlgebraObject(N);

# print cohomXn;
# print confAlgebraRn;


class dgBasisTuple(object):
    """docstring for dgBasisTuple"""
    def __init__(self, element, coeff = 1):
        self.element = element
        self.coeff = coeff
        self.hpart = CohomologyBasisElement(element[0])
        self.gpart = ConfAlgebraElement([ConfAlgebraBasisElement(element[1])])

        i = self.hpart.degree
        j = self.gpart.degree
        self.multi_degree = (i,j)

    def __str__(self):
        return str(self.coeff) + "*" + str(self.element)

    def __repr__(self):
        return str(self)

    def __mul__(a,b):
        # We multiply some ([(d_1,k_1),...,(d_i,k_i)],[(a_1, b_1), ..., (a_j, b_j)]) with ([(d'_1,k'_1),...,(d'_i',k'_i')],[(a'_1, b'_1), ..., (a'_j', b'_j')]).
        # First commute algebra generators from a with cohomology generators from b.
        # Calculate sign
        sign = (-1)**(a.multi_degree[1] * b.multi_degree[0])

        # Multiply cohomology generators
        hpart = a.hpart * b.hpart
        # Multiply algebra generators
        
        gpart = a.gpart * b.gpart

        # Get base coefficient
        coeff = sign * a.coeff * b.coeff * hpart.coeff 
        
        
        # Replace a[j]*G_ij with a[i]*G_ij
        helm = hpart.element
        elements = []
        if coeff != 0:
            for gsummand in gpart.elements:
                if gsummand.coeff != 0:
                    out_coeff = coeff * gsummand.coeff
                    if gsummand.element != 1:
                        #print gsummand
                        for g in reversed(gsummand.element):
                            (i,j) = (g[0], g[1])
                            if helm[j] != (0,0):
                                hsub = CohomologyBasisElement([helm[i]]) * CohomologyBasisElement([helm[j]])

                                # Factor in sign from "switching" places
                                sub_sign = 1
                                for k in range(i+1,j):
                                    sub_sign *= (-1)**(helm[j][0] * helm[k][0])

                                helm[i] = hsub.element[0]
                                helm[j] = (0,0)
                                out_coeff *= sub_sign * hsub.coeff
                    if out_coeff != 0:
                        elements.append(dgBasisTuple((helm, gsummand.element), out_coeff))

        
        return dgElement(elements)

    def __rmul__(self, other):
        return dgBasisTuple(self.element, self.coeff * other)


class dgElement(object):
    """docstring for dgElement"""
    def __init__(self, elements):
        self.elements = elements;
        
        if elements == []:
            self.multi_degree = (0,0)
        else:
            if all(x.multi_degree == elements[0].multi_degree for x in elements):
                self.multi_degree = elements[0].multi_degree
            else:
                raise ValueError("Basiselements must be of the same bi-degree.")

    def __add__(a,b):
        if a.multi_degree == b.multi_degree:
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

    def __str__(self):
        return "ConfSpaceCohomComplexCurve object with top degrees: " + str(self.top_degree);

    def __repr__(self):
        return str(self)

    def __init__(self, npoints, genus):
        self.npoints = npoints
        self.genus = genus
        self.complex_dimension = 1
        
        self.cohomObject = ComplexCurveCohomPower(genus,npoints)
        self.algebraObject = ConfAlgebraObject(npoints)
        self.top_degree = (self.cohomObject.top_degree, self.algebraObject.top_degree)

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
        (i,j) = elems.multi_degree
        out_vector = 0
        for elem in elems:
            indx = self.E2[i][j].index(elem.element)
            out_vector += elem.coeff * self.E2Vector[i][j].basis()[indx] 
        return out_vector

    def GetElementFromVector(self, vect, i, j):
        out_elements = []

        for indx, coeff in enumerate(self.E2Vector[i][j].coordinates(vect)):
            if coeff != 0:
                out_elements.append(dgBasisTuple(self.E2[i][j][indx],coeff))

        return dgElement(out_elements)

    def relGH1(self, x, y):
        if y != 1:
            for gelm in y:
                k = gelm[1]
                if x[k] != (0,0):
                    return False
        return True     

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
                vectorCombination.append(basis_elm.coeff * out_vector)

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
                hElement = dgElement([dgBasisTuple(self.E2[i][0][hIndex],basis_elm.coeff)])
            
                #print inElement, hElement, g0Element, grElement

                #Calculate d(G_0) * (G_1 * ... * G_j-1)
                d1Element = self.GetElementFromVector(self.d_on_vector(g0Vector,0,1),2,0) * grElement
                d1Vector = self.GetVectorFromElement(d1Element)
                
                #Calculate G_0 * d(G_1 * ... * G_j-1)
                d2Element = g0Element * self.GetElementFromVector(self.d_on_vector(grVector,0,j-1),2,j-2) 
                d2Vector = self.GetVectorFromElement(d2Element)

                # Sign comes from d(H * G) = d(H) * G + (-1)**(deg(H)) H * d(G) = (-1)**(deg(H)) H * d(G),
                # since d(h) = 0 for all h in H^*(X^n)
                sign = (-1)**(hElement.multi_degree[0])
                
                # Calculate (-1)**(deg(H)) H * d(G), where 
                # d(G) = d(G_0) * (G_1 * ... * G_j-1) + (-1)^(deg(G_0)) G_0 * d(G_1 * ... * G_j-1)
                #      = d1Element - d2Element
                outElement =  (sign * hElement) * self.GetElementFromVector(d1Vector - d2Vector, 2, j-1) ;
                outVector = self.GetVectorFromElement(outElement)

                vectorCombination.append(outVector)
            
            #Sum vectors and return
            return reduce(lambda x,y : x+y, vectorCombination)

    def d(self,i,j):

        if i+2 > self.top_degree[0] or j-1 < 0:
            return linear_transformation(self.E2Vector[i][j], VectorSpace(QQ,0), lambda elm : self.d_on_vector(elm,i,j) );
        if i < 0 or j > self.top_degree[1]:
            return linear_transformation(VectorSpace(QQ,0), self.E2Vector[i+2][j-1], lambda elm : self.d_on_vector(elm,i,j));
        else:
            return linear_transformation(self.E2Vector[i][j], self.E2Vector[i+2][j-1], lambda elm : self.d_on_vector(elm,i,j));

    def _FormCohomology(self):
        self._ECohom = []
        self._BettiNumbers = []
        for i in xrange(0,self.top_degree[0]+1):
            self._ECohom.append([])
            self._BettiNumbers.append([])
            for j in xrange(0,self.top_degree[1]+1):
                #print i,j;
                self._ECohom[i].append(testObject.d(i,j).kernel().quotient(testObject.d(i-2,j+1).image()))
                self._BettiNumbers[i].append(self._ECohom[i][j].dimension())

    def GetBettiNumbers(self):
        if not hasattr(self, "_BettiNumbers"):
            self._FormCohomology()
        return self._BettiNumbers

    def GetBettiNumber(self,i,j):
        return self.GetBettiNumbers()[i][j]

    def GetECohom(self,i,j):
        if not hasattr(self, "_ECohom"):
            self._FormCohomology()
        return self._ECohom[i][j]

    def PrintBettiNumbers(self):
        out_str = "Betti numbers \nN=" + str(self.npoints) + ", G=" + str(self.genus) + ": \n" 
        #print self.GetBettiNumbers()
        for j in xrange(self.top_degree[1],-1,-1):
            for i in xrange(0,len(self.GetBettiNumbers())):
                out_str += str(i) + "," + str(j) + ": " + str(self.GetBettiNumbers()[i][j]) + "\t"
            out_str += "\n"
        return out_str

    def PrintE2Dimensions(self):
        out_str = "E2 Dimensions \nN=" + str(self.npoints) + ", G=" + str(self.genus) + ": \n" 
        #print self.GetBettiNumbers()
        for j in xrange(self.top_degree[1],-1,-1):
            for i in xrange(0,self.top_degree[0]+1):
                out_str += str(i) + "," + str(j) + ": " + str(len(self.E2[i][j])) + "\t"
            out_str += "\n"
        return out_str

testObject = CohomConfSpaceComplexCurve(N,G)
print testObject
print testObject.cohomObject
print testObject.algebraObject
print testObject.PrintE2Dimensions()
print testObject.PrintBettiNumbers()

# print "In: " + str(testObject.GetElementFromVector(vector((1,0)),0,2))
# vect = testObject.d_on_vector(vector((1,0)),0,2)
# elem = testObject.GetElementFromVector(vect,2,1)
# print "Out1: " + str(elem)

# vect2 = testObject.d_on_vector(vect,2,1)
# elem2 = testObject.GetElementFromVector(vect2,4,0)
# print "Out2: " + str(elem2)

# sum_vect = 0
# for indx, coeff in enumerate(testObject.E2Vector[2][1].coordinates(vect)):
#     if coeff != 0:
#         print coeff, testObject.E2[2][1][indx]
#         vect2 = testObject.d_on_vector(coeff*testObject.E2Vector[2][1].basis()[indx],2,1)
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


#print testObject.GetTensor(0,3)

# for i in xrange(2*N+1):
#     for j in xrange(round((N/2)**2)+1):
#         print i,j;
        
#         # print testObject.d(i,j).kernel();
#         # print testObject.d(i-2,j+1).image();
#         print testObject.d(i,j).kernel().quotient(testObject.d(i-2,j+1).image()).dimension();
#         # print testObject.d(i,j).kernel().dimension() - testObject.d(i-2,j+1).image().dimension();
