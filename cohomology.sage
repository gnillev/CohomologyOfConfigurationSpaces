from collections import defaultdict

class CohomologyBasisElement(object):
    """docstring for CohomologyBasisElement"""

    #Static variable, overwritten in use
    _genus = None;

    def __init__(self, element, coeff = 1):
        self.element = element
        self.degree = reduce(lambda x,y : x+y, map(lambda x: x[0], element))
        self.coeff = coeff

    def __str__(self):
        return str(self.coeff) + "*" + str(self.element)

    def __repr__(self):
        return str(self)

    def __mul__(a,b):
        if len(a.element) != len(b.element):
            raise ValueError("Different number of tensor coeffinates.")
        else:
            length = len(a.element)

        if length == 1:
            # Implement relations between generators of H^*(X). 
            # Assuming H^*(X) is of dimensions [1, 2*genus, 1]:

            if a.degree == 0:
                return CohomologyBasisElement(b.element, a.coeff * b.coeff)
            elif b.degree == 0:
                return CohomologyBasisElement(a.element, a.coeff * b.coeff)
            elif a.degree == 1 and b.degree == 1:
                if  a[0][1] == (b[0][1] - a._genus):
                    return CohomologyBasisElement([(2,0)], a.coeff * b.coeff)
                if a[0][1] == (b[0][1] + a._genus):
                    return CohomologyBasisElement([(2,0)], -1*a.coeff * b.coeff)
                else:
                    return CohomologyBasisElement([(2,0)],0)
            else:
                return CohomologyBasisElement([(0,0)],0)
        elif length > 1:

            # Implement tensor product multiplication recursively
            a0 = CohomologyBasisElement([a[0]])
            a1 = CohomologyBasisElement(a[1:])
            b0 = CohomologyBasisElement([b[0]])
            b1 = CohomologyBasisElement(b[1:])

            result_element = (a0 * b0).element + (a1 * b1).element
            result_coeff = a.coeff * b.coeff * (-1)**(a1.degree * b0.degree) * (a0 * b0).coeff * (a1 * b1).coeff

            return CohomologyBasisElement(result_element, result_coeff)
        else:
            raise ValueError("Cannot multiply 'empty' elements")

    def __rmul__(self, other):
        return CohomologyBasisElement(self.element, self.coeff * other)

    def __getitem__(self, item):
        return self.element[item];

class CohomologyElement(object):
    """docstring for CohomologyElement"""
    def __init__(self, elements):
        if all(x.degree == elements[0].degree for x in elements):
            self.degree = elements[0].degree
        else:
            raise ValueError("Basiselements must be of the same degree.")

        self.elements = []

        tally = defaultdict(list)
        for i,x in enumerate(elements):
            tally[tuple(x.element)].append(i)
        for elm, positions in tally.items():
            coeff = 0
            for pos in positions:
                coeff += elements[pos].coeff
            self.elements.append(CohomologyBasisElement(list(elm),coeff))

    def __str__(self):
        return str(self.elements)

    def __repr__(self):
        return str(self)

    def __mul__(a,b):
        result_elements = [];
        for aelm in a.elements:
            for belm in b.elements:
                result_elements.append(aelm * belm)
        return CohomologyElement(result_elements)


    def __rmul__(self, other):
        return CohomologyElement(map(lambda elm : elm * other, self.elements))

    def __getitem__(self, item):
        return self.elements[item];

#Tests:
#CohomologyBasisElement._genus = 1
elem = CohomologyBasisElement([(1, 0), (1, 1), (0, 0)]);
elem2 = CohomologyBasisElement([(0, 0), (0, 0), (1, 0)]);
#print elem * elem2;
#elemList = CohomologyElement([elem, elem2])
# print elemList;



class CohomologyObject(object):
    """Implements a cohomology ring as a collection of vector spaces.
    Assumes all H^i(X;R) are free over R.
    Initial input: (List of dimensions of H^*(X;R), base ring (R) = QQ)
    Methods: 
    self[i] returns the (Sage) FreeModule object of degree i.
    
    self**n returns the CohomologyObject of H^*(X^n;R), calculated by means of the Kunneth formula, assuming no torsion.
    H^k(X^n;R) = H^(k_1) (X; R) @ ... @ H^(k_n) (X; R)
    where SUM(k_1, ..., k_n) = k.

    """
    def __init__(self, dimensions, power=1, base_ring = QQ):
        self.H = []
        self.dimensions = dimensions
        self.top_degree = len(self.dimensions)-1
        for indx, dimension in enumerate(dimensions):
            self.H.append(FreeModule(base_ring, dimension))

        if power > 1:
            self.base_cohomology = self
            self.power = power
            self._cohomology = self.ExplicitPower(self.power)
            self.dimensions = self._cohomology.dimensions
            self.top_degree = self._cohomology.top_degree
            self.H = self._cohomology.H
            self._BasisRepresentation = self._cohomology.BasisRepresentation()
            self._TensorRepresentation = self._cohomology.TensorRepresentation()


    def __str__(self):
        strReturn =  "CohomologyObject of dimensions: "
        for i,d in enumerate(self.dimensions):
            strReturn += "H^" + str(i) + ": " + str(d) + ", "

        return strReturn

    def __getitem__(self, degree):
        return self.TensorRepresentation()[degree]

    def __pow__(self, value):
        # top_degree = len(self.dimensions)-1;
        # return_dimensions = []
        # for i in xrange(0,top_degree*value+1):
        #   outdim = 0
        #   for partition in IntegerVectors(i, length=value, max_part=len(self.dimensions)-1).list():
        #       indim = 1
        #       for dim in partition:
        #           indim *= self.dimensions[dim]
        #       outdim += indim
        #   return_dimensions.append(outdim)
        return self.ExplicitPower(value)

    def BasisRepresentation(self):
        if not hasattr(self, "_BasisRepresentation"):
            basisList = [];
            for i,d in enumerate(self.dimensions):
                basisList.append([]);
                for j in xrange(d):
                    basisList[i].append( (i,j) )
            self._BasisRepresentation = basisList
        return self._BasisRepresentation

    def TensorRepresentation(self):
        if not hasattr(self,"_TensorRepresentation"):
            return self.BasisRepresentation()
        else:
            return self._TensorRepresentation

    def ExplicitPower(self, value):
        result_tensor = []
        result_dimensions = []
        for k in xrange(0, self.top_degree*value + 1):
            result_tensor.append([])

            # Iterate through a list of possible non-negative integer partitions of k, having length = value, and highest integer summand equal to the top dimension.
            for partition in IntegerVectors(k, length=value, max_part=self.top_degree).list():
                basisList = []

                # Store the possible basis elements from the degrees taken from the integer partition
                for i in partition:
                    basisList.append(self.BasisRepresentation()[i])

                # A little cryptic, but the basic idea is to form all possible products of elements from the lists of basis elements, 
                # in degrees determined by the integer partition:

                # Initialize the list of elements to be added to degree k for this partition
                lstVals = [[]]
                # Convert the lists of basis elements to immutable tuples
                pools = map(tuple, basisList)
                for pool in pools:
                    # For each tuple, add the "products" of all the existing elements in lstVals and all the elements from the tuple,
                    # and store in lstVals (to be used in the next iteration)
                    lstVals = [x+[y] for x in lstVals for y in pool]

                # Add the list of "tensor products" to the degree k elements
                result_tensor[k]+=lstVals
            # Add the number of distinct tensor products of basis elements as dimension of degree k
            result_dimensions.append(len(result_tensor[k]))
        result = CohomologyObject(result_dimensions)
        result._TensorRepresentation = result_tensor
        result.power = value
        result.base_cohomology = self
        return result

    def GetVector(self, cohomElement):
        tensorBasis = self.TensorRepresentation()[cohomElement.degree]
        position = tensorBasis.index(cohomElement.element)
        return cohomElement.coeff * self.H[cohomElement.degree].basis()[position]


#Tests:
#testCohom = CohomologyObject([1,2,1],"h",3)
#print testCohom.GetVector(CohomologyElement([(1,1),(1,0),(0,0)]))
#print testCohom.TensorRepresentation();
#print (testCohom**2).TensorRepresentation();
#print testCohom.BasisRepresentation()
#print testCohom.ExplicitPower(2)[2]
#print IntegerVectors(4, 4, max_part = 2).list();
#print IntegerVectors(4, length=2).list();

    

class ComplexCurveCohomPower(CohomologyObject):
    def __init__(self, genus,power = 1, base_ring = QQ):
        self.genus = genus
        super(ComplexCurveCohomPower, self).__init__([1,2*self.genus,1],power, base_ring)
        #self.cohomObject = CohomologyObject([1,2*self.genus,1],"h",self.power)

    def Diagonal(self,i,j):
        """ p^*_{i,j}([D]) where [D] \in H^n(X^2;R) is the class of the diagonal for an oriented manifold of genus g """
        
        diagonal_sum = [];
        diagonal_sum += [CohomologyBasisElement([(0,0),(2,0)])]

        for k in xrange(0,self.genus):
            diagonal_sum += [CohomologyBasisElement([(1,k),(1,k+self.genus)],-1), CohomologyBasisElement([(1,k+self.genus),(1,k)],1)]

        diagonal_sum += [CohomologyBasisElement([(2,0),(0,0)])]

        if self.power > 1:
            
            out_sum = []
            for x in diagonal_sum:
                out_elm = []
                for k in xrange(0,self.power):
                    if k == i:
                        out_elm.append(x[0])
                    elif k == j:
                        out_elm.append(x[1])
                    else:
                        out_elm.append((0,0))
                out_sum.append(CohomologyBasisElement(out_elm, x.coeff))

            return out_sum
        else:
            return diagonal_sum


#Tests:
# testComplexCurveCohom = ComplexCurveCohomPower(1,3)
# # print testComplexCurveCohom.TensorRepresentation()
# #print map(lambda c : c.degree, testComplexCurveCohom.Diagonal(0,2))
# elem2 = CohomologyBasisElement([(1,0),(0,0),(0,0)],1)
# elem1 = CohomologyBasisElement([(0,0),(1,0),(0,0)],1)
# print elem1 * elem2

