import itertools

class ConfAlgebraBasisElement(object):
    """docstring for ConfAlgebraBasisElement"""

    def __init__(self, element, coeff = 1):
        self.element = element

        if element == 1:
            self.degree = 0
        else:
            self.degree = len(self.element)
        self.coeff = coeff
        

    def __str__(self):
        return str(self.coeff) + "*" + str(self.element)

    def __repr__(self):
        return str(self)

    def __mul__(a,b):
        if a.element == 1:
            res = ConfAlgebraBasisElement(b.element, a.coeff * b.coeff)
            return res
        elif b.element == 1:
            return ConfAlgebraBasisElement(a.element, a.coeff * b.coeff)
        else:
            return ConfAlgebraBasisElement(a.element + b.element, a.coeff * b.coeff)

    def __rmul__(self, other):
        return ConfAlgebraBasisElement(self.element, self.coeff * other)

    def __getitem__(self, item):
        if self.element = 1:
            return 1
        else:
            return self.element[item];


class ConfAlgebraElement(object):
    """docstring for ConfAlgebraElement"""
    def __init__(self, elements):

        if elements != []:

            if all(x.degree == elements[0].degree for x in elements):
                self.degree = elements[0].degree
            else:
                raise ValueError("Basiselements must be of the same degree.")
            
            expand_elements = []

            for elm in elements:
                expand_elements += self.expandRelations(elm.coeff, elm.element)

            self.elements = []
            tally = defaultdict(list)
            for i,x in enumerate(expand_elements):
                
                if x.element == 1:
                    tally[x.element].append(i)
                else:
                    tally[tuple(x.element)].append(i)

            for elm, positions in tally.items():
                coeff = 0
                for pos in positions:
                    coeff += expand_elements[pos].coeff
                if elm == 1:
                    self.elements.append(ConfAlgebraBasisElement(elm,coeff))
                else:
                    self.elements.append(ConfAlgebraBasisElement(list(elm),coeff))
        else:
            self.degree = 0
            self.elements = []

        

    def __str__(self):
        return str(self.elements)

    def __repr__(self):
        return str(self)

    def __mul__(a,b):
        result_elements = [];

        for aelm in a.elements:
            for belm in b.elements:
                result_elements.append(aelm * belm)

        return ConfAlgebraElement(result_elements)

    def __rmul__(self, other):
        return ConfAlgebraElement(map(lambda elm : elm * other, self.elements))

    def __getitem__(self, item):
        return self.elements[item]

    def reorderList(self,element):
        # Handle zero-case
        if element == 1:
            return (1,1);

        # Check for duplicates, i.e. (G_{i,j})^2 factors
        if not len(set(element)) == len(element):
            return (0,0);

        elementSorted = sorted(element);
        permutationSorted = Permutation(map(lambda x : element.index(x)+1, elementSorted));
        return (permutationSorted.sign(), elementSorted);

    def expandRelations(self, sign, element):
        reorder = self.reorderList(element);
        (sign, element) = (sign * reorder[0], reorder[1])

        if element == 1:
            return [ConfAlgebraBasisElement(element,sign)]
        if element == 0:
            return []

        for indx, key in enumerate(element):
            triangleMatch = [ x for x in element[indx+1:] if x[0]==key[1] ]
            if len(triangleMatch) > 0:
                key2 = triangleMatch[0]
                restList = list(element)
                restSign = (-1)**indx
                restList.remove(key)
                restSign = (-1)**restList.index(key2) * restSign
                restList.remove(key2)
                out_elements = []
                out_elements += self.expandRelations(sign * restSign, [(key[0],key[1]),(key[0], key2[1])] + restList)
                out_elements += self.expandRelations(sign * restSign, [(key[0],key2[1]),(key2[0], key2[1])] + restList)
                return out_elements
        return [ConfAlgebraBasisElement(element,sign)]

#testElm = ConfAlgebraBasisElement([(0,1),(0,2),(1,3)]);
#testElm2 = ConfAlgebraBasisElement([(1,2),(1,4),(3,4)]);
#print testElm2;
#testElms = ConfAlgebraElement([testElm]);
#testElms2 = ConfAlgebraElement([testElm2]);
#print (testElms2*testElms)[0].coeff;
#testElms3 = ConfAlgebraElement([ConfAlgebraBasisElement([(0, 2),(0, 3),(1, 2), (1, 3), (1, 4), (3, 4)])]);
#print testElms3


class ConfAlgebraObject(object):
    """docstring for ConfAlgebraObject"""
    def __init__(self, npoints, base_ring = QQ):
        self.npoints = npoints
        self.base_ring = base_ring

        self.rels = [self._relConf1, self._relConf2]
        self.G = []
        self.G.append([1])
        self.G.append([])
        for i in xrange(0, self.npoints-1):
            for j in xrange(i+1, self.npoints):
                self.G[1] += [[(i,j)]]

        k = 2
        while len(self.G[k-1])>0:
            self.G.append([])
            for x in self.G[1]:
                for y in self.G[k-1]:
                    if all(f(x,y) for f in self.rels):
                        self.G[k].append(x+y)
            k += 1
        del self.G[-1]

        self.GVector = []
        for k in xrange(0,len(self.G)):
            self.GVector.append(FreeModule(self.base_ring,len(self.G[k])))

        # Public attributes
        self.top_degree = len(self.G)-1

    def _relConf1(self,x,y):
        return x[0] < y[0]

    def _relConf2(self,x,y):
        return all(z[0] != x[0][1] for z in y)

    def GetBasis(self,degree):
        # if degree >= len(self.G):
        #     return []
        # else:
        return self.G[degree]

    def GetVectorSpace(self,degree):
        # if degree >= len(self.G):
        #     return FreeModule(self.base_ring,0)
        # else:
        return self.GVector[degree]

    def __str__(self):
        strReturn =  "ConfAlgebraObject with " + str(self.npoints) + " points: "
        for i,d in enumerate(self.G):
            strReturn += "G(" + str(i) + "): " + str(len(d)) + ", "
        return strReturn

    def __getitem__(self, degree):
        return self.GetBasis(degree)


# Tests:
# testConfAlgebraObject = ConfAlgebraObject(3);
# print testConfAlgebraObject.G;

