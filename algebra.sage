import itertools

class AlgebraBasisElement(object):
    """docstring for AlgebraBasisElement"""

    def __init__(self, element, coord = 1):
        self.element = element

        if element == 1:
            self.degree = 0
        else:
            self.degree = len(self.element)
        self.coord = coord
        

    def __str__(self):
        return str(self.coord) + "*" + str(self.element)

    def __repr__(self):
        return str(self)

    def __mul__(a,b):
        if a.element == 1:
            res = AlgebraBasisElement(b.element, a.coord * b.coord)
            return res
        elif b.element == 1:
            return AlgebraBasisElement(a.element, a.coord * b.coord)
        else:
            return AlgebraBasisElement(a.element + b.element, a.coord * b.coord)

    def __rmul__(self, other):
        return AlgebraBasisElement(self.element, self.coord * other)

    def __getitem__(self, item):
        return self.element[item];


class AlgebraElement(object):
    """docstring for AlgebraElement"""
    def __init__(self, elements):

        if elements != []:
            
            if all(x.degree == elements[0].degree for x in elements):
                self.degree = elements[0].degree
            else:
                raise ValueError("Basiselements must be of the same degree.")
            
            expand_elements = []

            for elm in elements:
                expand_elements += self.expandRelations(elm.coord, elm.element)

            self.elements = []
            tally = defaultdict(list)
            for i,x in enumerate(expand_elements):
                
                if x.element == 1:
                    tally[x.element].append(i)
                else:
                    tally[tuple(x.element)].append(i)

            for elm, positions in tally.items():
                coord = 0
                for pos in positions:
                    coord += expand_elements[pos].coord
                if elm == 1:
                    self.elements.append(AlgebraBasisElement(elm,coord))
                else:
                    self.elements.append(AlgebraBasisElement(list(elm),coord))
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

        return AlgebraElement(result_elements)

    def __rmul__(self, other):
        return AlgebraElement(map(lambda elm : elm * other, self.elements))

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
            return [AlgebraBasisElement(element,sign)]
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
        return [AlgebraBasisElement(element,sign)]

#testElm = AlgebraBasisElement([(0,1),(0,2),(1,3)]);
#testElm2 = AlgebraBasisElement([(1,2),(1,4),(3,4)]);
#print testElm2;
#testElms = AlgebraElement([testElm]);
#testElms2 = AlgebraElement([testElm2]);
#print (testElms2*testElms)[0].coord;
#testElms3 = AlgebraElement([AlgebraBasisElement([(0, 2),(0, 3),(1, 2), (1, 3), (1, 4), (3, 4)])]);
#print testElms3


class ConfAlgebraObject(object):
    """docstring for ConfAlgebraObject"""
    def __init__(self, dimension):
        self.dimension = dimension
        self.rels = [self.relConf1, self.relConf2]
        self.G = []
        self.G.append([1])
        self.G.append([])
        for i in xrange(0, self.dimension-1):
            for j in xrange(i+1, self.dimension):
                #self.G[1] += ['g{'+str(i)+','+str(j)+'}']
                self.G[1] += [[(i,j)]]

        for i in xrange(2,round((self.dimension/2)**2)+1):
            self.G.append([])
            for x in self.G[1]:
                for y in self.G[i-1]:
                        if all(f(x,y) for f in self.rels):
                            self.G[i].append(x+y)

    def relConf1(self,x,y):
        return x[0] < y[0]

    def relConf2(self,x,y):
        return all(z[0] != x[0][1] for z in y)

    def __str__(self):
        strReturn =  "ConfAlgebraObject with " + str(self.dimension) + " points: "
        for i,d in enumerate(self.G):
            strReturn += "G(" + str(i) + "): " + str(len(d)) + ", "
        return strReturn

    def __getitem__(self, degree):
        return self.G[degree]


# Tests:
#testConfAlgebraObject = ConfAlgebraObject(3);
#print testConfAlgebraObject.G;

