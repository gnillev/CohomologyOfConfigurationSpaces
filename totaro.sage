import pprint
load("cohomology.sage")
load("algebra.sage")


G = 1; # Genus of complex curve
N = 3; # Number of points in configuration space.
L = 1; # Complex dimension; No impact at the moment.


# cohomXn = CohomologyObject([1,2*G,1],"h")**N
# confAlgebraRn = ConfAlgebraObject(N);

# print cohomXn;
# print confAlgebraRn;




def relGH1(x, y):
	if y != 1:
		for gelm in y:
			for j in xrange(0,N):
				if x[j] != (0,0) and gelm[1] == j:
					return False
	return True
	
    #return not any(gelm[1] == helm or gelm[1]+N == helm for gelm in gkey for helm in hkey);

class CohomConfSpaceComplexCurve(object):
	"""docstring for ConfSpaceCohomComplexCurve"""
	def __init__(self, npoints, genus):
		self.npoints = npoints
		self.genus = genus
		self.complex_dimension = 1

		self.cohomObject = CohomologyObject([1,2*genus,1],"h")**npoints
		self.confEuclid = ConfAlgebraObject(npoints)

		self.E2 = [];
		self.E2Vector = [];

		for i,h in enumerate(self.cohomObject):
			self.E2.append([]);
			self.E2Vector.append([]);
			for j,g in enumerate(self.confEuclid):
				self.E2[i].append([]);
				for x in h:
					for y in g:
						if relGH1(x,y):
							self.E2[i][j].append((x,y));

				self.E2Vector[i].append(FreeModule(QQ,len(self.E2[i][j])));
				#print str(i) + "," + str(j) + ": " + str(self.E2[i][j]);

	def d_on_vectorbasis(self, x, i, j):


		
testObject = CohomConfSpaceComplexCurve(N,G);
print testObject.E2Vector;

class dgBasisTuple(object):
    """docstring for dgBasisTuple"""
    def __init__(self, basisTuple):
        self.basisTuple = basisTuple

        if basisTuple[0] == 1:
            i = 0;
        else:
            i = len(basisTuple[0]);

        if basisTuple[1] == 1:
            j = 0;
        else:
            j = len(basisTuple[1]);

        self.multiIndex = (i,j);