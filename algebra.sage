import itertools



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

