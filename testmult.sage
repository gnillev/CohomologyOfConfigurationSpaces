# Test of multiplication defined outside scope of class:

def mult_function1(a,b):
	return TestClass(str(a.arg)+ " mult " + str(b.arg))

class TestClass(object):
	"""docstring for TestClass"""
	mult_function = None;

	def __init__(self, arg):
		super(TestClass, self).__init__()
		self.arg = arg
		print arg

	def DefineMult(self,mult_function):
		self.mult_function = mult_function

	def __mul__(a, b):
		return a.mult_function(b)



def mult_function3(a,b):
	return TestClass(str(a.arg)+ " mult3 " + str(b.arg))

# TestClass.mult_function = mult_function3

# testObject = TestClass("Obj1")
# testObject2 = TestClass("Obj2")
# result = testObject * testObject2

print latex(IntegerVectors(2, length=3, max_part=2).list());
