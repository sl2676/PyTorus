import PyTorus as p
import math
'''
def func(x):
	return x * x -2

def test_rtsafe():
	try:
		root = p.rtsafe(func, 0, 2, 1e-5)
		ex_root = math.sqrt(2)	
		print(ex_root)
		print(abs(root-ex_root))
	except Exception as e:
		print(e)

test_rtsafe()

'''

print(p.WDabs(-1))
print(p.WDabs(-1.5))
print(p.WDabs(complex(-3, 5)))
