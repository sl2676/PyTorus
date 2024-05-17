import PyTorus as p
import math
import numpy as np

a_test = p.Matrix(3, 3)
maddie_mat = p.Matrix(3, 3)
vec = p.Vector(3, 2, 3)
maddie_mat.set_row(0, vec)
vec = p.Vector(0, 5, 4)
maddie_mat.set_row(1, vec)
vec = p.Vector(0, 4, 6)
maddie_mat.set_row(2, vec)
vex = p.Vector(1, 0, 0)
a_test.set_row(0, vex)
vex = p.Vector(0, 1, 0)
a_test.set_row(1, vex)
vex = p.Vector(0, 0, 1)
a_test.set_row(2, vex)
print(maddie_mat)
print(a_test)

vec = p.Vector(0, 0, 0.0)

print(vec)
print("--------------------------")
print("first gauss\n")

val = p.GaussJordan(maddie_mat, 3, a_test, 3)
print("a_test")
print(a_test)
print("maddie_mat")
print(maddie_mat)
print(val)

print("--------------------------")
print("sec gauss")

a_test = p.Matrix(3, 3)
maddie_mat = p.Matrix(3, 3)
vec = p.Vector(3, 2, 3)
maddie_mat.set_row(0, vec)
vec = p.Vector(0, 5, 4)
maddie_mat.set_row(1, vec)
vec = p.Vector(0.0, 4, 6)
maddie_mat.set_row(2, vec)
vex = p.Vector(1, 0, 0)
a_test.set_row(0, vex)
vex = p.Vector(0, 1, 0)
a_test.set_row(1, vex)
vex = p.Vector(0, 0, 1)
a_test.set_row(2, vex)
print(maddie_mat)
print(a_test)

vec = p.Vector(1, 1, 1)
print(vec)
val = p.GaussJordan(maddie_mat, 3, vec)
print("a_test")
print(a_test)

print("maddie_mat")
print(maddie_mat)
print("vec")
print(vec)

print("---------------------")
vec = p.Vector(0, 0, 0)
a_test = p.Matrix(3, 3)
maddie_mat = p.Matrix(3, 3)
vec = p.Vector(3, 2, 3)
maddie_mat.set_row(0, vec)
vec = p.Vector(0, 5, 4)
maddie_mat.set_row(1, vec)
vec = p.Vector(0, 4, 6)
maddie_mat.set_row(2, vec)
vex = p.Vector(1, 0, 0)
a_test.set_row(0, vex)
vex = p.Vector(0, 1, 0)
a_test.set_row(1, vex)
vex = p.Vector(0, 0, 1)
a_test.set_row(2, vex)

print(maddie_mat)

val = p.GaussBack(maddie_mat, 3, vec)

print(maddie_mat)
print(vec)
print(val)
print()
print("---------------------")
print("qbulrisch")

def eq_func(x):
	return math.exp(-x * x)
a = 0.0
b = 1.0
eps = 1e-6
err = 0.0
result = p.qbulir(eq_func, a, b, eps, err)

print(f"Approximated integral: {result}")
print(f"Estimated error: {err}\n")

def sinFunc(x):
    return math.sin(x)

a = 0.0
b = math.pi
eps = 1e-6
err = 0.0

result = p.qbulir(sinFunc, a, b, eps, err)

print(f"Approximated integral of sin(x) from {a} to {b}: {result}")
print(f"Estimated error: {err}\n")


def expFunc(x):
    return math.exp(x)

a = 0.0  
b = 1.0  
eps = 1e-12  
err = 0.0

result = p.qbulir(expFunc, a, b, eps, err)

print(f"Approximated integral of exp(x) from {a} to {b}: {result}")
print(f"Estimated error: {err}\n")

print("---------------------")
print("GaussLegendre")
n = 5
x = p.Vector(0, 0, 0, 0, 0)
w = p.Vector(0, 0, 0, 0, 0)

val = p.GaussLegendre(x, w, n)

print("Roots (x):", x)
print("Weights (w):", w)
print()

print("---------------------")
print("dLegendrePeven")

x_value = 0.5
np_value = 5

pec = p.Vector(0, 0, 0, 0, 0)

p.LegendrePeven(pec, x_value, np_value)

print("Legendre Polynomials (even):", pec)
print()

print("---------------------")
print("dLegendrePeven")

x_value = 0.5
np_value = 5

pec = p.Vector(0, 0, 0, 0, 0)
d = p.Vector(0, 0, 0, 0, 0)

p.dLegendrePeven(pec, d, x_value, np_value)

print("Legendre Polynomials (even):", pec)
print("Derivatives of Legendre Polynomials (even):", d)
print()

print("---------------------")
print("zbrent")

def sample_function(x):
    return np.cos(x) - x

def sample_function_wrapper(x):
    return sample_function(x)

x1 = 0.0
x2 = 1.0
tol = 1e-6

root = p.zbrent(sample_function_wrapper, x1, x2, tol)

print(f"Root of the function in the interval [{x1}, {x2}] with tolerance {tol}: {root}")
print()

print("---------------------")
print("heap_index")

A = p.Vector(3.2, 1.5, 4.8, 2.9, 0.6)
indx = p.Vector(0, 0, 0, 0, 0)

p.heap_index(A, 5, indx)

print("Original array A:", A)
print("Index array indx sorted by values in A:", indx)
sorted_A = [A[indx[i]] for i in range(5)]
print("Sorted array A based on indx:", sorted_A)
print()

def my_func(x):
    return x ** 2

n = 10
indx = p.Vector(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
p.heap_index(my_func, n, indx)

print("Sorted indices according to my_func:")
print(indx)
print()

print("---------------------")
print("qsplin")

x = p.Vector(0.0, 1.0, 2.0, 3.0, 4.0)
y = p.Vector(1.0, 2.0, 0.0, 3.0, 4.0)
y2 = p.Vector(0.0, -0.5, 0.0, 0.5, 0.0)

n = 5
al = 1.0  
x1 = 0.5 
x2 = 3.5  

result = p.qsplin(x, y, y2, n, al, x1, x2)

print("Result of qsplin integration:", result)
print()

print("---------------------")
print("CholeskyDecomposition")
a = p.Matrix(3, 3)
vec = p.Vector(4.0, 12.0, -16.0)
a.set_row(0, vec)
vec = p.Vector(12.0, 37.0, -43.0)
a.set_row(1, vec)
vec = p.Vector(-16.0, -43.0, 98.0)
a.set_row(2, vec)

n = 3
print(a)
result = p.CholeskyDecomposition(a, n)
print("Result of Cholesky Decomposition:", result)
print("Decomposed Matrix L (lower triangular):")
print(a)
print()

print("---------------------")
print("CholeskySolution")
a = p.Matrix(3, 3)
vec = p.Vector(4.0, 12.0, -16.0)
a.set_row(0, vec)
vec = p.Vector(12.0, 37.0, -43.0)
a.set_row(1, vec)
vec = p.Vector(-16.0, -43.0, 98.0)
a.set_row(2, vec)
b = p.Vector(1.0, 1.0, 1.0)
print(a)

n = 3
p.CholeskyDecomposition(a, n)

p.CholeskySolution(a, n, b)

print("Solution vector b after Cholesky Decomposition and solving:")
print(b)

print("---------------------")
print("CholeskyInvertL")

a = p.Matrix(3, 3)
vec = p.Vector(4.0, 12.0, -16.0)
a.set_row(0, vec)
vec = p.Vector(12.0, 37.0, -43.0)
a.set_row(1, vec)
vec = p.Vector(-16.0, -43.0, 98.0)
a.set_row(2, vec)
b = p.Vector(1.0, 1.0, 1.0)
print(a)
n = 3

p.CholeskyDecomposition(a, n)
p.CholeskyInvertL(a, n)

print("CholeskyInvertL: ")

print(a)

print("---------------------")
print("CholeskyInvertF")

a = p.Matrix(3, 3)
vec = p.Vector(4.0, 12.0, -16.0)
a.set_row(0, vec)
vec = p.Vector(12.0, 37.0, -43.0)
a.set_row(1, vec)
vec = p.Vector(-16.0, -43.0, 98.0)
a.set_row(2, vec)
b = p.Vector(1.0, 1.0, 1.0)
print(a)

n = 3
p.CholeskyDecomposition(a, n)

p.CholeskyInvertL(a, n)

p.CholeskyInvertF(a, n)
print("Inverse of the original matrix A:")
print(a)
print()
'''
print("---------------------")
print("LUDecomposition")
a = p.Matrix(3, 3)
vec = p.Vector(4.0, 12.0, -16.0)
a.set_row(0, vec)
vec = p.Vector(12.0, 37.0, -43.0)
a.set_row(1, vec)
vec = p.Vector(-16.0, -43.0, 98.0)
a.set_row(2, vec)
b = p.Vector(1.0, 1.0, 1.0)
print(a)

indx = p.Vector(0, 0, 0)
n = 3
d = 1
p.LUDecomposition(a, n, indx, d)
#p.LUSolution(a, n, indx, b)
print("LU Decomposed Matrix A:")
print()
print(a)

print("Index Vector indx:")
print(indx)

print("Parity d:")
print(d)
print()
'''
print("---------------------")
print("LUDet3")
a = p.Matrix(3, 3)
vec = p.Vector(4.0, 12.0, -16.0)
a.set_row(0, vec)
vec = p.Vector(12.0, 37.0, -43.0)
a.set_row(1, vec)
vec = p.Vector(-16.0, -43.0, 98.0)
a.set_row(2, vec)
b = p.Vector(1.0, 1.0, 1.0)
print(a)

det_A = p.LUDet3(a)
print("Determinant of the matrix A:")
print(det_A)
print()

b = p.Matrix(3, 3)
vec = p.Vector(6, 1, 1)
b.set_row(0, vec)
vec = p.Vector(4, -2, 5)
b.set_row(1, vec)
vec = p.Vector(2, 8, 7)
b.set_row(2, vec)
print(b)

det_B = p.LUDet3(b)
print(det_B)

c = p.Matrix(3, 3)
vec = p.Vector(1, 1, 1)
for i in range(3):
	c.set_row(i, vec)

print(c)

det_C = p.LUDet3(c)
print(det_C)

print("---------------------")
print("LUDecomposition")

a = p.Matrix(3, 3)
vec = p.Vector(4, 12, -16)
a.set_row(0, vec)
vec = p.Vector(12, 37, -43)
a.set_row(1, vec)
vec = p.Vector(-16, -43, 98)
a.set_row(2, vec)
print(a)
b = p.Vector(1, 1, 1)
indx = p.Vector(0, 0, 0)

n = 3
d = 1
g = p.LUDecomposition(a, n, indx, d)
#p.LUSolution(a, n, indx, b)
print(a)
print(g)
print("Vector indx")
print(indx)

print("Parity")
print(d)

print("Solution Vector")
print(b)


