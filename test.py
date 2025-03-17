import PyTorus as torus

A = torus.Matrix(3, 3)
a_data = [[4.0, 12.0, -16.0], [12.0, 37.0, -43.0], [-16.0, -43.0, 98.0]]
N = 3

for i in range(len(a_data)):
  A.set_row(i, torus.Vector(a_data[i]))
print(A)

result = torus.CholeskyDecomposition(A, N)
print("Result of Cholesky Decomposition: ", result)
print("Decomposed Matrix L (lower triangular):\n", A)
