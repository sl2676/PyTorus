import PyTorus as p
print("\\***************\\")
print("|TESTING MATRICES|")
print(dir(p.Matrix))
print("\\****************\\")

vi1 = p.Vector(1,2,3)
vi2 = p.Vector(3, 2, 1)

vf1 = p.Vector(1.5, 2.5, 3.5)
vf2 = p.Vector(3.5, 2.5, 1.5)

vc1 = p.Vector(complex(1, 5), complex(2, 5), complex(3, 5))
vc2 = p.Vector(complex(3, 5), complex(2, 5), complex(1, 5))

vex = p.Vector(10.0, 10, 10)


print("\ninteger vectors")
print(vi1)
print(vi2)

print("\nfloat vectors")
print(vf1)
print(vf2)

print("\ncomplex vectors")
print(vc1)
print(vc2)

matrix_int = p.Matrix(3, 3)
matrix_float = p.Matrix(3, 3)
matrix_complex = p.Matrix(3, 3)
test_matrix = p.Matrix(3, 3)

print("\nint matrix")
matrix_int.set_column(0, vi1)
matrix_int.set_column(1, vf2)
matrix_int.set_row(2, vex)
tmatrix = matrix_int * matrix_int
print(matrix_int)
print(tmatrix)

print("\ntest matrix")
print(test_matrix)
test_matrix = matrix_int
ex = test_matrix[1]
print(ex)
print("\n")
print(test_matrix)

print("\nmatrix float")
matrix_float.set_column(0, vf1)
matrix_float.set_column(1, vf2)
matrix_float.set_column(2, vi1)
print(matrix_float)

print("\n matrix complex")
matrix_complex.set_column(0, vc1)
matrix_complex.set_column(1, vi1)
matrix_complex.set_column(2, vf1)
print(matrix_complex)

print("---------------------")
print("\nmatrix_int")
print(matrix_int)
print("\nmatrix_float")
print(matrix_float)
matrix_int+=matrix_float
print("\nresult matrix")
print(matrix_int)

print("----------------------")
print("\nscalar matrix")
sc_matrix = p.Matrix(3, 3)
print("\nmatrix int")
print(matrix_int)
sc_matrix=matrix_int
sc_matrix /= 45
print("\nsc matrix")
print(sc_matrix)
print("\nmatrix int")
print(matrix_int)

print("----------------------")
print("\nbinary operations with matrix")
bmatrix = p.Matrix(3, 3)
print("\nmatrix int")
print(matrix_int)
print("\nmatrix float")
print(matrix_float)
bmatrix = matrix_int + complex(10)
print(bmatrix)

print("----------------------")
print("\ntest comparison")
print(matrix_float == matrix_float)
print(matrix_float != matrix_float)

print("-----------------------")
print("\nmatrix functions")
print(matrix_float)
print("\n")
matrix_float.fill_column(0, 30)
matrix_float.fill_row(0, 50)
matrix_float.multiply_column(0, 2)
matrix_float.multiply_row(0, 4)
matrix_float.set_column(0, vex)
matrix_float.set_row(0, vf1)
matrix_float.add_to_column(0, vf2)
matrix_float.add_to_row(0, vf2)
matrix_float.subtract_from_column(0, vc1)
matrix_float.subtract_from_row(0, vex)
print(matrix_float)

print("------------------------")
print("maddie")
ronan = p.Matrix(10, 10)
ronan.set_to_unity()
print(ronan)
vex = p.Vector(complex(1, 1), 0, complex(2, 1), complex(5, 1))
maddie = p.Matrix(4, 4)
maddie.set_column(0, vex)
vex = p.Vector(complex(5, 1), complex(2, 1), complex(6, 1), complex(2, 1))
maddie.set_column(1, vex)
vex = p.Vector(complex(7, 1), complex(6, 1), complex(3, 1), complex(3, 1))
maddie.set_column(2, vex)
vex = p.Vector(complex(9, 1), complex(8, 1), complex(7, 1), complex(4, 1))
maddie.set_column(3, vex)
print(maddie)
print(maddie.trace())

print("-----------------------")
print("matrix multiply")
vex = p.Vector(1, 2, 3)
matrix_a = p.Matrix(3, 3)
matrix_b = p.Matrix(3, 3)
matrix_c = p.Matrix(3, 3)
for i in range(3):
	matrix_a.set_row(i, vex)
	matrix_b.set_row(i, vex)
for i in range(3):
	for i in range(3):
		matrix_a[i][i] = 0
matrix_b[0][1] = 0
matrix_b[1][0] = 0
matrix_b[1][2] = 0
print(matrix_a)
print(matrix_b)

p.multiply(matrix_a, matrix_b, matrix_c)
print(matrix_c)

print("------------------------")
print("matrix multiplyz")
print(matrix_a)
print(matrix_b)
scalar = 2.99
p.multiplyZ(matrix_a, matrix_b, matrix_c)
print(matrix_c)
print("------------------------")
print("matrix multiplyz by scalar")
for i in range(3):
	matrix_a.set_row(i, vex)
	matrix_b.set_row(i, vex)
print(matrix_a)
print(matrix_b)
scalar = 2.5
p.multiplyZ(matrix_a, scalar, matrix_c)

print(matrix_c)
print("----------------------")
print("apl_ml")
vec = p.Vector(0, 0, 0)
scalar = 2.99
for i in range(3):
	matrix_a.set_row(i, vex)
	matrix_b.set_row(i, vex)
	matrix_c.set_row(i, vec)

print(scalar)
print(matrix_a)
print(matrix_b)
print(matrix_c)
p.apl_ml(matrix_a, scalar, matrix_c)
print(matrix_c)
print("----------------------")
print("as_ml_ml")
print(matrix_a)
print(matrix_b)
print(matrix_c)
scalar = 2
p.as_ml_ml(matrix_a, matrix_b, scalar, matrix_c)
print(matrix_c)
print("----------------------")
print("test_matrix operations")
a_test = p.Matrix(3, 3)
b_test = p.Matrix(3, 3)
a_vec = p.Vector(1, 2, 3)

test_vec = p.Vector(0, 0, 0)
test_vec2 = p.Vector(0, 0, 0)
for i in range(3):
	a_test.set_row(i, a_vec)
	b_test.set_row(i, a_vec)
test_vec2 = a_test * a_vec
test_vec = a_test.mul_left(a_vec)

print(b_test.transpose())
print(test_vec)
print(test_vec2)

print("--------------------")
print("test matrix inverse")
print("a_test")
char_mat = p.Matrix(3, 3)
vec = p.Vector(1, 2, 3)
a_test = p.Matrix(3, 3)
a_test.set_row(0, vec)
char_mat.set_row(0, vec)
vec = p.Vector(3, 2, 1)
a_test.set_row(1, vec)
char_mat.set_row(1, vec)
vec = p.Vector(2, 1, 3)
a_test.set_row(2, vec)
char_mat.set_row(2, vec)
print(a_test)
p.GaussInvert(a_test)
print(a_test)

print("b_test")
vec = p.Vector(1, 2, 3)
b_test = p.Matrix(3, 3)
b_test.set_row(0, vec)
vec = p.Vector(0, 5, 4)
b_test.set_row(1, vec)
vec = p.Vector(0, 4, 6)
b_test.set_row(2, vec)

print(b_test)
p.GaussInvert(b_test)
print(b_test)

print("-----------------")

test_mat = p.Matrix(3, 3)
maddie_mat = p.Matrix(3, 3)
print(char_mat)
p.invert(char_mat, test_mat)
print(char_mat)
print(test_mat)
print("inverse function")
print("b_test")
print("maddie_mat")
print(maddie_mat)
maddie_mat = p.inverse(char_mat)
print(char_mat)
print("maddie_mat")
print(maddie_mat)

print("------------------")
vec = p.Vector(1, 2, 3)
vex = p.Vector(3, 2, 1)

maddie_mat = vex%vec
print(maddie_mat)
