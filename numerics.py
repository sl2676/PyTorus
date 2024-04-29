import PyTorus as p

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


for i in range(3):
	maddie_mat.set_row(i, vex)
print(maddie_mat)
val = p.GaussJordan(maddie_mat, 3, a_test, 3)
print("GaussJordan1")
print(val)
print("GaussJordan2")

vec = p.Vector(0, 1, 0)
vex = p.Vector(3, 3, 3)
'''
for i in range(3):
	maddie_mat.set_row(i, vex)
val = p.GaussJordan(maddie_mat, 3, vec)
print(val)
'''
'''
val = p.GaussJordan(maddie_mat, 3, a_test, 3)
print(val)
'''
