import time
import PyTorus as p

st = time.time()
print(st)
maddie_mat = p.Matrix(10, 10)
# 1
vec = p.Vector(1, 3, 4, 5, 6, 56,4,5,1,36)
maddie_mat.set_row(0, vec)
# 2
vec = p.Vector(6,5,8,4,33, 8,45,1,2,3)
maddie_mat.set_row(1, vec)
# 3
vec = p.Vector(8,45,1,2,3, 6,5,8,4,33)
maddie_mat.set_row(2, vec)
# 4
vec = p.Vector(38,6,43,43,2, 8,45,1,2, 3)
maddie_mat.set_row(3, vec)
# 5
vec = p.Vector(56,4,5,1,36,38,6,43,43,2)
maddie_mat.set_row(4, vec)

#6
vec = p.Vector(1,36,38,6,43, 1,2,3, 6,5)
maddie_mat.set_row(5, vec)
#7
vec = p.Vector(30, 45,1,2,3, 6, 5,8,4,33)
maddie_mat.set_row(6, vec)
#8
vec = p.Vector(4, 5, 6, 56,4, 8,45,1,2,3)
maddie_mat.set_row(7, vec)
#9
vec = p.Vector(38,6,43,43,2, 2, 8,45,1,2)
maddie_mat.set_row(8, vec)
#10
vec = p.Vector(5,8,4,33, 8, 43,43,2, 2, 8)
maddie_mat.set_row(9, vec)

print(maddie_mat)

p.GaussInvert(maddie_mat)

print(maddie_mat)

et = time.time()
elapsed_time = et - st
print(et)
print('Execution time:', elapsed_time, 'seconds')

