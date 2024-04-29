import math as m
import PyTorus as p

vc1 = p.Vector(complex(1.2, 3), complex(1, 0), complex(4, 2))
vc2 = p.Vector(complex(3), complex(2), complex(1))
vi1 = p.Vector(round(m.sqrt(3)), 12), round(1/(m.sqrt(3)), 12), round(1/(m.sqrt(3)))
vi1 = p.Vector(5, 5, 5)
vi2 = p.Vector(1, 1, 1)
vf1 = p.Vector(1.0, 1.0, 1.0)
vf2 = p.Vector(3.4, 2.0, 1.0)
ve = p.Vector(complex(1), 2.3, 3.0)
# dot product $
# cross product $
# norm function
# arithmetic perator with scalar on left -> + $, - $, * $, \ $ 
# compress in and output
# division by zero function
# arithemtic operators with assign $
# all arithmetic operators 
# real
# imag
# arg
# conj 

print("ve vector")
print(ve)
print("-------")
print(dir(p))
print("-------")
print(dir(p.Matrix))
print("-------")
print(vc1)
print(vc2)
print("-------")
print("int vec")
print(vi1)
print(vi2)
print("------")
print("float vec")
print(type(vf1[1]))
print(vf2)
print("-------")

resdot1 = vc1 * vi2
resdot2 = vi1 * vi2
resdot3 = vf1 * vf2
print("resdot")
print(resdot1)
print(resdot2)
print(resdot3)

print("-------")
print("rescross")
rescross1 = vc1 ^ vc2
rescross2 = vi1 ^ vi2
rescross3 = vf1 ^ vf2

print(rescross1)
print(rescross2)
print(rescross3)

print("-------")
print("resscalara")
val = 0
resscalara1 = vc1 / val 
resscalara2 = vi1 / val
resscalara3 = vf1 / val
print(resscalara1)
print(resscalara2)
print(resscalara3)

print("-------")
print("norm")
normres = vi1.norm()
print(normres)


print("-------")
print("multiply elements")
result = vf1.multiply_elements(vf2)
print(result)

print("-------")
print("arg")
print(vc1.arg())
print("real")
print(vc1.real())
print("imag")
print(vc1.imag())
print("conj")
print(vc1.conj())

print("--------")

