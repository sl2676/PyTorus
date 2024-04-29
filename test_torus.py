import PyTorus as p

print("dir pytorus")
print(dir(p))

print("\n--------------------")
print("\ndir pymatrix")
print(dir(p.Matrix))

print("\n--------------------")
print("\ndir pyvector")
print(dir(p.Vector))

print("\n--------------------")
print("\ndir cheby")
print(dir(p.Cheby))

print("\n-------------------")
print("\ndir toruserr")
print(dir(p.TorusExcept))


print("\n------------------")
print("src\\utils\n")
src_util = ["~CHB", "~Compress", "-Constants", "-Err", "FreeMemory", "~Matrix", "Numerics", "PJMCoords", "PJMNum", "PJM_cline", "PJM_utils", "-PJMebf", "-Pi", "Pspline", "Random", "Types", "-Units", "-Vector", "~WDMath", "~Jjb_utils"]
src_util_dict  = {}

for index, element in enumerate(src_util):
	src_util_dict[index] = element

print(src_util_dict)
for index in src_util_dict:
	print(f"{index} " + src_util_dict[index])
