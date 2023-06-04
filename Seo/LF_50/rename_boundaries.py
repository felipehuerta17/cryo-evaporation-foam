# rename_boundaries.py
# Renames the boundary patches to wall and wedge, respectively
# For the cylindrical 2D mesh

file = open("constant/polyMesh/boundary", "r")
a = file.readlines()
w = a.copy()

# Define the lines where "patch" will be replaced to be consistent with the nomenclature
p2wedge_lines = [22, 23, 50, 51]
p2wall_lines = [29, 30, 36, 37, 43, 44]

for i in p2wedge_lines:
    w[i] = a[i].replace("patch", "wedge")

for i in p2wall_lines:
    w[i] = a[i].replace("patch", "wall")

file.close()

file = open("constant/polyMesh/boundary", "w")

# Write the new files
for line in w:
    file.write(line)
file.close()
