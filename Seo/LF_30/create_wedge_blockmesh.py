## Find the coordinates to create an uniform wedge
import numpy as np

# Ingress cylinder properties
diameter = 0.200871  # m
height = 0.213 * 0.7  # m
alpha = np.pi * 5 / 180  # Default 5 degrees mesh good practice
# Calculate vertices
x = (diameter / 2) * np.cos(alpha / 2)
y_plus = (diameter / 2) * np.sin(alpha / 2)
y_minus = -y_plus
# Input vertex lines
lines = [21, 22, 23, 24, 25, 26]
# Read the file and store in an array
file = open("system/blockMeshDict", "r")
a = file.readlines()

a[lines[0]] = "    (0 0 0)\n"
a[lines[1]] = (
    "    (" + format(x, ".5f") + " " + format(y_plus, ".5f") + " " + str(0) + ")\n"
)
a[lines[2]] = (
    "    ("
    + format(x, ".5f")
    + " "
    + format(y_plus, ".5f")
    + " "
    + format(height, ".5f")
    + ")\n"
)
a[lines[3]] = "    (" + str(0) + " " + str(0) + " " + str(height) + ")\n"
a[lines[4]] = (
    "    (" + format(x, ".5f") + " " + format(y_minus, ".5f") + " " + str(0) + ")\n"
)
a[lines[5]] = (
    "    ("
    + format(x, ".5f")
    + " "
    + format(y_minus, ".5f")
    + " "
    + format(height, ".5f")
    + ")\n"
)

file.close()
file = open("system/blockMeshDict", "w")

# Write the new files
for line in a:
    file.write(line)
file.close()
