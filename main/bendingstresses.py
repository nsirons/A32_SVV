import matplotlib.pyplot as plt
from math import *
import numpy as np

# finding array of coordinates
def calculate_stringer_positions(stiffener):
    height = 0.173
    chord = 0.484
    length_flat_skin = sqrt((height / 2) ** 2 + (chord - height / 2) ** 2)
    length_circular_skin = pi * height / 2
    angle_flat_skin = atan2(height / 2, chord - height / 2)

    length_total_skin = 2 * length_flat_skin + length_circular_skin
    distance_between_stringers = length_total_skin / (stiffener + 1)

    stringer_positions = []
    current_length = 0

    for i in range(1, stiffener + 1):
        current_length += distance_between_stringers
        position = []
        if 0 <= current_length <= length_flat_skin:
            x = chord - current_length * cos(angle_flat_skin)
            y = current_length * sin(angle_flat_skin)
            position = [x, y]
        elif length_flat_skin < current_length <= (length_flat_skin + length_circular_skin):
            current_angle = pi / 2 + (current_length - length_flat_skin) / (2 * pi * height / 2) * 2 * pi
            x = height / 2 + height / 2 * cos(current_angle)
            y = height / 2 * sin(current_angle)
            position = [x, y]
        elif (length_flat_skin + length_circular_skin) < current_length <= length_total_skin:
            modified_current_length = length_total_skin - current_length
            x = chord - modified_current_length * cos(angle_flat_skin)
            y = modified_current_length * sin(-angle_flat_skin)
            position = [x, y]
        else:
            position = [-1, -1]
        stringer_positions.append(position)
    return stringer_positions

# x <= la:
# x = float(input("enter a x-coordinate: "))

def find_bending_stresses(x, n, la, x1, x2, x3, xa, d1, d3,
                          Izz, Iyy, Izy, ybar, zbar,
                          My, Mz):
   
    xlst = []
    ylst = []
    zlst = []
    sigmaxlst = []

    positions = calculate_stringer_positions(n)

    for i in range(len(positions)):
        z = zbar - positions[i][0]
        y = positions[i][1]
       
        M_y = My(x)
        M_z = Mz(x)

        sigmax = ((-M_z * Iyy - M_y * Izy) * y + (M_y * Izz + M_z * Izy) * z) / (Izz * Iyy - Izy ** 2)
        # print(z,y,sigmax)
        sigmaxlst.append(sigmax)
        xlst.append(x)
        ylst.append(y)
        zlst.append(z)
    return xlst, ylst, zlst, sigmaxlst
    # x = x + dx

# print(sigmaxlst)
# print(sigmaxlst)
# print(zlst)
# print(len{zlst))
# plt.scatter(sigmaxlst, zlst)
# plt.axis([0.15, -0.5, -0.5, 0.5])
# plt.plot([x for x,y in positions], [y for x,y in positions])
# plt.show()
