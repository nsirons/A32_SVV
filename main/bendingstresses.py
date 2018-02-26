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
                          Fy1, Fy2, Fy3, Fx1, Fx3, Fz1, FzI, P, q):

    x2a = x2 - xa / 2  # x-location actuator 1
    x2b = x2 + xa / 2  # x-location actuator 2

    # --------------------------------------BENDING------------------------------------------------

    # return lists
    xlst = []
    ylst = []
    zlst = []
    sigmaxlst = []

    # Steps (switching terms on/off)
    s1 = 1
    s2 = 1
    s2a = 1
    s2b = 1
    s3 = 1

    # x = 0
    # dx = 0.1 #discretization along length aileron
    # n = 50 #discretization cross-section

    positions = calculate_stringer_positions(n)

    for i in range(len(positions)):
        z = zbar - positions[i][0]
        y = positions[i][1]
        if x > x3:
            s3 = 0
            s2 = 0
            s2a = 0
            s2b = 0
            s1 = 0

        elif x > x2b and x < x3:
            s2 = 0
            s2a = 0
            s2b = 0
            s1 = 0

        elif x > x2 and x < x2b:
            s2 = 0
            s2a = 0
            s1 = 0

        elif x > x2a and x < x2:
            s2a = 0
            s1 = 0

        elif x > x1 and x < x2a:
            s1 = 0

        Mz = Fy3 * (x3 - x) * s3 + Fy2 * (x2 - x) * s2 + Fy1 * (x1 - x) * s1 - Fx3 * d3 * s3 \
             - Fx1 * d1 * s1 - ((la / 2) - x2) * q * la
        My = Fz1 * (x1 - x) * s1 + P * (x2b - x) * s2b + FzI * (x2a - x) * s2a
        # Mzlst.append(Mz)
        # Mylst.append(My)

        sigmax = ((Mz * Iyy - My * Izy) * y + (My * Izz - Mz * Izy) * z) / (Izz * Iyy - Izy ** 2)
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
