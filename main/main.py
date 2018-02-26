# constants/geometry
from reaction_forces import reaction_forces
from bendingstresses import find_bending_stresses
from shearstressestorshe import find_shear_stresses
import numpy as np
from math import *
import matplotlib.pyplot as plt
from aileron import aileron
from MOI import calculate_inertia_yy, calculate_inertia_zz, calculate_rotated_inertia

# -- Geometry --
C_a = 0.484  # chord length aileron
l_a = 1.691  # span aileron
x_1 = 0.149  # x-location hinge 1
x_2 = 0.554  # x-location hinge 2
x_3 = 1.541  # x-location hinge 3
x_a = 27.2e-2  # distance actuator 1 - 2
h_a = 17.3e-2  # aileron height
t_sk = 1.1e-3  # spar thickness
t_sp = 1.2e-3  # skin thickness
h_st = 1.4e-2  # Height of stiffener
w_st = 1.8e-2  # Width of stiffener
n_st = 13  # number of stringers
d_1 = 6.81e-2  # ver. displacement hinge 1
d_3 = 20.3e-2  # ver. displacement hinge 3
d_act1 = 0  # ? # iterative*
d_act2 = 0  # ? # iterative*
theta = 26 * np.pi/180  # Maximum upward deflection

# -- Loads --
P = 37.9e3  # load in actuator 2
q = 27.1e3  # net aerodynamic load

# -- Material Property
# source: http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t3
E = 73.1e9  # Young modulus aluminium 2024-T3
G = 28.e9  # Shear modulus aluminium 2024-T3


# MoI
aileron_obj = aileron()

inertia_uu = calculate_inertia_zz(aileron)
inertia_vv = calculate_inertia_yy(aileron)

Ivec = calculate_rotated_inertia(inertia_uu, inertia_vv, 0, -26)

I_zz = Ivec[0]
I_yy = Ivec[1]
I_zy = Ivec[2]
ybar = 0#None  # change to functions
zbar = 0.15#None  # change to functions

# Calculate reaction forces
x = reaction_forces(C_a, l_a, h_a, x_a, x_1, x_2, x_3, d_1, d_3, d_act1, d_act2, theta, E, P, q, I_yy, I_zz, I_zy)
F_y1 = x[0][0]
F_y2 = x[1][0]
F_y3 = x[2][0]
F_z1 = x[3][0]
F_z2 = x[4][0]
F_zI = x[5][0]
cst_1 = x[6][0]  # what is cst?
cst_2 = x[7][0]

# Calculate normal stress along the cross-sectional area and x direction
d = 50 #discretization along span
x = 0 #initial x-coordinate
dx = l_a/(d+1) #steps along span
n = 3  # number of discretized points
x_points = []
y_points = []
z_points = []
sigma_points = []
threedsigmamax = []

#stresses
sigmay = 0
sigmaz = 0
tauxy = 0
tauzx = 0

while x <= l_a:  

    for i in range(n):
        
        xpos = find_bending_stresses(x, n, l_a, x_1, x_2, x_3, x_a, d_1, d_3,
                                                     I_zz, I_yy, I_zy, ybar, zbar,
                                                     F_y1, F_y2, F_y3, 0, 0, F_z1, F_zI, P, q)[0][i]
        ypos = find_bending_stresses(x, n, l_a, x_1, x_2, x_3, x_a, d_1, d_3,
                                                     I_zz, I_yy, I_zy, ybar, zbar,
                                                     F_y1, F_y2, F_y3, 0, 0, F_z1, F_zI, P, q)[1][i]
        zpos = find_bending_stresses(x, n, l_a, x_1, x_2, x_3, x_a, d_1, d_3,
                                                     I_zz, I_yy, I_zy, ybar, zbar,
                                                     F_y1, F_y2, F_y3, 0, 0, F_z1, F_zI, P, q)[2][i]
        sigmax = find_bending_stresses(x, n, l_a, x_1, x_2, x_3, x_a, d_1, d_3,
                                                     I_zz, I_yy, I_zy, ybar, zbar,
                                                     F_y1, F_y2, F_y3, 0, 0, F_z1, F_zI, P, q)[3][i]
        tauyz = find_shear_stresses(x, n, l_a, x_1, x_2, x_3, x_a, d_1, d_3, C_a, h_a, G,
                          t_sp, t_sk, d_act1, d_act2, I_zz, I_yy, I_zy, ybar, zbar,
                          theta, F_z2, F_y1, F_y2, F_y3, 0, 0, F_z1,
                        F_zI, P, q)[i]
        print(sigmax)
        
        sigmamax = sqrt(1/2*((sigmax-sigmay)**2 + (sigmay-sigmaz)**2 + (sigmaz-sigmax)**2)) + sqrt(3*(tauxy**2 + tauyz**2 + tauzx**2))
        x_points.append(xpos)
        y_points.append(ypos)
        z_points.append(zpos)
        sigma_points.append(sigmamax)

        threedsigmamax.append([xpos,zpos,ypos,sigmamax])
    x = x + dx

print(threedsigmamax)
plt.plot(x_points, sigma_points)
plt.show()
# TODO: add plotters
# TODO: Add shear and Torsion
# TODO: Add Von Misses
# TODO: Ribs
# TODO: Numerical Part
# TODO: Verification and unit-tests for bending/shear/torsion
# TODO: Validation
