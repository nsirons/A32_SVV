from math import *
from bendingstresses import find_bending_stresses
from shearstressestorshe import find_shear_stresses
import os
import matplotlib.pyplot as plt
import numpy as np

#print(os.path.abspath('.'))

la = 1.691 #length aileron
n = 3 #discretization cross-section
d = 5 #discretization along span
x = 0 #initial x-coordinate
dx = la/(d+1) #steps along span

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


# TODO: Add MoI
I_yy = 0.0002  # change to functions
I_zz = 0.0002  # change to functions
I_zy = 0  # change to functions
ybar = 0  # change to functions
zbar = 0.15  # change to functions

# Calculate reaction forces
#x = 0#reaction_forces(C_a, l_a, h_a, x_a, x_1, x_2, x_3, d_1, d_3, d_act1, d_act2, theta, E, P, q, I_yy, I_zz, I_zy)
F_y1 = 1000
F_y2 = 1000
F_y3 = 1000
F_z1 = 1000
F_z2 = 1000
F_zI = 1000
cst_1 = 0  # what is cst?
cst_2 = 0

#stresses
sigmay = 0
sigmaz = 0
tauxy = 0
tauzx = 0

threedsigmamax = []
xlst = []
sigmamaxlst = []

while x <= la:   

    for i in range(n):
        
        x = find_bending_stresses(x, n, l_a, x_1, x_2, x_3, x_a, d_1, d_3,
                                                     I_zz, I_yy, I_zy, ybar, zbar,
                                                     F_y1, F_y2, F_y3, 0, 0, F_z1, F_zI, P, q)[0][i]
        xlst.append(x)
        y = find_bending_stresses(x, n, l_a, x_1, x_2, x_3, x_a, d_1, d_3,
                                                     I_zz, I_yy, I_zy, ybar, zbar,
                                                     F_y1, F_y2, F_y3, 0, 0, F_z1, F_zI, P, q)[1][i]
        z = find_bending_stresses(x, n, l_a, x_1, x_2, x_3, x_a, d_1, d_3,
                                                     I_zz, I_yy, I_zy, ybar, zbar,
                                                     F_y1, F_y2, F_y3, 0, 0, F_z1, F_zI, P, q)[2][i]

        sigmax = find_bending_stresses(x, n, l_a, x_1, x_2, x_3, x_a, d_1, d_3,
                                                     I_zz, I_yy, I_zy, ybar, zbar,
                                                     F_y1, F_y2, F_y3, 0, 0, F_z1, F_zI, P, q)[3][i]
        tauyz = find_shear_stresses(x, n, l_a, x_1, x_2, x_3, x_a, d_1, d_3, C_a, h_a, G,
                          t_sp, t_sk, d_act1, d_act2, I_zz, I_yy, I_zy, ybar, zbar,
                          theta, F_z2, F_y1, F_y2, F_y3, 0, 0, F_z1,
                        F_zI, P, q)[i]
        
        sigmamax = sqrt(1/2*((sigmax-sigmay)**2 + (sigmay-sigmaz)**2 + (sigmaz-sigmax)**2)) + sqrt(3*(tauxy**2 + tauyz**2 + tauzx**2))
        sigmamaxlst.append(sigmamax)

        threedsigmamax.append([x,z,y,sigmamax])
        
    x = x + dx

print(threedsigmamax)

plt.plot(xlst, sigmamaxlst)
plt.show()
