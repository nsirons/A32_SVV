# constants/geometry
from reaction_forces import reaction_forces
from bendingstresses import find_bending_stresses
import numpy as np


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
I_yy = None  # change to functions
I_zz = None  # change to functions
I_zy = None  # change to functions
ybar = None  # change to functions
zbar = None  # change to functions

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
x = np.linspace(0, l_a, 100)
n = 50  # number of discretized points
x_points = []
y_points = []
z_points = []
sigma_points = []
for x_pos in x:
    xlst, ylst, zlst, sigmax = find_bending_stresses(x_pos, n, l_a, x_1, x_2, x_3, x_a, d_1, d_3,
                                                     I_zz, I_yy, I_zy, ybar, zbar,
                                                     F_y1, F_y2, F_y3, 0, 0, F_z1, F_zI, P, q)
    x_points.append(xlst)
    y_points.append(ylst)
    z_points.append(zlst)
    sigma_points.append(sigmax)

# TODO: add plotters
# TODO: Add shear and Torsion
# TODO: Add Von Misses
# TODO: Ribs
# TODO: Numerical Part
# TODO: Verification and unit-tests for bending/shear/torsion
# TODO: Validation
