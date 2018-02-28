import unittest
from main.reaction_forces import reaction_forces
import numpy as np


class TestMoI(unittest.TestCase):
    # -- Geometry --
    C_a = 0.484  # chord length aileron
    l_a = 2  # span aileron
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
    theta = 26 * np.pi / 180  # Maximum upward deflection

    # -- Loads --
    P = 37.9e3  # load in actuator 2
    q = 27.1e3  # net aerodynamic load

    # -- Material Property
    # source: http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t3
    E = 73.1e9  # Young modulus aluminium 2024-T3
    G = 28.e9  # Shear modulus aluminium 2024-T3

    I_zz = 1e-6
    I_yy = 2e-6
    I_zy = 0

    # x = reaction_forces(self.C_a, self.l_a, self.h_a, self.x_a, self.x_1, self.x_2, self.x_3,
    #                     self.d_1, self.d_3, self.d_act1, self.d_act2, self.theta,
    #                     self.E, self.P, self.q, self.I_zz, self.I_yy, self.I_yz)
    def test_zero_U(self):
        x = reaction_forces(self.C_a,self.l_a,self.h_a,self.x_a,self.x_1,self.x_2,self.x_3,
                            0.0,0.0,0,0, self.theta * 180/np.pi,
                            self.E,self.P,self.q,self.I_zz,self.I_yy,self.I_zy)

        F_y1 = x[0][0]
        F_y2 = x[1][0]
        F_y3 = x[2][0]
        F_z1 = x[3][0]
        F_z2 = x[4][0]
        F_zI = x[5][0]
        cst_1 = x[6][0]
        cst_2 = x[7][0]

        self.assertAlmostEqual(F_y1, 0)
        self.assertAlmostEqual(F_y2, 0)
        self.assertAlmostEqual(F_y3, 0)
        self.assertAlmostEqual(F_z1, 0)
        self.assertAlmostEqual(F_z2, 0)
        self.assertAlmostEqual(F_zI, 0)
        self.assertAlmostEqual(cst_1, 0)
        self.assertAlmostEqual(cst_2, 0)

    def other(self):
        x = reaction_forces(self.C_a, self.l_a, self.h_a, self.x_a, self.x_1, self.x_2, self.x_3,
                            self.d_1, self.d_3, self.d_act1, self.d_act2, self.theta,
                            self.E, self.P, self.q, self.I_zz, self.I_yy, self.I_yz)