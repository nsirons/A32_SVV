import numpy as np
# TODO: Unit-tests

def f2_local(Fz, Fy, theta):
    return np.reshape(np.array([-Fz, Fy]) @ np.array([[np.cos(theta), -np.sin(theta)],
                                                      [np.sin(theta), np.cos(theta)]]), (2, 1))

def calc_shear(Ca, ha, t_skin, t_spar, theta, G,
               Izz, Iyy, Izy,
               P, q, Fz, Fy):

    reaction_vec = f2_local(-Fz, Fy, theta)
    q_vec = f2_local(0, q, theta)
    p_vec = f2_local(P, 0, theta)
    
    Sx, Sy = q_vec + p_vec + reaction_vec

    # for local coordinate system
    k1 = -(Sx * Izz - Sy * Izy) / (Izz * Iyy - Izy ** 2)
    k2 = - (Sy * Iyy - Sx * Izy) / (Izz * Iyy - Izy ** 2)
    l_23 = np.hypot((Ca - ha / 2), ha / 2)

    # qb as a function of theta/s for each wall, q12_I init for I, q_23 init for II, others q2+q1(end)
    q12_I = lambda theta: k1 * t_skin * (ha / 2) ** 2 * np.cos(theta) + k2 * t_skin * (ha / 2) ** 2 * np.sin(theta)

    q21 = lambda s: k2 * t_spar * (-ha / 2 * s + 0.5 * s ** 2) + q12_I(np.pi / 2)

    q23 = lambda s: k1 * t_skin * (Ca - ha / 2) / l_23 * s ** 2 / 2 + \
                    k2 * t_skin * (-ha / 2 * s + (ha / 2) / l_23 * s ** 2 / 2)

    q31 = lambda s: k1 * t_skin * ((Ca - ha / 2) * s - (Ca - ha / 2) / l_23 * s ** 2 / 2) + \
                    k2 * t_skin * (ha / 2) / l_23 * s ** 2 / 2 + q23(l_23)

    q12_II = lambda s: k2 * t_spar * (ha / 2 * s - s ** 2 / 2) + q31(l_23)

    # Calculate Moment around hinge, caused by Forces
    q_arm = np.reshape(np.array([0.25 * Ca - ha / 2, 0]), (1, 2))
    p_arm = np.reshape(np.array([-ha / 2, ha / 2]), (1, 2))
    left_side = q_arm @ q_vec + p_arm @ p_vec

    # Area of cells
    A_I = np.pi * (ha / 2) ** 2 / 2
    A_II = (Ca - ha / 2) * ha / 2 / 2

    # Calculate Moment around hinge, caused by shear flow
    p0_2 = (Ca - ha / 2) * ha / 2 / l_23  # altitude of right triangle

    theta = np.pi / 2  # along the semicircle
    s_skin = np.hypot((Ca - ha / 2), ha / 2)  # along the wall_23 = wall_31
    s_spar = ha  # spar height
    M1 = ha ** 3 * t_skin * (k1 * np.sin(theta) - k2 * np.cos(theta)) / 8 + 0  # p0 = ha/2

    M2 = -t_skin * (((ha - 2 * Ca) * k1 - 3 * ha * k2) * s_skin + 6 * l_23 * (-Ca * k1 + ha * k2)) * s_skin ** 2 / (
            12 * l_23)  # from Maple

    # Solve System: 1st&2nd are dtheta/dz = int of q, 3rd Moment caused by forces + shearflow+torque = 0
    A = np.array([[(theta + s_skin) / (A_I * G), -s_spar / (A_I * G), -1],
                  [-s_spar / (A_II * G), (2 * s_skin + s_spar) / (A_II * G), -1],
                  [2 * A_I, 2 * A_II, 0]])

    # minuse sign, because it moved from left side to right side
    b = np.array([[-1 / (A_I * G) * (0.25 * t_skin * ha ** 2 * (k1 * np.sin(theta) - k2 * np.cos(theta) + k2 * s_spar) +
                                     k2 * t_spar * (
                                             -0.25 * ha * s_spar ** 2 + 1 / 6 * s_spar ** 3) - 0.5 * k1 * t_skin * (
                                             Ca - ha / 2) * s_skin ** 2 * s_spar / l_23 -
                                     k2 * t_skin * (-0.5 * ha * s_skin / 2 + 0.25 * ha * s_skin ** 2 / l_23) -
                                     k1 * t_skin * (Ca * s_skin - 0.5 * (Ca - ha / 2) * s_skin ** 2 / l_23) -
                                     0.25 * k2 * t_skin * ha * s_skin ** 2 * s_spar / l_23 - k2 * t_spar * (
                                             0.25 * ha * s_spar ** 2 - 1 / 6 * s_spar ** 3))],

                  [-1 / (A_II * G) * k1 * t_skin * (Ca - (1 / 2) * ha) * s_skin ** 3 / (3 * l_23) + 2 * k2 * t_skin * (
                          -(1 / 4) * ha * s_skin ** 2 + ha * s_skin ** 3 / (12 * l_23))
                   + k1 * t_skin * ((1 / 2) * Ca * s_skin ** 2 - (Ca - (1 / 2) * ha) * s_skin ** 3 / (6 * l_23)) +
                   t_skin * s_skin ** 3 * ha * k2 / (12 * l_23) + k1 * t_skin * (Ca - (1 / 2) * ha) * s_skin ** 2 * s_spar / (2 * l_23) +
                   k2 * t_skin * (-(1 / 2) * ha * s_skin + ha * s_skin ** 2 / (4 * l_23)) * s_spar +
                   k1 * t_skin * (Ca * s_skin - (Ca - (1 / 2) * ha) * s_skin ** 2 / (2 * l_23)) * s_spar + k2 * t_skin * ha * s_skin ** 2 * s_spar / (4 * l_23) +
                   k2 * t_spar * ((1 / 4) * ha * s_spar ** 2 - (1 / 6) * s_spar ** 3) - (1 / 4) * k2 * t_skin * ha ** 2 * s_spar -
                   k2 * t_spar * (-(1 / 4) * ha * s_spar ** 2 + (1 / 6) * s_spar ** 3)],

                  [left_side - M1 - M2]])

    x = np.linalg.solve(A, b)

    # cut_1, cut_2, dtheta/dz, shear flow function
    return x, q12_I, q21, q23, q31, q12_II
