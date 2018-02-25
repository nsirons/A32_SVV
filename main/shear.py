import numpy as np


def calc_shear(Ca,ha,t_skin, t_spar, theta, G,
               Izz, Iyy, Izy,
               P, q, Fz, Fy):

    # TODO: Unit-tests
    reaction_vec = np.reshape(np.array([-Fz, Fy]) @ np.array([[np.cos(theta), -np.sin(theta)],
                                                              [np.sin(theta), np.cos(theta)]]), (2, 1))

    q_vec = np.reshape(-q * np.array([np.sin(theta), np.sin(theta)]), (2, 1))
    p_vec = np.reshape(P * np.array([np.cos(theta), -np.sin(theta)]), (2, 1))
    Sx, Sy = q_vec + p_vec + reaction_vec

    # for local coordinate system
    k1 = -(Sx * Izz - Sy * Izy) / (Izz * Iyy - Izy ** 2)
    k2 = - (Sy * Iyy - Sx * Izy) / (Izz * Iyy - Izy ** 2)

    # qb as a function of theta/s for each wall, q12_I init for I, q_23 init for II, others q2+q1(end)
    q12_I = lambda theta: k1*t_skin*(ha/2)**2*np.cos(theta) + k2*t_skin*(ha/2)**2*np.sin(theta)

    q21 = lambda s: k2*t_spar*(-ha/2*s+0.5*s**2) + q12_I(np.pi/2)

    q23 = lambda s: k1*t_skin*(Ca-ha/2)/np.hypot((Ca-ha/2), ha/2)*s**2 / 2 + \
                    k2*t_skin*(-ha/2*s + (ha/2)/np.hypot((Ca-ha/2), ha/2)*s**2/2)

    q31 = lambda s: k1*t_skin*((Ca-ha/2)*s -(Ca-ha/2)/np.hypot((Ca-ha/2), ha/2)*s**2 / 2) + \
                    k2*t_skin*(ha/2)/np.hypot((Ca-ha/2), ha/2)*s**2/2 + q23(np.hypot((Ca-ha/2), ha/2))

    q12_II = lambda s: k2*t_spar*(ha/2*s-s**2 / 2) + q31(np.hypot((Ca-ha/2), ha/2))


    # Calculate Moment around hinge, caused by Forces
    q_arm = np.reshape(np.array([0.25*Ca-ha/2, 0]), (1, 2))
    p_arm = np.reshape(np.array([-ha/2, ha/2]), (1, 2))
    left_side = q_arm @ q_vec + p_arm @ p_vec

    # Area of cells
    A_I = np.pi * (ha/2)**2 / 2
    A_II = (Ca-ha/2)*ha/2 / 2

    # Calculate Moment around hinge, caused by shear flow
    p0_2 = (Ca-ha/2)*ha/2 / np.hypot(Ca-ha/2, ha/2)  # altitude of right triangle
    M1 = ha**3*t_skin*(k1*np.sin(np.pi/2)-k2*np.cos(np.pi/2)) / 8 + 0  # p0 = ha/2
    s = np.hypot((Ca-ha/2), ha/2)
    d = np.hypot((Ca-ha/2), ha/2)
    M2 = -p0_2*t_skin * (s * (-2 * Ca * k1 + ha* k2) * np.sqrt(4 * Ca ** 2 - 4 * Ca * ha + 2 * ha ** 2) + (k1 + k2) * ha ** 3 -
                         4 * Ca * (k1 + (1 / 2) * k2) * ha ** 2 +
                         ((6 * k1 + 2 * k2) * Ca ** 2 -
                          4 * k2 * s ** 2 * (1 / 3)) * ha - 4 * Ca ** 3 * k1) * s / \
         (4 * np.sqrt(4 * Ca ** 2 - 4 * Ca * ha + 2 * ha ** 2))

    # Solve System: 1st&2nd are dtheta/dz = int of q, 3rd Moment caused by forces + shearflow+torque = 0
    # TODO: Check&Change values of s1 and s2
    A = np.array([[(theta+s)/(A_I*G), -s/(A_I*G), -1],
                  [-s/(A_II*G), 3*s/(A_II*G), -1],
                  [2*A_I, 2*A_II, 0]])

    b = np.array([[-1/(A_I*G)*(0.25*t_skin*ha**2 * (k1*np.sin(theta)-k2*np.cos(theta)) + 0.25*k2*t_skin*ha**2*s+
                   k2*t_spar*(-0.25*ha*s**2+1/6*s**3) - 0.5*k1*t_skin*(Ca-ha/2)*d*s -k1*t_skin*(Ca*d-0.5*(Ca-ha/2)*d)*s -
                   -k2*t_spar*(0.25*ha*s**2-1/6*s**3))],
                  [-1/(A_II*G)*k1*t_skin*(Ca-(1/2)*ha)*s**3/(3*np.sqrt(4*Ca**2-4*Ca*ha+2*ha**2))+
                   k2*t_skin*(-(1/4)*ha*s**2+ha*s**3/(6*np.sqrt(4*Ca**2-4*Ca*ha+2*ha**2)))+
                   (1/2)*k1*t_skin*(Ca-(1/2)*ha)*np.sqrt(4*Ca**2-4*Ca*ha+2*ha**2)*s-
                   (1/8)*k2*t_skin*ha*np.sqrt(4*Ca**2-4*Ca*ha+2*ha**2)*s+
                   k1*t_skin*((1/2)*Ca*s**2-(Ca-(1/2)*ha)*s**3/(3*np.sqrt(4*Ca**2-4*Ca*ha+2*ha**2)))+
                   k2*t_skin*ha*s**3/(6*np.sqrt(4*Ca**2-4*Ca*ha+2*ha**2))+k1*t_skin*((1/2)*Ca*np.sqrt(4*Ca**2-4*Ca*ha+2*ha**2)-
                                                                                      (1/4)*(Ca-(1/2)*ha)*np.sqrt(4*Ca**2-4*Ca*ha+2*ha**2))*s+
                   k2*t_spar*((1/4)*ha*s**2-(1/6)*s**3)-(1/4)*k2*t_skin*ha**2*s-k2*t_spar*(-(1/4)*ha*s**2+(1/6)*s**3)],
                  [left_side - M1 - M2]])

    x = np.linalg.solve(A, b)

    # cut_1, cut_2, dtheta/dz
    return x