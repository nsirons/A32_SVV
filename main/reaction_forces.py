import numpy as np


def reaction_forces(C_a, l_a, h_a, x_a, x_1, x_2, x_3, d_1, d_3, d_act1, d_act2, theta, E, P, q, I_zz, I_yy, I_yz):

    v_bot = E * ((I_zz * I_yy) - (I_yz ** 2))
    y_p = (h_a / 2) * np.cos(theta)
    z_q = (C_a / 4 - h_a / 2) * np.cos(theta)

    a = np.array([[1, 1, 1, 0, 0, 0, 0, 0],
                  [0, 0, 0, 1, 1, 1, 0, 0],
                  [0, 0, 0, d_1, 0, (y_p + d_act1), 0, 0],
                  [0, 0, 0, (x_2 - x_1), 0, x_a / 2, 0, 0],
                  [-(x_2 - x_1), 0, (x_3 - x_2), 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, (1 / 6) * I_yz / v_bot, x_1 / v_bot, 1 / v_bot],
                  [(I_yy * (x_2 - x_1) ** 3) / (6 * v_bot), 0, 0, ((1 / 6) * I_yz * (x_2 - x_1) ** 3) / v_bot, 0,
                   ((1 / 6) * I_yz * (x_a / 2) ** 3) / v_bot, x_2 / v_bot, 1 / v_bot],
                  [(I_yy * (x_3 - x_1) ** 3) / (6 * v_bot), (I_yy * (x_3 - x_2) ** 3) / (6 * v_bot), 0,
                   (I_yz * (x_3 - x_1) ** 3) / (6 * v_bot), (I_yz * (x_3 - x_2) ** 3) / (6 * v_bot),
                   (I_yz * (x_3 - x_2 + (x_a / 2)) ** 3) / (6 * v_bot), x_3 / v_bot, 1 / v_bot]])

    b = np.array([[q * l_a],
                  [P],
                  [P * (y_p + d_act2) + (q * l_a * z_q)],
                  [-0.5 * x_a * P],
                  [(l_a / 2 - x_2) * q * l_a],
                  [q * I_yy * (x_1 ** 4) / (24 * v_bot)],
                  [((q * I_yy * (x_2 ** 4)) / (24 * v_bot)) + ((P * I_yz * (x_a / 2) ** 3) / (6 * v_bot)) - d_1],
                  [((q * I_yy * (x_3 ** 4)) / (24 * v_bot)) + (
                          (P * I_yz * (x_3 - x_2 - (x_a / 2)) ** 3) / (6 * v_bot)) - d_1 + d_3]])

    x = np.linalg.solve(a, b)

    return x
