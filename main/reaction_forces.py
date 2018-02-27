import numpy as np
from math import sin, cos, radians, sqrt


def reaction_forces(ca, la, ha, xa, x1, x2, x3, d1, d3, d_act1, d_act2, angle, E, p, q, Izz, Iyy, Izy):


    c = -1/(E*(Izz*Iyy - (Izy**2)))

    # x =        [Ay Az By Bz Cy Dz K1 K2]T
    A = np.array([[1, 0, 1, 0, 1, 0, 0, 0],
                 [0, 1, 0, 1, 0, 1, 0, 0],
                 [0, d1,0, 0, 0, sin(radians(135+angle))*ha/sqrt(2),0, 0],
                 [0, x2-x1,0,0,0,xa/2., 0, 0],
                 [-(x2-x1),0,0,0,(x3-x2),0,0,0],
                 [0, 0, 0, 0, 0, 0, x1*c, c],
                 [-c*((x2-x1)**3)/6.*Iyy, -c*((x2-x1)**3)/6.*Izy, 0, 0, 0, -c*((x2 - (x2 - xa/2.))**3)/6.*Izy, x2*c, c],
                 [-c*((x3-x1)**3)/6.*Iyy, -c*((x3-x1)**3)/6.*Izy, -c*((x3-x2)**3)/6.*Iyy, -c*((x3-x2)**3)/6.*Izy, 0, -c*((x3 - (x2 - xa/2.))**3)/6.*Izy, x3*c, c]])

    C = np.array([[-q*la],
                 [-p],
                 [-q*la*(ca/4 - ha/2)*cos(radians(angle)) - p*sin(radians(135+angle))*ha/sqrt(2)],
                 [p*xa/2],
                 [-q*la*(la/2 - x2)],
                 [c*q*(x1**4)/24.*Iyy],
                 [c*q*(x2**4)/24.*Iyy],
                 [c*q*(x3**4)/24.*Iyy + c*p/6*((x3-(x2+xa/2.))**3)*Izy]])

    B = np.array([[0],[0],[0],[0],[0],[0],[-d1],[-d1+d3]])
    T = B-C
    x=np.linalg.solve(A,T)
    return x