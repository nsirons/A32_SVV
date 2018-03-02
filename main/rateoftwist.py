import matplotlib.pyplot as plt
from math import *
import numpy as np
from bendingstresses import discretize_skin
from tools import heaviside

def rate_of_twist(x, rotated_discretized_skin_pos, la, x1, x2, x3, xa, d1, d3, Ca, ha, G,
                          tsp, tsk, dact1, dact2, Izz, Iyy, Izy, ybar, zbar,
                          theta, Fz2, Fy1, Fy2, Fy3, Fx1, Fx3, Fz1,
                        FzI, P, q, S_y, S_z):

    #to be replaced for deflection formula
    def def_y(x):
        if x>=0 and x<=x2:
            y = d1-(d1/(x2-0))*(x)
        elif x>x2 and x<=la:
            y = (d3/(la-x2))*(x-x2)
        return y

    #parameter simplifying formula
    c = sqrt((Ca-ha/2)**2+(ha/2)**2)
    yp = ha/2*cos(theta)

    #cell areas
    A1 = pi*ha**2/8
    A2 = ha*(Ca-ha/2)/2

    A = np.matrix([[2*A1, 2*A2, 0],
               [-ha/(G*tsp*2*A2), (2*c/tsk+ha/tsp)/(G*2*A2), -1],
               [(ha/tsp+pi*ha/(4*tsk))/(2*G*A1), -ha/(2*G*tsp*A1), -1]])

    #finding array of coordinates
    positions = rotated_discretized_skin_pos

    for i in range(len(positions)):

        Sz = S_z(x)
        Sy = S_y(x)

        # total qb for each wall, q12_I init for I, q_23 init for II, others q2+q1(end)
        k1 = -(Sz * Izz - Sy * Izy) / (Izz * Iyy - Izy ** 2)
        k2 = -(Sy * Iyy - Sz * Izy) / (Izz * Iyy - Izy ** 2)

        phi = pi
        s12 = ha
        s23 = c

        q12_Ib = k1*tsk*(ha/2)**2*np.cos(phi) + k2*tsk*(ha/2)**2*np.sin(phi)

        q21b = k2*tsp*(-ha/2*s12+0.5*s12**2) + q12_Ib*(np.pi/2)

        q23b = k1*tsk*(Ca-ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2 / 2 + \
                        k2*tsk*(-ha/2*s23 + (ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2/2)

        q31b = k1*tsk*((Ca-ha/2)*s23 -(Ca-ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2 / 2) + \
                        k2*tsk*(ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2/2 + q23b*(np.hypot((Ca-ha/2), ha/2))

        q12_IIb = k2*tsp*(ha/2*s12-s12**2 / 2) + q31b*(np.hypot((Ca-ha/2), ha/2))

        #b matrix (internal torque, torque due to variable shear flows, rate of twist due to variable shear flows (in cell 1 & 2))
        T = +q*cos(theta)*((Ca/4 - ha/2) * cos(theta))*x\
             + Fz1*(d1- def_y(x))*heaviside(x-x1) \
             + FzI*((def_y(x2 - xa/2) + ha/sqrt(2) *cos(radians(135)+theta) ) - def_y(x))*heaviside(x - (x2 - xa/2)) \
             + Fz2 * (0 - def_y(x) )*heaviside(x-x2) \
             + P*(def_y(x) - (def_y(x2 + xa/2) + ha/sqrt(2) *cos(radians(135)+theta)))*heaviside(x - (x2 + xa/2)) 

        Tvar = -(q23b+q31b)*(ha/2*sin(atan((Ca-ha/2)/(ha/2))))-2*q12_Ib*ha/2

        dthetadxc1 = (q12_Ib/tsk + q21b/tsp)/(2*A1*G)

        dthetadxc2 = (q23b/tsk + q31b/tsp + q12_IIb/tsp)/(2*A2*G)

        b = np.matrix([[T+Tvar],
                       [dthetadxc1],
                       [dthetadxc2]])

        sol=np.linalg.solve(A,b)

    return float(degrees(sol[2]))
