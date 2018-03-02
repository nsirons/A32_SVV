import matplotlib.pyplot as plt
from math import *
import numpy as np
from bendingstresses import discretize_skin
from tools import heaviside
from bendingdeflection import get_deflection_y

def find_shear_stresses(x, discretized_skin_pos, la, x1, x2, x3, xa, d1, d3, Ca, ha, G,
                          tsp, tsk, dact1, dact2, Izz, Iyy, Izy, ybar, zbar,
                          theta, Fz2, Fy1, Fy2, Fy3, Fx1, Fx3, Fz1,
                        FzI, P, q, S_y, S_z):

    #input variables
    x2a = x2-xa/2 #x-location actuator 1#
    x2b = x2+xa/2 #x-location actuator 2#
    phi = pi #angle wall 12outer

    #parameter simplifying formula
    c = sqrt((Ca-ha/2)**2+(ha/2)**2)
    yp = ha/2*cos(theta)

    #cell areas
    A1 = pi*ha**2/8
    A2 = ha*(Ca-ha/2)/2

    A = np.matrix([[2*A1, 2*A2, 0],
               [(ha/tsp+pi*ha/(4*tsk))/(2*G*A1), -ha/(2*G*tsp*A2), -1],
               [-ha/(G*tsp*2*A1), (2*c/tsk+ha/tsp)/(G*2*A2), -1]])

    #finding array of coordinates
    positions = discretized_skin_pos

    shearstress = []

    for i in range(len(positions)):
        #shear forces
        Sz = S_z(x)
        Sy = S_y(x)

        # qb as a function of theta/s for each wall, q12_I init for I, q_23 init for II, others q2+q1(end)
        k1 = -(Sz * Izz - Sy * Izy) / (Izz * Iyy - Izy ** 2)
        k2 = -(Sy * Iyy - Sz * Izy) / (Izz * Iyy - Izy ** 2)

        phi = pi
        s12 = ha
        s23 = c

        q12_Ib = k1*tsk*(ha/2)**2*np.cos(phi) + k2*tsk*(ha/2)**2*np.sin(phi)

        q21b = k2*tsp*(-ha/2*s12+0.5*s12**2) + q12_Ib

        q23b = k1*tsk*(Ca-ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2 / 2 + \
                        k2*tsk*(-ha/2*s23 + (ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2/2)

        q31b = k1*tsk*((Ca-ha/2)*s23 -(Ca-ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2 / 2) + \
                        k2*tsk*(ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2/2 + q23b*(np.hypot((Ca-ha/2), ha/2))

        q12_IIb = k2*tsp*(ha/2*s12-s12**2 / 2) + q31b

        #b matrix (internal torque, torque due to variable shear flows, rate of twist due to variable shear flows (in cell 1 & 2))
        T = +q*cos(theta)*((Ca/4 - ha/2) * cos(theta))*x\
             + Fz1*(d1- get_deflection_y(x_position, moment_y, moment_z, Izz, Iyy, Izy))*heaviside(x-x1) \
             + FzI*((get_deflection_y((x_position - xa/2), moment_y, moment_z, Izz, Iyy, Izy) + ha/sqrt(2) *cos(radians(135)+theta) ) - get_deflection_y(x_position, moment_y, moment_z, Izz, Iyy, Izy))*heaviside(x - (x2 - xa/2)) \
             + Fz2 * (0 - get_deflection_y(x_position, moment_y, moment_z, Izz, Iyy, Izy))*heaviside(x-x2) \
             - P*(get_deflection_y(x_position, moment_y, moment_z, Izz, Iyy, Izy) - (get_deflection_y((x_position + xa/2), moment_y, moment_z, Izz, Iyy, Izy) + ha/sqrt(2) *cos(radians(135)+theta)))*heaviside(x - (x2 + xa/2))

        Tvar = -(q23b+q31b)*(ha/2*sin(atan((Ca-ha/2)/(ha/2))))-2*q12_Ib*ha/2

        dthetadxc1 = (q12_Ib/tsk + q21b/tsp)/(2*A1*G)
 
        dthetadxc2 = (q23b/tsk + q31b/tsp + q12_IIb/tsp)/(2*A2*G)

        b = np.matrix([[T+Tvar],
                       [dthetadxc1],
                       [dthetadxc2]])
        sol=np.linalg.solve(A,b)

        #calculating variable shear flow at each position
        if positions[i][1] == 0:
            phi = 0
        else:
            phi = tan(positions[i][0]/positions[i][1]) #angle wall 12outer

        s23 = c - ((positions[i][0])**2 + (positions[i][1])**2)
        s31 = ((positions[i][0])**2 + (positions[i][1])**2)

        q12_Ibloc = k1*tsk*(ha/2)**2*np.cos(phi) + k2*tsk*(ha/2)**2*np.sin(phi)

        q23bloc = k1*tsk*(Ca-ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2 / 2 + \
                        k2*tsk*(-ha/2*s23 + (ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2/2)

        q31bloc = k1*tsk*((Ca-ha/2)*s31 -(Ca-ha/2)/np.hypot((Ca-ha/2), ha/2)*s31**2 / 2) + \
                        k2*tsk*(ha/2)/np.hypot((Ca-ha/2), ha/2)*s31**2/2 + q31b*(np.hypot((Ca-ha/2), ha/2))

        #adding constant shear flow in cell to local variable shear flow
        q12 = sol[0] + q12_Ibloc
        q23 = sol[1] + q23bloc
        q31 = sol[1] + q31bloc

        #dividing by correct thickness to get shear stress
        if positions[i][0] > ha/2 and positions[i][1] > 0:
            shearstress.append(float(q31)/tsk)
        if positions[i][0] > 0:
            shearstress.append(float(q12)/tsk)
        elif positions[i][0] < 0 and positions[i][1] < 0:
            shearstress.append(float(q23)/tsk)
        else:
            shearstress.append(0)
    return shearstress
