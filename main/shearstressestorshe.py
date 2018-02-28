import matplotlib.pyplot as plt
from math import *
import numpy as np
from bendingstresses import discretize_skin

def find_shear_stresses(x, discretized_skin_pos, la, x1, x2, x3, xa, d1, d3, Ca, ha, G,
                          tsp, tsk, dact1, dact2, Izz, Iyy, Izy, ybar, zbar,
                          theta, Fz2, Fy1, Fy2, Fy3, Fx1, Fx3, Fz1,
                        FzI, P, q):

    #input variables
    x2a = x2-xa/2 #x-location actuator 1#
    x2b = x2+xa/2 #x-location actuator 2#
    phi = pi #angle wall 12outer

    #discretizing
    #n = 20 #number of sections
    #dx = la/(n+1) #discretization step along length aileron
    #x = 0

    #parameter simplifying formula
    c = sqrt((Ca-ha/2)**2+(ha/2)**2)
    yp = ha/2*cos(theta)

    #cell areas
    A1 = pi*ha**2/8
    A2 = ha*(Ca-ha/2)/2

    #Steps (switching terms on/off)
    s1 = 1
    s2 = 1
    s2a = 1
    s2b = 1
    s3 = 1
    sa = 1
    sb = 1

    twistsec = 0

    twistseclst = []
    xlst = []

    #torsion
    #T = - P*(yp+dact2)*s2b + FzI*(yp+dact1)*s2a + Fz1*d1*s1 - ((Ca/4-ha/2)*cos(theta))*q*x

    A = np.matrix([[2*A1, 2*A2, 0],
               [-ha/(G*tsp*2*A2), (2*c/tsk+ha/tsp)/(G*2*A2), -1],
               [(ha/tsp+pi*ha/(4*tsk))/(2*G*A1), -ha/(2*G*tsp*A1), -1]])

    #b = np.matrix([[T],
               #[0],
               #[0]])

    #sol=np.linalg.solve(a,b)
    #print(sol[2])

    #finding array of coordinates
    positions = discretized_skin_pos

    shearstress = []

    
    for i in range(len(positions)):
        if x == la:
            s3 = 0
            s2 = 0
            s2a = 0
            s2b = 0
            s1 = 0
            sa = 0 
            sb = 0
        elif x>x3:
            s3 = 0
            s2 = 0
            s2a = 0
            s2b = 0
            s1 = 0
            sa = 0
            
        elif x>x2b and x<x3:
            s2 = 0
            s2a = 0
            s2b = 0
            s1 = 0
            sa = 0

        elif x>x2 and x<x2b:
            s2 = 0
            s2a = 0
            s1 = 0
            sa = 0

        elif x>x2a and x<x2:
            s2b = 0
            s2 = 0
            s3 = 0
            sb = 0

        elif x>x1 and x<x2a:
            s3 = 0
            s2 = 0
            s2a = 0
            s2b = 0
            sb = 0
            
        elif x<x1:
            s3 = 0
            s2 = 0
            s2a = 0
            s2b = 0
            s1 = 0
            sb = 0

        elif x==0:
            s3 = 0
            s2 = 0
            s2a = 0
            s2b = 0
            s1 = 0
            sb = 0
            sa = 0

        Sz = Fz1*s1 + FzI*s2a + Fz2*s2 - P*s2b
        Sy = Fy1*s1 + Fy2*s2 + Fy3*s3 - q*x*sa - q*(x3-x)*sb
        T = -P*(dact2)*s2b + FzI*(dact1)*s2a + Fz1*(d1-yp)*s1 - ((Ca/4)*cos(theta))*q*x*sa - ((Ca/4)*cos(theta))*q*(x3-x)*sb

        k1 = -(Sz * Izz - Sy * Izy) / (Izz * Iyy - Izy ** 2)
        k2 = -(Sy * Iyy - Sz * Izy) / (Izz * Iyy - Izy ** 2)

        s12 = ha
        s23 = c

        # qb as a function of theta/s for each wall, q12_I init for I, q_23 init for II, others q2+q1(end)
        q12_Ib = k1*tsk*(ha/2)**2*np.cos(phi) + k2*tsk*(ha/2)**2*np.sin(phi)

        q21b = k2*tsp*(-ha/2*s12+0.5*s12**2) + q12_Ib*(np.pi/2)

        q23b = k1*tsk*(Ca-ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2 / 2 + \
                        k2*tsk*(-ha/2*s23 + (ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2/2)

        q31b = k1*tsk*((Ca-ha/2)*s23 -(Ca-ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2 / 2) + \
                        k2*tsk*(ha/2)/np.hypot((Ca-ha/2), ha/2)*s23**2/2 + q23b*(np.hypot((Ca-ha/2), ha/2))

        q12_IIb = k2*tsp*(ha/2*s12-s12**2 / 2) + q31b*(np.hypot((Ca-ha/2), ha/2))

        dthetadxc1 = (q12_Ib*(pi*ha/2)/tsk + q21b*ha/tsp)/(2*A1*G)
        dthetadxc2 = (q23b*c/tsk + q31b*c/tsp + q12_IIb*ha/tsp)/(2*A2*G)
        Tvar = -(q23b+q31b)*(ha/2*sin(atan((Ca-ha/2)/(ha/2))))-2*q12_Ib*ha/2

        b = np.matrix([[T+Tvar],
                       [dthetadxc1],
                       [dthetadxc2]])
        sol=np.linalg.solve(A,b)

        #twistsec += sol[2]*dx*-1*180/pi

        q12 = sol[0] + q12_Ib
        q21 = (sol[0] - sol[1]) + (q21b-q12_IIb)
        q23 = sol[1] + q23b
        q31 = sol[1] + q31b
        
        if positions[i][0] > ha/2 and positions[i][1] > 0:
            shearstress.append(float(q31)/tsk)
        #if positions[i][0] = ha/2:
            #shearstress.append(float(q21))
        elif positions[i][0] < ha/2:
            shearstress.append(float(q12)/tsk)
        elif positions[i][0] > ha/2 and positions[i][1] < 0:
            shearstress.append(float(q23)/tsk)
        else:
            shearstress.append(0)
    return shearstress
    #print(sol[1])
    #print(q23)
    #print(twistsec)
    #xlst.append(x)
    #print(float(twistsec))
    #twistseclst.append(float(twistsec))
    #x = x + dx
    

#twistsecrad = twistsec*180/pi
#print(twistseclst)

#plt.plot(xlst, twistseclst)
#plt.axis([0.15, -0.5, -0.5, 0.5])
#plt.show()
