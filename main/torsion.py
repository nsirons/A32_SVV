import numpy as np

# --------------------------------------TORSION------------------------------------------------


def torsion(x, n, Ca, la, ha, theta, x1, x2, x3, xa, d1, d3, dact1, dact2, tsp, tsk,  G, Fy1, Fy2, Fy3, Fx1, Fx3, Fz1, FzI, P, q):
    x2a = x2 - xa / 2  # x-location actuator 1
    x2b = x2 + xa / 2  # x-location actuator 2

    # discretizing
    # n = 20  # number of sections
    dx = la/(n+1)  # discretization step along length aileron
    pos = 0

    #parameter simplifying formula
    b = np.sqrt((Ca-ha/2)**2+(ha/2)**2)
    yp = ha/2*np.cos(theta)

    #cell areas
    A1 = np.pi*ha**2/8
    A2 = ha*(Ca-ha/2)/2

    #Steps (switching terms on/off)
    s1 = 1
    s2 = 1
    s2a = 1
    s2b = 1
    s3 = 1

    twistsec = 0

    twistseclst = []
    poslst = []

    #torsion
    T = - P*(yp+dact2)*s2b + FzI*(yp+dact1)*s2a + Fz1*d1*s1 - ((Ca/4-ha/2)*np.cos(theta))*q*x

    a = np.matrix([[2*A1, 2*A2, 0],
                  [-ha/(G*tsp*2*A2), (2*b/tsk+ha/tsp)/(G*2*A2), -1],
                  [(ha/tsp+np.pi*ha/(4*tsk))/(2*G*A1), -ha/(2*G*tsp*A1), -1]])

    b = np.matrix([[T],
                  [0],
                  [0]])

    sol=np.linalg.solve(a,b)
    #print(sol[2])

    while pos < x:

       if pos>x3:
           s3 = 0
           s2 = 0
           s2a = 0
           s2b = 0
           s1 = 0

       elif pos>x2b and pos<x3:
           s2 = 0
           s2a = 0
           s2b = 0
           s1 = 0

       elif pos>x2 and pos<x2b:
           s2 = 0
           s2a = 0
           s1 = 0

       elif pos>x2a and pos<x2:
           s2a = 0
           s1 = 0

       elif pos>x1 and pos<x2a:
           s1 = 0

       T = -P*(yp+dact2)*s2b + FzI*(yp+dact1)*s2a + Fz1*d1*s1 - ((Ca/4-ha/2)*np.cos(theta))*q*dx

       b = np.matrix([[T],
                      [0],
                      [0]])

       sol=np.linalg.solve(a,b)

       twistsec += sol[2]*dx*-1*180/np.pi
       print(pos)
       print(sol[2])
       poslst.append(pos)
       print(float(twistsec))
       twistseclst.append(float(twistsec))
       pos = pos + dx

    #print(poslst)
    twistsecrad = twistsec*180/pi
    #print(twistseclst)

    # plt.plot(poslst, twistseclst)
    #plt.axis([0.15, -0.5, -0.5, 0.5])
    # plt.show()
