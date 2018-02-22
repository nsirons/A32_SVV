import matplotlib.pyplot as plt
from math import *
import numpy as np

# x <= la:

#x = float(input("enter a x-coordinate: "))

def find_bending_stresses(x,n):
    #input variables
    x1 = 0.149 #x-location hinge 1
    x2 = 0.554 #x-location hinge 2
    x3 = 1.541 #x-location hinge 3
    xa = 0.272 #distance actuator 1 - 2
    x2a = x2-xa/2 #x-location actuator 1
    x2b = x2+xa/2 #x-location actuator 2
    d1 = 0.0681 #ver. displacement hinge 1
    d3 = 0.203 #ver. displacement hinge 3
    la = 1.691 #span aileron
    q = 27100 #net aerodynamic load
    P = 37900 #load in actuator 2
    Ca = 0.84 #chord length aileron
    ha = 0.173 #aileron height
    G = 28E+9 #shear modulus aluminium 2024-T3
    tsp = 0.0025 #spar thickness
    tsk = 0.0011 #skin thickness
    dact1 = 0 #iterative*
    dact2 = 0 #iterative*
    theta = 26/180*pi #max. upward deflection
    Fy1 = 10000#from reaction force code*
    Fy2 = 10000#from reaction force code*
    Fy3 = 10000#from reaction force code*
    FzI = 10000#from reaction force code*
    Fz1 = 10000#from reaction force code*
    Fx1 = 10000#from reaction force code*
    Fx3 = 10000#from reaction force code*

    #--------------------------------------BENDING------------------------------------------------

    #MOIS
    Izy = 0
    Izz = 0.0002#output MOI code*
    Iyy = 0.0002#output MOI code*

    #centroid cross-section
    zbar = 0.15 #output MOI code*
    ybar = 0 

    #lists
    sigmaxlst = []
    zlst = []
    Mylst = []
    Mzlst = []

    #Steps (switching terms on/off)
    s1 = 1
    s2 = 1
    s2a = 1
    s2b = 1
    s3 = 1

    #x = 0
    #dx = 0.1 #discretization along length aileron
    #n = 50 #discretization cross-section

    #finding array of coordinates
    def calculate_stringer_positions(stiffener):
        height = 0.173
        chord = 0.484
        length_flat_skin = sqrt( (height/2)**2 + (chord - height/2)**2 )
        length_circular_skin = pi*height/2
        angle_flat_skin = atan2(height/2, chord - height/2)

        length_total_skin = 2*length_flat_skin + length_circular_skin
        distance_between_stringers = length_total_skin / (stiffener + 1)

        stringer_positions = []
        current_length = 0

    for i in range(1, stiffener+1):
        current_length += distance_between_stringers
        position = []
        if  0 <= current_length <= length_flat_skin:
            x = chord - current_length*cos(angle_flat_skin)
            y = current_length*sin(angle_flat_skin)
            position = [x,y]
        elif  length_flat_skin < current_length <= (length_flat_skin + length_circular_skin):
            current_angle = pi/2 + (current_length - length_flat_skin)/(2*pi*height/2)*2*pi
            x = height/2 + height/2*cos(current_angle)
            y = height/2*sin(current_angle)
            position = [x,y]
        elif (length_flat_skin + length_circular_skin) < current_length <= length_total_skin:
            modified_current_length = length_total_skin - current_length
            x = chord - modified_current_length*cos(angle_flat_skin)
            y = modified_current_length*sin(-angle_flat_skin)
            position = [x,y]
        else:
            position = [-1, -1]
        stringer_positions.append(position)
    return stringer_positions

    positions = calculate_stringer_positions(n)
    
    xlst = []
    ylst = []
    zlst = []
    sigmaxlst = []
    for i in range(len(positions)-1):

        z = zbar - positions[i][0] 
        y = positions[i][1]
        if x>x3:
            s3 = 0
            s2 = 0
            s2a = 0
            s2b = 0
            s1 = 0

        elif x>x2b and x<x3:
            s2 = 0
            s2a = 0
            s2b = 0
            s1 = 0

        elif x>x2 and x<x2b:
            s2 = 0
            s2a = 0
            s1 = 0

        elif x>x2a and x<x2:
            s2a = 0
            s1 = 0
            
        elif x>x1 and x<x2a:
            s1 = 0

        Mz = Fy3*(x3-x)*s3 + Fy2*(x2-x)*s2 + Fy1*(x1-x)*s1 - Fx3*d3*s3\
             - Fx1*d1*s1 - ((la/2)-x2)*q*la
        My = Fz1*(x1-x)*s1 + P*(x2b-x)*s2b + FzI*(x2a-x)*s2a
        #Mzlst.append(Mz)
        #Mylst.append(My)
        
        sigmax =((Mz*Iyy-My*Izy)*y+(My*Izz-Mz*Izy)*z)/(Izz*Iyy-Izy**2)
        #print(z,y,sigmax)
        sigmaxlst.append(sigmax)
        xlst.append(x)
        ylst.append(y)
        zlst.append(z)
    return xlst,ylst,zlst,sigmaxlst
    #x = x + dx

#print(sigmaxlst)
#print(sigmaxlst)
#print(zlst)
#print(len{zlst))
#plt.scatter(sigmaxlst, zlst)
#plt.axis([0.15, -0.5, -0.5, 0.5])
#plt.plot([x for x,y in positions], [y for x,y in positions])
#plt.show()

#--------------------------------------TORSION------------------------------------------------

#discretizing
##n = 20 #number of sections
##dx = la/(n+1) #discretization step along length aileron
##pos = 0
##
###parameter simplifying formula
##b = sqrt((Ca-ha/2)**2+(ha/2)**2)
##yp = ha/2*cos(theta)
##
###cell areas
##A1 = pi*ha**2/8
##A2 = ha*(Ca-ha/2)/2
##
###Steps (switching terms on/off)
##s1 = 1
##s2 = 1
##s2a = 1
##s2b = 1
##s3 = 1
##
##twistsec = 0
##
##twistseclst = []
##poslst = []
##
###torsion
##T = - P*(yp+dact2)*s2b + FzI*(yp+dact1)*s2a + Fz1*d1*s1 - ((Ca/4-ha/2)*cos(theta))*q*x
##
##a = np.matrix([[2*A1, 2*A2, 0],
##               [-ha/(G*tsp*2*A2), (2*b/tsk+ha/tsp)/(G*2*A2), -1],
##               [(ha/tsp+pi*ha/(4*tsk))/(2*G*A1), -ha/(2*G*tsp*A1), -1]])
##
##b = np.matrix([[T],
##               [0],
##               [0]])
##
##sol=np.linalg.solve(a,b)
###print(sol[2])
##
##while pos < x:
##    
##    if pos>x3:
##        s3 = 0
##        s2 = 0
##        s2a = 0
##        s2b = 0
##        s1 = 0
##
##    elif pos>x2b and pos<x3:
##        s2 = 0
##        s2a = 0
##        s2b = 0
##        s1 = 0
##
##    elif pos>x2 and pos<x2b:
##        s2 = 0
##        s2a = 0
##        s1 = 0
##
##    elif pos>x2a and pos<x2:
##        s2a = 0
##        s1 = 0
##        
##    elif pos>x1 and pos<x2a:
##        s1 = 0
##
##    T = -P*(yp+dact2)*s2b + FzI*(yp+dact1)*s2a + Fz1*d1*s1 - ((Ca/4-ha/2)*cos(theta))*q*dx
##
##    b = np.matrix([[T],
##                   [0],
##                   [0]])
##
##    sol=np.linalg.solve(a,b)
##
##    twistsec += sol[2]*dx*-1*180/pi
##    print(pos)
##    print(sol[2])
##    poslst.append(pos)
##    print(float(twistsec))
##    twistseclst.append(float(twistsec))
##    pos = pos + dx
##    
###print(poslst)
##twistsecrad = twistsec*180/pi
###print(twistseclst)
##
##plt.plot(poslst, twistseclst)
###plt.axis([0.15, -0.5, -0.5, 0.5])
##plt.show()




