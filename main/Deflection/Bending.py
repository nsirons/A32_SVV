from math import *
from aileron import aileron



def get_moment_y(x):
    return 0

def get_moment_z(x):
    return 0



def get_deflection(x_position, aileron):

    moment_y = get_moment_y(x_position)
    moment_z = get_moment_z(x_positon)

    deflection_y = get_deflection_y(x_position, moment_y, moment_z, aileron)
    deflection_z = get_deflection_z(x_position, moment_y, moment_z, aileron)

    return (deflection_y, deflection_z)

def get_deflection_y(x_position, moment_y, moment_z, aileron):

    #TODO: correct for deflection in hinge 2
    y = -( ( - moment_z * aileron.Iyy - moment_y * aileron.Izy) / (E * (aileron.Izz*aileron.Iyy - aileron.Izy**2)) )
  
    return y


def get_deflection_z(x_position, moment_y, moment_z, aileron):
    
    #TODO: correct for deflection in hinge 2
    z = - ( (-moment_z * aileron.Izy + moment_y * aileron.Izz) / (E * (aileron.Izz*aileron.Iyy - aileron.Izy**2)) )  

    return z




