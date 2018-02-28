from math import *
from aileron import aileron




def get_deflection(x_position, aileron, Izz, Iyy, Izy):

    moment_y = get_moment_y(x_position)
    moment_z = get_moment_z(x_positon)

    deflection_y = get_deflection_y(x_position, moment_y, moment_z, Izz, Iyy, Izy)
    deflection_z = get_deflection_z(x_position, moment_y, moment_z, Izz, Iyy, Izy)

    return (deflection_z, deflection_y)

def get_deflection_y(x_position, moment_y, moment_z, Izz, Iyy, Izy):

    #TODO: correct for deflection in hinge 2
    y = -( ( - moment_z * Iyy - moment_y * Izy) / (E * (Izz*Iyy - Izy**2)) )
  
    return y


def get_deflection_z(x_position, moment_y, moment_z, Izz, Iyy, Izy):
    
    #TODO: correct for deflection in hinge 2
    z = - ( (-moment_z * Izy + moment_y * Izz) / (E * (Izz*Iyy - Izy**2)) )  

    return z




