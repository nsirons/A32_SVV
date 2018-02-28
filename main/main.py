# constants/geometry
from reaction_forces import reaction_forces
from bendingstresses import find_bending_stresses
from shearstressestorshe import find_shear_stresses
from forcesmomentsfunc import int_shear, int_moment
import numpy as np
import logging
import sys
import argparse
from math import *
import matplotlib.pyplot as plt
from aileron import aileron
from MOI import calculate_inertia_yy, calculate_inertia_zz, calculate_rotated_inertia
from tools import *

# -- Geometry --
C_a = 0.484  # chord length aileron
l_a = 1.691  # span aileron
x_1 = 0.149  # x-location hinge 1
x_2 = 0.554  # x-location hinge 2
x_3 = 1.541  # x-location hinge 3
x_a = 27.2e-2  # distance actuator 1 - 2
h_a = 17.3e-2  # aileron height
t_sk = 1.1e-3  
t_sp = 2.5e-3  
t_st = 1.2e-3
h_st = 1.4e-2  # Height of stiffener
w_st = 1.8e-2  # Width of stiffener
n_st = 13  # number of stringers
d_1 = 6.81e-2/2.54  # ver. displacement hinge 1
d_3 = 20.3e-2/2.54  # ver. displacement hinge 3
d_act1 = 0  # ? # iterative*
d_act2 = 0  # ? # iterative*
theta = 26 * np.pi/180  # Maximum upward deflection

# -- Loads --
P = 37.9e3  # load in actuator 2
q = 27.1e3  # net aerodynamic load

# -- Material Property
# source: http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t3
E = 73.1e9  # Young modulus aluminium 2024-T3
G = 28.e9  # Shear modulus aluminium 2024-T3

# Calculate normal stress along the cross-sectional area and x direction
d = 500 #discretization along span
x = 0 #initial x-coordinate
dx = l_a/(d+1) #steps along span
n = 3  # number of discretized points
x_points = []
y_points = []
z_points = []
sigma_points = []
threedsigmamax = []

#stresses
sigmay = 0
sigmaz = 0
tauxy = 0
tauzx = 0


logger = logging.getLogger()
aileron_obj = aileron(C_a, h_a, t_sk, t_sp, t_st,h_st,w_st,n_st,theta)


def get_section_properties(aileron_obj):


    inertia_uu = calculate_inertia_zz(aileron_obj)
    inertia_vv = calculate_inertia_yy(aileron_obj)
    
    Ivec = calculate_rotated_inertia(inertia_uu, inertia_vv, 0, -degrees(aileron_obj.angle))

    I_zz = Ivec[0]
    I_yy = Ivec[1]
    I_zy = Ivec[2]
    ybar = 0#None  # change to functions
    zbar = 0.189#None  # change to functions

    return I_zz, I_yy, I_zy, ybar, zbar


def get_reaction_forces(I_zz, I_yy, I_zy):
    x = reaction_forces(C_a, l_a, h_a, x_a, x_1, x_2, x_3, d_1, d_3, d_act1, d_act2, degrees(theta), E, P, q, I_zz, I_yy, I_zy)
    reaction_forces_dict = {
        'Fy1': x[0][0],
        'Fz1': x[1][0],
        'Fy2': x[2][0],
        'Fz2': x[3][0],
        'Fy3': x[4][0],
        'FzI': x[5][0]
        }
    return reaction_forces_dict


def get_moment_functions(reaction_forces_dict):
    sf_y = int_shear([x_1, x_2, x_3, l_a], [reaction_forces_dict['Fy1'], reaction_forces_dict['Fy2'], reaction_forces_dict['Fy3'], 0], -q)
    sf_z = int_shear([x_1, (x_2-x_a/2), x_2, (x_2+x_a/2), l_a], [reaction_forces_dict['Fz1'], reaction_forces_dict['FzI'], reaction_forces_dict['Fz2'], -P, 0], -q)
    m_z = int_moment(sf_z)
    m_y = int_moment(sf_y)
    return (sf_y, sf_z, m_y, m_z)
    


def get_von_misses(sigmaz, tauyz):
    sigmamax = []
    sigmay = 0
    sigmax = 0
    tauzx = 0
    tauxy = 0
    for i in range(len(sigmaz)):
        sigmamax.append(sqrt(1/2*((sigmax-sigmay)**2 + (sigmay-sigmaz[i])**2 + (sigmaz[i]-sigmax)**2)) + sqrt(3*(tauxy**2 + tauyz[i]**2 + tauzx**2)))
    return sigmamax

def write_header(fp):
    fp.write("-------------------------------------------------\n")
    fp.write("|                 SVV A32 Results               |\n")
    fp.write("-------------------------------------------------\n")


def write_program_settings(fp):
    fp.write("-------------------------------------------------\n")
    fp.write("|          Discretization Parameters            |\n")
    fp.write("| Amount of segments x-dir: " + str(d)+ "\t\t\t|\n")
    fp.write("| Step Size x-dir: " + str(dx)+ "\t|\n")
    fp.write("| Number of discretization points along skin: " + str(n) + "\t|\n")
    fp.write("-------------------------------------------------\n")


def write_section_properties(fp, Izz, Iyy, Izy):
    fp.write("-------------------------------------------------\n")
    fp.write("|          Discretization Parameters            |\n")
    fp.write("| Inertia zz: " + str(Izz)+ "\t\t|\n")
    fp.write("| Inertia yy: " + str(Iyy)+ "\t\t|\n")
    fp.write("| Inertia zy: " + str(Izy) + "\t\t|\n")
    fp.write("-------------------------------------------------\n")

def write_forces(fp, reaction_forces_dict):
    fp.write("-------------------------------------------------\n")
    fp.write("|               Reaction Forces [N]             |\n")
    fp.write("| Fy1: " + str(int(reaction_forces_dict['Fy1'])) + "\t\t\t Fz1: " + str(int(reaction_forces_dict['Fz1'])) + "\t|\n")
    fp.write("| Fy2: " + str(int(reaction_forces_dict['Fy2'])) + "\t\t\t Fz2: " + str(int(reaction_forces_dict['Fz2'])) + "\t|\n")
    fp.write("| Fy3: " + str(int(reaction_forces_dict['Fy3'])) + "\t\t\t FzI: " + str(int(reaction_forces_dict['FzI'])) + "\t|\n")
    fp.write("-------------------------------------------------\n")


def main(args):
    
    parser = argparse.ArgumentParser(description="SVV Program 2018 A32")
    parser.add_argument("-o", "--out", help="Specify output file", default=sys.stdout, type=argparse.FileType('w'))
    parser.add_argument("-p", "--plots", help="Show graphical plots", default=0)
    parser.add_argument("-v", "--verbose", help="Give more updates during runnning", action='count', default=0)
    parser.add_argument("-rg", "--reaction_forces_generate", help="Generates a JSON file containing manual entered reaction forces", type=argparse.FileType('w'))
    parser.add_argument("-r", "--reaction_forces", help="Reads a JSON file containing manual entered reaction forces", type=argparse.FileType('r'))
    
    arguments = parser.parse_args(args)
    logger.setLevel(arguments.verbose*10)

    context = (logger, aileron_obj)

    #Inertia/cross section properties
    I_zz, I_yy, I_zy, ybar, zbar = get_section_properties(context[1])
    
    #Reaction forces
    reaction_forces_dict = {}

    if(arguments.reaction_forces_generate != None):
        generate_reaction_forces_file(arguments.reaction_forces_generate)
        print("Reaction Forces file written")
    elif(arguments.reaction_forces != None):
        reaction_forces_dict = read_reaction_forces(arguments.reaction_forces)
    else:
        reaction_forces_dict = get_reaction_forces(I_zz, I_yy, I_zy)

    F_y1, F_y2, F_y3, F_z1, F_z2, F_zI = convert_reaction_forces_dict(reaction_forces_dict)

    F_y, F_z, M_y, M_z = get_moment_functions(reaction_forces_dict)

    ##plt.plot(np.linspace(0,l_a, 50), F_y(np.linspace(0,l_a,50)))
    ##plt.plot(np.linspace(0,l_a, 50), F_z(np.linspace(0,l_a,50)))
    ##plt.plot(np.linspace(0,l_a, 50), M_y(np.linspace(0,l_a,50)))
    ##plt.plot(np.linspace(0,l_a, 50), M_z(np.linspace(0,l_a,50)))
    #plt.show()


    #start loop over length

    x_lst = []
    y_lst = []
    z_lst = []
    sigma_z_lst = []
    tau_yz_lst = []
    sigma_max_lst = []

    for current_distance in np.arange(0, l_a, dx):
        print(current_distance)

        #Bending stresses
        x_pos, y_pos, z_pos, sigma_z = find_bending_stresses(current_distance, n, l_a, x_1, x_2, x_3, x_a, d_1, d_3, I_zz, I_yy, I_zy, ybar, zbar, M_y, M_z)

        #Shear stresses
        tau_yz = find_shear_stresses(current_distance, n, l_a, x_1, x_2, x_3, x_a, d_1, d_3, C_a, h_a, G, t_sp, t_sk, d_act1, d_act2, I_zz, I_yy, I_zy, ybar, zbar, theta, F_z2, F_y1, F_y2, F_y3, 0, 0, F_z1, F_zI, P, q)
        
        #Von Misses
        sigma_max = get_von_misses(sigma_z, tau_yz)        

        #Deflection
        #deflection_z, deflection_y = get_aileron_deflections()


        #Store local section results
        x_lst.append(x_pos)
        y_lst.append(y_pos)
        z_lst.append(z_pos)
        sigma_z_lst.append(sigma_z)
        tau_yz_lst.append(tau_yz)
        sigma_max_lst.append(sigma_max)


    plt.plot([x[0] for x in x_lst], [s[0] for s in sigma_max_lst])
    plt.show()

    #Write Output
    write_header(arguments.out)
    write_program_settings(arguments.out)
    write_section_properties(arguments.out, I_zz, I_yy, I_zy)
    write_forces(arguments.out, reaction_forces_dict)

   

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))