import numpy as np
from plotter import plotter

def read_node_loc(path):
    # Element label - X - Y - Z
    return np.genfromtxt(path, comments='*', delimiter=',', max_rows=3205 - 9, names=["index", "x", "y", "z"])

def read_reaction_forces(path):
    return np.genfromtxt(path, skip_header=3231, max_rows=3247-3231)

def read_node_U(path):
    # node label - U Magnitude - Ux - Uy - Uz
    return np.genfromtxt(path, skip_header=19, usecols=(0, 5, 6, 7, 8), max_rows=3205 - 9)


def read_node_stress(path):
    # for first part?
    # TODO : check with other files
    region_1 = np.unique(np.genfromtxt(path, usecols=(1, 2), skip_header=23, max_rows=2703 - 23), axis=0)
    region_2 = np.unique(np.genfromtxt(path, usecols=(1, 2), skip_header=2724, max_rows=4332-2724), axis=0)
    region_3 = np.unique(np.genfromtxt(path, usecols=(1, 2), skip_header=4353, max_rows=8641 - 4353), axis=0)
    region_4 = np.unique(np.genfromtxt(path, usecols=(1, 2), skip_header=8662, max_rows=12950 - 8662), axis=0)
    return np.unique(np.vstack((region_1, region_2, region_3, region_4)), axis=0)

    # return np.genfromtxt(path, skip_header=23, max_rows=3205 - 9)


path_node = '../data/CRJ700n.inp'
path_U = '../data/CRJ700n_UR1.rpt'
path_stress = '../data/CRJ700n_SR2.rpt'

node_data = read_node_loc(path_node)
node_U = read_node_U(path_U)
node_S = read_node_stress(path_stress)

unique, counts = np.unique(node_S[:,0], return_counts=True)
node_S_upd = np.unique(np.array(list(map(lambda x: x if counts[int(x[0])-1] == 1 else np.array([x[0], 0]), node_S))), axis=0)
plotter(node_data, node_stress=node_S_upd, node_U=node_U, show_stress=True)


def check_reaction_forces():
    pass

def diff(node_data, node_U):
    # run main for specific points
    simulation_U = []
    for node in node_data:
        # calc_deflection(node coordinates)
        simulation_U.append(None)
    # plot it to see
    return np.array(simulation_U) - node_U[:, 1:]
