import numpy as np
from plotter import plotter
import os


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


def check_reaction_forces(path_forces_valid, path_forces_sim):
    data_sim = np.genfromtxt(path_forces_sim, delimiter=',')
    data_valid = read_reaction_forces(path_forces_valid)
    print('Y-Axis')
    for i, j in [[0, 0], [2, 2], [3, 4]]:
        sim = data_sim[i, 1]
        val = data_valid[j, 3]*1e3
        print('{:.3f}\t-\t{:.3f}\t=\t{:.3f}\t = {}%'.format(sim/10**3, val/10**3, (sim-val)/10**3, round(abs((sim-val)/val)*100, 2)))
    print('Z-Axis')
    for i, j in [[0, 0], [1, 1], [2, 2]]:
        sim = data_sim[i, 2]
        val = data_valid[j, 4]*1e3
        print('{:.3f}\t-\t{:.3f}\t=\t{:.3f}\t = {}%'.format(sim/10**3, val/10**3, (sim-val)/10**3, round(abs((sim-val)/val)*100, 2)))


def check_deflections(path_node, node_U, path_deflection, name):
    # idea is to compare deflections with closes node
    defl_data = np.genfromtxt(path_deflection, delimiter=',')
    node_data = np.genfromtxt(path_node, comments='*', delimiter=',', max_rows=3205 - 9)
    defl_error = []
    for point in defl_data:
        closest_point = node_data[np.argmin(np.sum((node_data[:,1:] - point[:3]*1e3)**2, axis=1)), :]
        defl_error.append((np.abs(point[3:]*1e3 -node_U[int(closest_point[0]-1)][-2:])/node_U[int(closest_point[0]-1)][-2:]*100))
    import matplotlib.pyplot as plt
    fig = plt.figure()
    plt.hist(list(filter(lambda x: 0 <= x < 200, np.array(defl_error)[:, 0])), bins=20, label=['y displacement'], histtype='step', fill=False, linewidth=2.5)
    plt.hist(list(filter(lambda x: 0 <= x < 200, np.array(defl_error)[:, 1])), bins=20, label=['z displacement'], histtype='step', fill=False, linewidth=2.5)
    plt.legend()
    # plt.hist(list(filter(lambda x: 0 <= x < 200, np.array(defl_error)[:, 1])), bins=20)
    plt.xlabel('Relative error [%]')
    plt.ylabel('Number of points')
    # plt.show()
    # plt.show()
    plt.savefig('defl_error_{}.png'.format(name))
    return defl_error


def check_stress(path_node, node_S, path_stress):
    stress_data = np.genfromtxt(path_stress, delimiter=',')
    node_data = np.genfromtxt(path_node, comments='*', delimiter=',', max_rows=3205 - 9)
    stress_error = []
    for point in stress_data:
        closest_point = node_data[np.argmin(np.sum((node_data[:,1:] - point[:3]*1e3)**2, axis=1)), :]
        stress_data.append((np.abs(point[3:]*1e3 -node_S[int(closest_point[0]-1)][-2:])/node_S[int(closest_point[0]-1)][1]*100))
    return stress_error

cases = ['R1', 'R2', 'LC1', 'LC2']
path_node = '../data/CRJ700n.inp'
node_data = read_node_loc(path_node)

# plotter(node_data, node_stress=None, node_U=read_node_U('../data/CRJ700n_ULC1.rpt'), show_stress=False)

for case in cases:
    path_forces_sim = '../data/reaction_force_{}.csv'.format(case)
    path_deflection = '../data/deflection_{}.csv'.format(case)

    if 'reaction_force_{}.csv'.format(case) in os.listdir('../data'):
        path_U = '../data/CRJ700n_U{}.rpt'.format(case)
        path_stress = '../data/CRJ700n_S{}.rpt'.format(case)
        node_U = read_node_U(path_U)
        node_S = read_node_stress(path_stress)
        unique, counts = np.unique(node_S[:, 0], return_counts=True)
        node_S_upd = np.unique(np.array(list(map(lambda x: x if counts[int(x[0]) - 1] == 1 else np.array([x[0], 0]), node_S))), axis=0)

        print('Case : {}'.format(case))
        check_reaction_forces(path_U, path_forces_sim)
        check_deflections(path_node, node_U, path_deflection, case)
        # check_stress(path_node, node_S, path_stress)
        print('----------------------------\n')
    else:
        print('----------------------------')
        print('Case {} not generated yet'.format(case))
        print('----------------------------')

