import json

def generate_reaction_forces_file(fp):
    reaction_forces = {
        'Fy1': 0,
        'Fz1': 0,
        'Fy2': 0,
        'Fz2': 0,
        'Fy3': 0,
        'FzI': 0
        }

    json.dump(reaction_forces, fp)
    return


def read_reaction_forces(fp):
    return json.load(fp)

def convert_reaction_forces_dict(dict):
    F_y1 = dict['Fy1']
    F_y2 = dict['Fy2']
    F_y3 = dict['Fy3']
    F_z1 = dict['Fz1']
    F_z2 = dict['Fz2']
    F_zI = dict['FzI']

    return F_y1, F_y2, F_y3, F_z1, F_z2, F_zI
