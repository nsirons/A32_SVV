import plotly
import plotly.figure_factory as FF

import numpy as np
from scipy.spatial import Delaunay
from scipy.interpolate import griddata


def gen_mesh(x,y):
    u = np.arange(0, x, 1)
    v = np.arange(0, y, 1)
    u, v = np.meshgrid(u, v)
    u = u.flatten()
    v = v.flatten()
    return u,v


def plot(node_data, node_stress, dx, ds):

   #+1 from point at TE
    ds = ds+2
    u = np.arange(0, ds, 1)
    v = np.arange(0, dx, 1)

    #u along cross-section, v along length
    u,v = np.meshgrid(u,v)

    u = u.flatten()
    v = v.flatten()

    x=[]
    y=[]
    z=[]

    xxx = u + v

    print(node_data)

    for i, j in zip(u,v):
        xx = i + j
        xx = node_data[j*ds + i][1]*1e3
        yy = node_data[j*ds + i][2]*1e3
        zz = node_data[j*ds + i][3]*1e3
        x.append(xx)
        y.append(yy)
        z.append(zz)


    points2D = np.vstack([u,v]).T
    tri = Delaunay(points2D)
    simplices = tri.simplices

    nodes = len(node_data)
    stresses = node_stress[:,1]
    xs = np.array([node_data[:,1]]).transpose()
    ys =  np.array([node_data[:,2]]).transpose()
    zs =  np.array([node_data[:,3]]).transpose()
    
    
    

    def f(x,y,z):
        print(x/(1e3),y/(1e3),z/(1e3))
        xq, yq, zq = np.meshgrid(x/(1e3),y/(1e3),z/(1e3))
        vq = griddata(np.concatenate((xs,ys,zs), axis=1),stresses, (xq, yq, zq))
        print(vq)
        return vq

    #f = RegularGridInterpolator((node_data[1], node_data[2], node_data[3]),stresses) 

  

    camera = dict(
        up=dict(x=0, y=1, z=0),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1, y=0.1, z=0.1)
    )
    fig1 = FF.create_trisurf(x=x, y=y, z=z,
                                simplices=simplices,
                                title="aileron",
                                colormap="Portland",
                                show_colorbar=True,
                                color_func=f,
                                aspectratio=dict(x=1, y=(max(y)-min(y)) / 1691, z=(max(z)-min(z)) / 1691),
                                )

    import plotly.graph_objs as go
    trace = go.Scatter3d()
    # fig1['layout'].update(
    #     scene=dict(camera=camera))
    plotly.offline.plot(fig1, filename="Aileron")
    #plotly.offline.plot([fig1.data[0],fig2.data[0], trace], filename="Aileron.html",)

def plotter(node_data,*, node_stress=None, node_U=None, show_stress=True):
    # TODO : theta as function of x
    # So ugly, but it works
    theta = 26*np.pi/180
    deformation = 0 if node_U is None else 1.5
    rot_mat = np.array([[-np.cos(theta), np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    # coord_to_node = dict([[(i[1], i[2], i[3]), node_U[int(i[0])-1, 1]] for i in node_data])

    if node_stress is None and node_U is not None:

        coord_to_node_rot = dict([[(i[1],
                                np.sin(theta)*(i[3]+node_U[int(i[0]) - 1, 4]*deformation)+np.cos(theta)*(i[2]+node_U[int(i[0]) - 1, 3]*deformation),
                                -np.cos(theta) * (i[3]+node_U[int(i[0]) - 1, 4]*deformation) + np.sin(theta) * (i[2]+node_U[int(i[0]) - 1, 3]*deformation)), node_U[int(i[0]) - 1, 1]] for i in node_data])

    elif node_stress and node_U is None:
        coord_to_node_rot = dict([[(i[1],
                                np.sin(theta)*i[3]+np.cos(theta)*i[2],
                                -np.cos(theta) * i[3] + np.sin(theta) * i[2]), node_stress[int(i[0]) - 1, 1]] for i in node_data])

    elif node_U and node_stress and show_stress:
        pass
    elif node_U and node_stress and not show_stress:
        pass
    else:
        raise AttributeError('stress and displacement is not provided')
    print(max(coord_to_node_rot.values()), min(coord_to_node_rot.values()))
    n_per_section = 47
    x_sections = 68

    # for looping through all elements
    u, v = gen_mesh(n_per_section + 1, x_sections)

    # for plotting aileron skin
    uu, vv = gen_mesh(43, x_sections)

    # for plotting spar
    uuu, vvv = gen_mesh(7, x_sections)

    x = []
    y = []
    z = []

    xxx = []
    yyy = []
    zzz = []
    fix = 0  # fix is required to adjust for missing nodes
    node_data.sort(order=["x", "z", "y"])
    print(node_data[:47][:])
    for i, j in zip(u, v):
        # start from one side, and continue in backward order
        if i < int(n_per_section / 2) + 1:
            if node_data[j * n_per_section + 2 * i - fix][3] != 0 or \
                    abs(node_data[j * n_per_section + 2 * i - fix][2]) == 86.5:
                x.append(node_data[j * n_per_section + 2 * i - fix][1])
                yy = node_data[j * n_per_section + 2 * i - fix][2] + float(node_U[int(node_data[j * n_per_section + 2 * i - fix][0])-1, 3]*deformation)
                zz = node_data[j * n_per_section + 2 * i - fix][3] + float(node_U[int(node_data[j * n_per_section + 2 * i - fix][0])-1, 4]*deformation)
                zz,yy = rot_mat @ np.array([zz, yy])
                y.append(yy)
                z.append(zz)
            else:
                fix = 1

        elif i != 47:
            if node_data[j * n_per_section + n_per_section * 2 - 2 * i + fix - 1][3] != 0 or \
                    abs(node_data[j * n_per_section + n_per_section * 2 - 2 * i + fix - 1][2]) == 86.5:
                x.append(node_data[j * n_per_section][1])
                yy = node_data[j * n_per_section + n_per_section * 2 - 2 * i + fix - 1][2] + float(node_U[int(node_data[j * n_per_section + n_per_section * 2 - 2 * i + fix - 1][0])-1, 3]*deformation)
                zz = node_data[j * n_per_section + n_per_section * 2 - 2 * i + fix - 1][3] + float(node_U[int(node_data[j * n_per_section + n_per_section * 2 - 2 * i + fix - 1][0])-1, 4]*deformation)
                zz, yy = rot_mat @ np.array([zz, yy])
                y.append(yy)
                z.append(zz)
            else:
                fix = 0
        else:
            x.append(node_data[j * n_per_section][1])
            yy = node_data[j * n_per_section][2] + float(node_U[int(node_data[j * n_per_section][0]-1), 3]*deformation)
            zz = node_data[j * n_per_section][3] + float(node_U[int(node_data[j * n_per_section][0]-1), 4]*deformation)
            zz, yy = rot_mat @ np.array([zz, yy])
            y.append(yy)
            z.append(zz)

    for i,j in zip(u,v):
        if i != 47:
            if node_data[j * n_per_section + i][3] == 0 and True:
                    # abs(node_data[j * n_per_section + i][2]) != 86.5:
                    xxx.append(node_data[j * n_per_section + i][1])
                    yy = node_data[j * n_per_section + i][2] + float(node_U[int(node_data[j * n_per_section + i][0] - 1), 3] * deformation)
                    zz = node_data[j * n_per_section + i][3] + float(node_U[int(node_data[j * n_per_section + i][0] - 1), 4] * deformation)
                    zz, yy = rot_mat @ np.array([zz, yy])
                    yyy.append(yy)
                    zzz.append(zz)

    points2D_skin = np.vstack([uu, vv]).T
    points2D_spar = np.vstack([uuu, vvv]).T
    tri_skin = Delaunay(points2D_skin)
    tri_spar = Delaunay(points2D_spar)
    simplices_skin = tri_skin.simplices
    simplices_spar = tri_spar.simplices
    import random
    def stress_color(x, y, z):

        if (x, y, z) in coord_to_node_rot:
            # print(x, y, z)
            return coord_to_node_rot[(x, y, z)]

        return random.random()  # hopefully, it would not

    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1, y=0.1, z=0.1)
    )
    fig1 = FF.create_trisurf(x=x, y=y, z=z,
                             simplices=simplices_skin,
                             title="aileron",
                             colormap="Portland",
                             show_colorbar=True,
                             color_func=stress_color,
                             aspectratio=dict(x=1, y=(max(y)-min(y)) / 1691, z=(max(z)-min(z)) / 1691),
                             )

    fig2 = FF.create_trisurf(x=xxx, y=yyy, z=zzz,
                             simplices=simplices_spar,
                             # title="aileron",
                             # colormap="Portland",
                             show_colorbar=True,
                             color_func=stress_color,
                             aspectratio=dict(x=1, y=(max(y)-min(y)) / 1691, z=(max(z)-min(z)) / 1691),
                             )

    print(simplices_skin)

    import plotly.graph_objs as go
    trace = go.Scatter3d()
    # fig1['layout'].update(
    #     scene=dict(camera=camera))
    plotly.offline.plot(fig1, filename="Aileron")
    #plotly.offline.plot([fig1.data[0],fig2.data[0], trace], filename="Aileron.html",)

