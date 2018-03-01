from scipy.interpolate import interp1d
from itertools import chain
import scipy.integrate as integrate
import numpy as np


def int_shear(x, y, q=0):
    dx = 0.00000000000001  # small offset, for step jump

    x_shear = []
    y_shear = []
    for i in range(len(x)-1):
        x_shear.append([x[i],x[i+1]-dx])
        y_shear.append([sum(y[:i+1])+x[i]*q, sum(y[:i+1])+(x[i+1]-dx)*q])
    x_shear.append(x[-1])
    y_shear.append(sum(y)+q*x[-1])

    f_get_val = lambda val: interp1d(x_shear[[val <= i[1] for i in x_shear[:-1]].index(True)],
                                y_shear[[val <= i[1] for i in x_shear[:-1]].index(True)], kind='linear')(val) if val != x_shear[-1] else y_shear[-1]

    import collections
    return lambda xx : [f_get_val(j) for j in xx] if isinstance(xx, collections.Iterable) else f_get_val(xx)


def int_moment(f_shear):
    # should be just integrating over shear(force) function
    import collections
    return lambda x: [integrate.quad(f_shear, 0, i)[0] for i in x] if isinstance(x, collections.Iterable) else integrate.quad(f_shear, 0, x)[0]


def int_moment_I(f_shear, max):
    import collections

    dt=0.00001
    x =  np.arange(0,max+dt,dt)
    lst_i  = [f_shear(0)*dt]
    for i in range(1, len(x)):
        lst_i.append(f_shear(x[i])*dt + lst_i[i-1])
    f = interp1d(x, lst_i, kind='linear')
    return lambda y: [f(i) for i in y] if isinstance(y, collections.Iterable) else f(y)

