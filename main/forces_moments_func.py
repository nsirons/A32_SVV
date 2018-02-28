from scipy.interpolate import interp1d
from itertools import chain
import scipy.integrate as integrate
import numpy as np


def int_shear(x, y, q=0):
    dx = 0.00000000001  # small offset, for step jump

    if x[0] != 0:
        x = [0] + x
        y = [0] + y

    x_shear = list(chain(*[[x[i], x[i + 1] - dx] for i in range(len(x) - 1)])) + [x[-1]]
    y_shear = [
        sum(y[:x.index(x_new) + 1]) + q * x[x.index(x_new)] if x_new in x else sum(y[:x.index(x_new + dx)]) + q * x_new
        for x_new in x_shear]

    return interp1d(x_shear, y_shear, kind='linear')


def int_moment(f_shear):
    # should be just integrating over shear(force) function
    import collections
    return lambda x: [integrate.quad(f_shear, 0, i)[0] for i in x] if isinstance(x, collections.Iterable) else integrate.quad(f_shear, 0, x)[0]

# example
# s_f = int_shear([pos],[forces], distributed load])
# s_f = int_shear([0, 10, 20], [10, 20, -30], 0)
# m_f = int_moment(s_f)
# print(s_f(3))
# print(s_f([1,2,3]))
# print(m_f(3))
# print(m_f([1,2,3]))