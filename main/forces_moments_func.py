from scipy.interpolate import interp1d
from itertools import chain
import scipy.integrate as integrate


def int_shear(x, y, q=0):
    dx = 0.0001  # small offset, for step jump

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
    return lambda x: [integrate.quad(f_shear, 0, i)[0] for i in x]
