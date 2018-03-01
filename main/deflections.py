import numpy as np
from tools import heaviside_defl


def get_deflections_func(x1,x2,x3,xa, E, Izz,Iyy,Izy,Ay, Az, By,Bz, Cy, Dz, K1,K2,p,q):

    c = -1/(E*(Izz*Iyy-(Izy**2)))
    Mz = lambda x:   ( ((heaviside_defl(x-x1)*(x-x1)) ** 3) / 6. * Ay \
                       +((heaviside_defl(x-x2)*(x-x2)) ** 3) / 6. * By \
                       +((heaviside_defl(x-x3)*(x-x3)) ** 3) / 6. * Cy \
                       -q*(x**4)/24.)
    My = lambda x:   ( ((heaviside_defl(x-x1)*(x-x1)) ** 3) / 6. * Az \
                       +((heaviside_defl(x-x2)*(x-x2)) ** 3) / 6. * Bz \
                       +((heaviside_defl(x-(x2-xa/2.))*(x-(x2-xa/2.))) ** 3) / 6. * Dz \
                       -((heaviside_defl(x-(x2+xa/2.))*(x-(x2+xa/2.))) ** 3)*p/6.)

    v_f = lambda x: c*(-Mz(x)*Iyy - My(x)*Izy + K1*x + K2)


    A = np.array([[c*x1, c*1], [c*x2, c*1]])
    b = np.array([[c*(My(x1)*Izz + Mz(x1)*Izy)],
                  [c*(My(x2)*Izz + Mz(x2)*Izy)]])

    D1, D2 = np.linalg.solve(A,-b)

    u_f = lambda x: c*(My(x)*Izz + Mz(x)*Izy + D1*x + D2)

    return v_f, u_f

