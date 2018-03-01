import numpy as np
from tools import heaviside_defl


def get_deflections_func(x1,x2,x3,xa, E, Izz,Iyy,Izy,Ay, Az, By,Bz, Cy, Dz, K1,K2,p,q):

    c = -1/(E*(Izz*Iyy-Izy**2))
    v_f = lambda x: \
        -c * ((heaviside_defl(x-x1)*(x-x1))**3)/6.*Iyy * Ay \
        -c * ((heaviside_defl(x-x1)*(x-x1))**3)/6.*Izy * Az \
        -c * ((heaviside_defl(x-x2)*(x - x2)) ** 3) / 6. * Iyy * By \
        -c * ((heaviside_defl(x - x2)*(x-x2)) ** 3) / 6. * Izy * Bz \
        -c * ((heaviside_defl(x - x3) * (x - x3)) ** 3) / 6. * Iyy * Cy \
        - c * ((heaviside_defl(x - (x2 - xa / 2.))*(x - (x2 - xa / 2.))) ** 3) / 6. * Izy * Dz \
        +x*c*K1 + c*K2 \
        + c*q*(x**4)/24.*Iyy + c*p/6*((heaviside_defl(x-(x2+xa/2.))*(x-(x2+xa/2.)))**3)*Izy

    A = np.array([[c*x1, c*1], [c*x2, c*1]])
    b = np.array([[-c*q*x1**4/24],
                  [-c*q*x2**4/24 + c*(x2-x1)**3 / 6 * (-Izy)*Ay /
                                + c*(x2-x1)**3 / 6 * (-Izz)*Az /
                                + c*(x2 - (x2 - xa / 2.)) ** 3 / 6 * (-Izz) * Dz]])

    D1, D2 = np.linalg.solve(A,b)

    u_f = lambda x: \
        -c * ((heaviside_defl(x-x1)*(x-x1))**3)/6.*(-Izy) * Ay \
        -c * ((heaviside_defl(x-x1)*(x-x1))**3)/6.*(-Izz) * Az \
        -c * ((heaviside_defl(x-x2)*(x - x2)) ** 3) / 6. * (-Izy) * By \
        -c * ((heaviside_defl(x - x2)*(x-x2)) ** 3) / 6. * (-Izy) * Bz \
        -c * ((heaviside_defl(x - x3) * (x - x3)) ** 3) / 6. * (-Izy) * Cy \
        - c * ((heaviside_defl(x - (x2 - xa / 2.))*(x - (x2 - xa / 2.))) ** 3) / 6. * (-Izz) * Dz \
        +x*c*D1 + c*D2 \
        + c*q*(x**4)/24.*(-Izy) + c*p/6*((heaviside_defl(x-(x2+xa/2.))*(x-(x2+xa/2.)))**3)*(-Izz)

    return v_f, u_f

