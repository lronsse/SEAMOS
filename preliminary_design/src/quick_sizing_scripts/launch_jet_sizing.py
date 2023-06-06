"""
A quick script for sizing the launch jet of the Buffin
Its messy, but thats fine, you shouldnt need it anyways
"""


##  params
from math import pi, sqrt
import matplotlib.pyplot as plt
rho_w = 1000. # kg/m3
p_a = 1e5 # Pa
p_r = 8e5 # Pa
V_r = 0.001 # m3
V_c = 0.001 # m3
A_t = (0.02/2)**2 * pi # m2
dt = 0.001
m_test = 16  #kg

plot_t = True
plot_imp = True

if plot_t:
    t = 0.
    V = 0.
    p = p_r

    parr = []
    Ftarr = []
    tarr = []
    while t < .25:
        p = p_r * V_r / (V_r+V)
        if p > p_a and V < V_c:
            Vdot = A_t * sqrt(2/rho_w*(p-p_a))
        else:
            Vdot = 0.
        V += dt * Vdot
        parr.append(p)
        if Vdot != 0.:
            Ft = (Vdot*rho_w)**2/A_t/rho_w + A_t*(p-p_a)
        else:
            Ft = 0.
        Ftarr.append(Ft)
        t += dt
        tarr.append(t)

    I_tot = sum(Ftarr) * dt
    print('total impulse:', I_tot, '[Ns]')
    m_test = 16.
    print(f'dV for {m_test}kg: {I_tot / m_test} m/s')


    plt.plot(tarr, parr)
    plt.show()
    plt.plot(tarr, Ftarr)
    plt.show()

if plot_imp:
    i = 0
    for vc in [0.0005, 0.001, 0.002, 0.005]:
        Iarr = []
        prarr = [5e5, 15e5, 30e5, 60e5, 120e5, 240e5]
        V_c = vc
        for p_r in prarr:
            t = 0.
            V = 0.
            p = p_r

            parr = []
            Ftarr = []
            tarr = []
            while t < .25:
                p = p_r * V_r / (V_r + V)
                if p > p_a and V < V_c:
                    Vdot = A_t * sqrt(2 / rho_w * (p - p_a))
                else:
                    Vdot = 0.
                V += dt * Vdot
                parr.append(p)
                if Vdot != 0.:
                    Ft = (Vdot * rho_w) ** 2 / A_t / rho_w + A_t * (p - p_a)
                else:
                    Ft = 0.
                Ftarr.append(Ft)
                t += dt
                tarr.append(t)

            I_tot = sum(Ftarr) * dt
            Iarr.append(I_tot)

        print(Iarr)
        varr = [I_i/m_test for I_i in Iarr]
        prarr = [5, 15, 30, 60, 120, 240]
        style = ['-', '--', '-.', ':']
        plt.plot(prarr, varr, label=f'{V_c*1000} L', linestyle=style[i])
        i += 1

    plt.legend()
    plt.xscale('log')
    plt.ylabel('Delivered dV [m/s]')
    plt.xlabel('Initial Pressurant Pressure [Bar]')
    plt.show()

