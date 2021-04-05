import numpy as np 
import matplotlib.pyplot as plt

from scipy.integrate import odeint

# from values import T, p, R
# from values import A1, A2, E1, E2, E1a, E2a

from values import *


from utility_functions import *



def k(A, E, Ea, T):
    R = 8.314   # [J/K mol] 
    return A * np.exp(-(E - Ea)/(R*T))

k_ox = k(A1, E1, E1a, T)
k_r = k(A2, E2, E2a, T)


def A_surf(X_CaO):
    return S_0 * (1-X_CaO) *C_CaO * (1-e)/e


def O2_conc(t):

    return (0.21 - 0.19/t_bed*t)*p/(R*T)


def dNH3_dt(t, conc_NH3, conc_NO):
    return 0

def dNO_dt(t, conc_NH3, conc_NO):
    return 0


def model(X, t, X_N, X_CaO):
    NH3 = X[0]
    NO = X[1]

    A = A_surf(X_CaO)

    O2 = O2_conc(t)

    if t<=t_bed:
        dNH3_dt = 0.6/t_bed*X_N - A*k_14*NH3 - A*k2*NH3*O2 - A*2/3*k3*NH3*NO - k_ox*NH3 - k_r*NH3*NO
        dNO_dt = A*k2*NH3*O2 - A*k3*NH3*NO + k_ox*NH3 - k_r*NH3*NO
    else:
        dNH3_dt = 0.6/t_bed*X_N - k_ox*NH3 - k_r*NH3*NO
        dNO_dt = k_ox*NH3 - k_r*NH3*NO

    return [dNH3_dt, dNO_dt]


def NH3_out(X_N, X_CaO):


    # INITIAL CONDITIONS
    conc_initial = [0, 0] 

    # Time points
    n = 3000

    t = np.linspace(0, t_free+t_bed , n)

    # Store sol
    NH3_conc = np.empty_like(t)
    NO_conc = np.empty_like(t)

    # Add initial condition
    NH3_conc[0] = conc_initial[0]
    NO_conc[0] = conc_initial[1]

    ### SOLVING ###
    for i in range(1, n):
        tspan = [t[i-1], t[i]] 

        X = odeint(model, conc_initial, tspan, args=(X_N, X_CaO))

        NH3_conc[i] = X[1][0]
        NO_conc[i] = X[1][1]

        # New init conditions
        conc_initial = X[1]

    return NO_conc[-1]



X_N = [1/100, 1.5/100, 2/100]

def X_N_label(X_N):
    return "$X_N$ = " + str(X_N)

for i in range(len(X_N)):

    number_X_CaO_points = 21
    X_CaO_list = np.linspace(0, 1, number_X_CaO_points)
    NO_conc_plot = np.zeros(number_X_CaO_points)

    for j in range(number_X_CaO_points):
        NO_conc_plot[j] = conc_to_ppmv(NH3_out(X_N[i], X_CaO_list[j]))
    
    print(NO_conc_plot)

    plt.plot(X_CaO_list, NO_conc_plot, label=X_N_label(X_N[i]))


plt.legend()
plt.xlabel("$X_{CaO}$ [-]")
plt.ylabel("NO emmision [ppmv]")
plt.show()









