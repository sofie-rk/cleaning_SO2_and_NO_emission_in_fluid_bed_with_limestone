import numpy as np 
import matplotlib.pyplot as plt

from scipy.integrate import odeint


from values import *


from utility_functions import *


# FIXED VALUES
X_CaO = 0.2
X_N = 1.5/100


def k(A, E, Ea, T):
    R = 8.314   # [J/K mol] 
    return A * np.exp(-(E - Ea)/(R*T))

k_ox = k(A1, E1, E1a, T)
k_r = k(A2, E2, E2a, T)



def A_surf(X_CaO):
    return S_0 * (1-X_CaO) * C_CaO * (1-e)/e


def O2_conc(t):

    if t>t_bed:
        return 0.02 * p/(R*T)

    else:
        return (0.21 - 0.19/t_bed*t)*p/(R*T)


def model(X, t):
    NH3 = X[0]
    NO = X[1]

    A = A_surf(X_CaO)

    O2 = O2_conc(t)

    if t<=t_bed:
        dNH3_dt = 0.6/t_bed*X_N - A*k_14*NH3 - A*k2*NH3*O2 - A*2/3*k3*NH3*NO - k_ox*NH3 - k_r*NH3*NO
        dNO_dt = A*k2*NH3*O2 - A*k3*NH3*NO + k_ox*NH3 - k_r*NH3*NO
    else:
        dNH3_dt = -k_ox*NH3 - k_r*NH3*NO
        dNO_dt = k_ox*NH3 - k_r*NH3*NO


    return [dNH3_dt, dNO_dt]




# INITIAL CONDITIONS
conc_initial = [0, 0] 

# Time points
n = 5000

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

    X = odeint(model, conc_initial, tspan)

    NH3_conc[i] = X[1][0]
    NO_conc[i] = X[1][1]

    # New init conditions
    conc_initial = X[1]


plt.plot(t, conc_to_ppmv(NH3_conc), label="NH3")
plt.plot(t, conc_to_ppmv(NO_conc), label="NO")
plt.xlabel("Time t [s]")
plt.ylabel("Concentration [ppmv]")
#plt.plot([t_bed, t_bed], [0, 300], '--', color = "black")
plt.legend()
plt.show()

print(conc_to_ppmv(NO_conc[-1]))
# x-axis: X_CaO
# will vary X_N




