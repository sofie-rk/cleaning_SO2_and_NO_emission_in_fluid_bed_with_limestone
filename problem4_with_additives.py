from values import A1, A2, E1, E2
from values import R
from values import t_cyclone

from utility_functions import molar_volume

from numpy import exp
import numpy as np
import matplotlib.pyplot as plt

xlabel_str = "$y_{C_2H_6}$ [ppmv]"

T_4_2 = 800 +273    # [K] Temperature
y_NO_in = 500       # [ppmv] inlet concentration of NO

def outlet_concentrations(beta, y_C2H6):

    E1a = (-0.59938*y_C2H6**2 + 268.6491*y_C2H6)
    E2a = (-0.43523*y_C2H6**2 + 179.98454*y_C2H6)

    # Temperature dependent variables
    Vm = molar_volume(T_4_2)
    k_ox = A1 * exp(-(E1-E1a)/(R*T_4_2))
    k_r = A2 * exp(-(E2-E2a)/(R*T_4_2)) / (Vm * 10**6)

    a = t_cyclone*k_r + 2*k_r*k_ox*t_cyclone**2
    b = k_r*t_cyclone*y_NO_in - beta*y_NO_in*t_cyclone*k_r + t_cyclone*k_ox + 1
    c = -beta*y_NO_in

    if (b**2 - 4*a*c < 0):
        print("Complex value!")
        return 0, 0

    y_NH3_out = (-b + (b**2 - 4*a*c)**(1/2))/(2*a)

    if (y_NH3_out < 0):
        y_NH3_out = (-b - (b**2 - 4*a*c)**(1/2))/(2*a)
    
    y_NO_out = (beta*y_NO_in - t_cyclone*k_ox*y_NH3_out - y_NH3_out)/(k_r*t_cyclone*y_NH3_out)

    return y_NO_out, y_NH3_out 



beta_list = [1, 2, 3, 4]


fig1, (NO_ax1, NH3_ax1) = plt.subplots(1, 2)

for beta in beta_list:

    y_C2H6_points = 20

    y_C2H6_list = np.linspace(0, 200, y_C2H6_points)

    NO_conc_list = np.zeros(y_C2H6_points)
    NH3_conc_list = np.zeros(y_C2H6_points)

    for j in range(y_C2H6_points):
        NO_conc_list[j], NH3_conc_list[j] = outlet_concentrations(beta, y_C2H6_list[j])

    NO_ax1.plot(y_C2H6_list, NO_conc_list, label=r"$\beta$ = " + str(beta))
    NH3_ax1.plot(y_C2H6_list, NH3_conc_list, label=r"$\beta$ = " + str(beta))

NO_ax1.set(xlabel=xlabel_str, ylabel="NO concentration [ppmv]")
NH3_ax1.set(xlabel=xlabel_str, ylabel="NH3 concentration [ppmv]")
NO_ax1.set_title("NO concentration")
NH3_ax1.set_title("NH3 concentration")
NO_ax1.legend(loc="upper right")
NH3_ax1.legend()
plt.show()



