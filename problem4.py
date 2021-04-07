from values import A1, A2, E1, E2
from values import R
from values import t_cyclone

from utility_functions import molar_volume

from numpy import exp
import numpy as np
import matplotlib.pyplot as plt

xlabel_str1 = "Temperature [$^o$C]"
xlabel_str2 = r"$y_{NO,in}$ [ppmv]"

def outlet_concentrations(T, beta, y_NO_in):

    # Temperature dependent variables
    Vm = molar_volume(T)
    k_ox = A1 * exp(-E1/(R*T))
    k_r = A2 * exp(-E2/(R*T)) / (Vm * 10**6)

    a = t_cyclone*k_r + 2*k_r*k_ox*t_cyclone**2
    b = k_r*t_cyclone*y_NO_in - beta*y_NO_in*t_cyclone*k_r + t_cyclone*k_ox + 1
    c = -beta*y_NO_in

    if (b**2 - 4*a*c < 0):
        print("Complex value!")
        return 0, 0

    y_NH3_out = (-b + (b**2 - 4*a*c)**(1/2))/(2*a)

    y_NH3_out1 = (-b + (b**2 - 4*a*c)**(1/2))/(2*a)
    y_NH3_out2 = (-b - (b**2 - 4*a*c)**(1/2))/(2*a)

    if (y_NH3_out1 > 0 and y_NH3_out2 > 0):
        print("!!!")

    if (y_NH3_out < 0):
        y_NH3_out = (-b - (b**2 - 4*a*c)**(1/2))/(2*a)
    
    y_NO_out = (beta*y_NO_in - t_cyclone*k_ox*y_NH3_out - y_NH3_out)/(k_r*t_cyclone*y_NH3_out)

    return y_NO_out, y_NH3_out 

print(outlet_concentrations(850+273, 2, 200))
    



    
beta_list = [1, 2, 3, 4]


fig1, (NO_ax1, NH3_ax1) = plt.subplots(1, 2)

for beta in beta_list:

    y_NO_in = 200

    T_points = 51
    T_list = np.linspace(750, 1200, T_points)

    NO_conc_list = np.zeros(T_points)
    NH3_conc_list = np.zeros(T_points)

    for j in range(T_points):
        NO_conc_list[j], NH3_conc_list[j] = outlet_concentrations(T_list[j] + 273, beta, y_NO_in)
    
    NO_ax1.plot(T_list, NO_conc_list, label=r"$\beta$ = " + str(beta))
    NH3_ax1.plot(T_list, NH3_conc_list, label=r"$\beta$ = " + str(beta))


NO_ax1.set(xlabel=xlabel_str1, ylabel="NO concentration [ppmv]")
NH3_ax1.set(xlabel=xlabel_str1, ylabel="NH3 concentration [ppmv]")
NO_ax1.set_title("NO concentration")
NH3_ax1.set_title("NH3 concentration")
fig1.suptitle("Keeping $y_{NO, in}$ = 200 ppmv")
NO_ax1.legend()
NH3_ax1.legend()
plt.show()




fig2, (NO_ax2, NH3_ax2) = plt.subplots(1, 2)

temperatures = [750, 900, 1000, 1200] # [C]

for T in temperatures:

    beta = 2

    y_NO_in_points = 51
    y_NO_in_list = np.linspace(25, 400, y_NO_in_points)

    NO_conc_list = np.zeros(y_NO_in_points)
    NH3_conc_list = np.zeros(y_NO_in_points)

    for j in range(y_NO_in_points):
        NO_conc_list[j], NH3_conc_list[j] = outlet_concentrations(T+273, beta, y_NO_in_list[j])

    print("T: ", T, " NO: ", NO_conc_list[-1])
    NO_ax2.plot(y_NO_in_list, NO_conc_list, label="T = " + str(T) + " $^o$C")
    NH3_ax2.plot(y_NO_in_list, NH3_conc_list, label="T = " + str(T) + " $^o$C")

NO_ax2.set_title("NO concentration")
NH3_ax2.set_title("NH3 concentration")
NO_ax2.set(xlabel=xlabel_str2, ylabel="NO concentration [ppmv]")
NH3_ax2.set(xlabel=xlabel_str2, ylabel="NH3 concentration [ppmv]")
fig2.suptitle("Keeping beta = 2")
NO_ax2.legend()
NH3_ax2.legend()
plt.show()

















