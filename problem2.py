from numpy import log, exp # natural logarithm
import numpy as np
import timeit

import matplotlib.pyplot as plt

from scipy.integrate import quad

from values import rho_over_De

from utility_functions import ppmv_to_molm3



### VALUES SPECIFIC FOR THIS EXERCISE ###
T1 = 850 + 273  # [K] temperature

def y_SO2_label(C_SO2):
    return "$y_{SO_2,0}=$ " + str(C_SO2) + "ppmv" 

def alpha_label(alpha):
    return r"$\alpha$ = " + str(alpha)

def dp_title(Rp):
    if len(Rp) == 1:
        return "$d_p$ = " + str(Rp[0]*2*10**3) + " mm"
    else:
        return "$d_p$ mixed" 

desulph_str = "Desulphurization degree"


def integrand(u, t_mean, tau):
    
    return u**3/t_mean * exp(-tau/t_mean *(1 - 3*u**2 + 2*u**3)) * 6*tau*u*(u-1)



def find_X_CaO(alpha, C_SO2_0, t_mean, Rp):

    # Default value for conversion in case numerical method does not converge
    X_CaO = 1/alpha

    step_size = 0.0001

    X_CaO_iter = 1/alpha - step_size
    error = 1
    
    while X_CaO_iter > 0 and error>0.0005:

        X_SO2 = alpha * X_CaO_iter

        C_SO2_out = (1-X_SO2) * C_SO2_0

        mean_C_SO2 = (C_SO2_0-C_SO2_out) / log(C_SO2_0/C_SO2_out)

        tau = rho_over_De * Rp**2 / (6*mean_C_SO2)

        X_CaO = 1 - quad(integrand, 1, 0, args=(t_mean,tau))[0]

        error = np.abs(X_CaO - X_CaO_iter)
    
        X_CaO_iter -= step_size

    if X_CaO_iter <= 0:
        return 1/alpha
    else:
        return X_CaO


    
def desulphurization_degree(alpha, X):
    return alpha*X


def D_vs_y_SO2_plot(Rp, ax):

    print("\n***Plotting desulpurization degree versus C_SO2_0***\n")

    start_time_alpha = timeit.default_timer()

    # x-axis : concentration of C_SO2_0
    # y-axis: desulphurization degree
    # keeping t_mean fixed
    # varying alpha

    alpha_list = [1, 2, 3, 4]
    t_mean_fixed = 4

    for alpha in alpha_list:

        conc_points = 20

        C_SO2_0 = np.linspace(200, 1400, conc_points)

        D = np.ones(conc_points)

        # Assuming same amount of the particles with different diameters 
        w = 1/len(Rp)

        for i in range(conc_points):
            X_CaO = 0
            for R in Rp:
                X_CaO += w * find_X_CaO(alpha, ppmv_to_molm3(C_SO2_0[i], T1), t_mean_fixed, R)
            D[i] = desulphurization_degree(alpha, X_CaO)
        ax.set_title(dp_title(Rp))
        ax.plot(C_SO2_0, D, label=alpha_label(alpha))


        print("Done with alpha = ", alpha, " after t = ", timeit.default_timer() - start_time_alpha, " seconds")

    ax.set(xlabel="$y_{SO_2, 0}$ [ppmv]", ylabel=desulph_str)
    ax.set_ylim([0.3, 1.05])
    ax.legend(loc="lower right")


def D_vs_t_mean(Rp, ax):

    print("\n***Plotting desulhpurization degree versus t_mean***\n")

    start_time_t = timeit.default_timer()

    alpha_fixed = 2
    conc_list = [200, 500, 1000, 1400]

    for C_SO2_0 in conc_list:

        t_points = 20

        t = np.linspace(0.001, 8, t_points)

        D = np.ones(t_points)

        # Assuming same amount of particles of each diameter
        w = 1/len(Rp)
        

        for i in range(t_points):
            X_CaO = 0
            for R in Rp:
                X_CaO += w*find_X_CaO(alpha_fixed, ppmv_to_molm3(C_SO2_0, T1), t[i], R)
            D[i] = desulphurization_degree(alpha_fixed, X_CaO)

        ax.set_title(dp_title(Rp))
        
        ax.plot(t, D, label=y_SO2_label(C_SO2_0))
        
        print("Done with y_SO2 = ", C_SO2_0, " after t = ", timeit.default_timer() - start_time_t)

    ax.set(xlabel=r"$\overline{t} $ [h]", ylabel=desulph_str)
    ax.set_ylim([0, 1.05])
    ax.legend(loc="lower right", fontsize="small")


Rp_list = [[3*10**(-3)/2], [2*10**(-3)/2], [1*10**(-3)/2], [3*10**(-3)/2, 2*10**(-3)/2, 1*10**(-3)/2]]


fig1, axes1 = plt.subplots(1, len(Rp_list))
fig2, axes2 = plt.subplots(1, len(Rp_list))

for i in range(len(Rp_list)):
    #D_vs_y_SO2_plot(Rp_list[i], axes1[i])
    D_vs_t_mean(Rp_list[i], axes2[i])


for ax in axes1.flat:
    ax.label_outer()

for ax in axes2.flat:
    ax.label_outer()

plt.show()


