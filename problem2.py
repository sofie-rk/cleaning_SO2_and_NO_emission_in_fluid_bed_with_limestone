from numpy import log, exp # natural logarithm
import numpy as np
import timeit

import matplotlib.pyplot as plt

from scipy.integrate import quad

from values import rho_over_De

from utility_functions import ppmv_to_molm3



### VALUES SPECIFIC FOR THIS EXERCISE ###
T1 = 850 + 273  # [K] temperature

def C_SO2_label(C_SO2):
    return "$C(SO_2)_0 =$ " + str(C_SO2) + " ppmv" 

def alpha_label(alpha):
    return r"$\alpha$ = " + str(alpha)

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

def generating_plots(Rp):

    fig, (plot1ax, plot2ax) = plt.subplots(1, 2)

    print("\n***Plotting desulpurization degree versus C_SO2_0***\n")

    start_time_alpha = timeit.default_timer()

    # x-axis : concentration of C_SO2_0
    # y-axis: desulphurization degree
    # keeping t_mean fixed
    # varying alpha

    alpha_list = [1, 2, 3, 4]
    t_mean_fixed = 4

    for alpha in alpha_list:

        conc_points = 50

        C_SO2_0 = np.linspace(200, 1400, conc_points)

        D = np.zeros(conc_points)

        for i in range(conc_points):
            X_CaO = find_X_CaO(alpha, ppmv_to_molm3(C_SO2_0[i], T1), t_mean_fixed, Rp)
            D[i] = desulphurization_degree(alpha, X_CaO)

        plot1ax.plot(C_SO2_0, D, label=alpha_label(alpha))

        print("Done with alpha = ", alpha, " after t = ", timeit.default_timer() - start_time_alpha, " seconds")



    print("\n***Plotting desulpurization degree versus t_mean***\n")

    start_time_t = timeit.default_timer()

    alpha_fixed = 2
    conc_list = [200, 500, 1000, 1400]

    for C_SO2_0 in conc_list:

        t_points = 50

        t = np.linspace(0.001, 8, t_points)

        D = np.zeros(t_points)

        for i in range(t_points):
            X_CaO = find_X_CaO(alpha_fixed, ppmv_to_molm3(C_SO2_0, T1), t[i], Rp)
            D[i] = desulphurization_degree(alpha_fixed, X_CaO)
        
        plot2ax.plot(t, D, label=C_SO2_label(C_SO2_0))
        
        print("Done with C_SO2 = ", C_SO2_0, " after t = ", timeit.default_timer() - start_time_t)




    #plot1ax.plot([200, 1400], [1, 1], '--', color="black")
    

    plot1ax.set(xlabel="Inlet concentration of SO2 $C(SO_2)_0$ ppmv", ylabel=desulph_str)
    plot1ax.set_title("Keeping " + r"$\overline{t}_{mean}$ = " + str(t_mean_fixed))

    plot2ax.set(xlabel="Mean residence time " + r"$\overline{t}$ [h]", ylabel=desulph_str)
    plot2ax.set_title("Keeping " + r"$\alpha = $" + str(alpha_fixed))

    plot1ax.legend()
    plot2ax.legend()

    fig.suptitle("With diameter = " + str(Rp*2*10**3) + " [mm]")

    plt.show()



generating_plots(3*10**(-3)/2)


