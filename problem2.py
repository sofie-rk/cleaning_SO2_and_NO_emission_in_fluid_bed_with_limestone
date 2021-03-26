from numpy import log, exp # natural logarithm
import numpy as np

import matplotlib.pyplot as plt

from scipy.integrate import quad

from values import rho_over_De

from utility_functions import ppmv_to_molm3


def integrand(u, t_mean, tau):
    
    return u**3/t_mean * exp(-tau/t_mean *(1 - 3*u**2 + 2*u**3)) * 6*tau*u*(u-1)



def find_X_CaO(alpha, C_SO2_0, t_mean, Rp):

    # Default value for conversion in case numerical method does not converge
    X_CaO = 1/alpha

    step_size = 0.0005

    X_CaO_iter = 1/alpha - step_size
    error = 1
    
    while X_CaO_iter > 0 and error>0.01:

        X_SO2 = alpha * X_CaO_iter

        C_SO2_out = (1-X_SO2) * C_SO2_0

        mean_C_SO2 = (C_SO2_0-C_SO2_out) / log(C_SO2_0/C_SO2_out)

        tau = rho_over_De * Rp**2 / (6*mean_C_SO2)

        X_CaO = 1 - quad(integrand, 1, 0, args=(t_mean,tau))[0]

        error = np.abs(X_CaO - X_CaO_iter)
    
        X_CaO_iter -= step_size

    return X_CaO
    
def desulphurization_degree(alpha, X):
    return alpha*X

### WITH d_p = 3mm
def plot_desulpu_degree_vs_C_SO2_0(Rp):

    # x-axis : concentration of C_SO2_0
    # y-axis: desulphurization degree
    # keeping t_mean fixed
    # varying alpha

    print("\n***Plotting desulpurization degree versus C_SO2_0***\n")

    t_mean_fixed = 4
    number_of_points = 200

    # x-axis is C_SO2_0
    C_SO2_0 = np.linspace(200, 1400, number_of_points)

    # Empty desulphurization arrays
    D1 = np.zeros(number_of_points)
    D2 = np.zeros(number_of_points)
    D3 = np.zeros(number_of_points)
    D4 = np.zeros(number_of_points)

    for i in range(number_of_points):

        # alpha = 1
        X1 = find_X_CaO(1, ppmv_to_molm3(C_SO2_0[i]), t_mean_fixed, Rp)
        D1[i] = desulphurization_degree(1, X1)

        # alpha = 2
        X2 = find_X_CaO(2, ppmv_to_molm3(C_SO2_0[i]), t_mean_fixed, Rp)
        D2[i] = desulphurization_degree(2, X2)

        # alpha = 3
        X3 = find_X_CaO(3, ppmv_to_molm3(C_SO2_0[i]), t_mean_fixed, Rp)
        D3[i] = desulphurization_degree(3, X3)

        # alpha = 4
        X4 = find_X_CaO(4, ppmv_to_molm3(C_SO2_0[i]), t_mean_fixed, Rp)
        D4[i] = desulphurization_degree(4, X4)

    plt.plot(C_SO2_0, D1, label=r"$\alpha = 1$") 
    plt.plot(C_SO2_0, D2, label=r"$\alpha = 2$") 
    #plt.plot(C_SO2_0, D3, label="alpha = 3")
    #plt.plot(C_SO2_0, D4, label="alpha = 4")
    plt.xlabel("Inlet concentration of SO2 $C(SO_2,0) [ppmv]$")
    plt.ylabel("Desulpurization degree")
    plt.legend()
    plt.show()

plot_desulpu_degree_vs_C_SO2_0(3*10**(-3)/2)


