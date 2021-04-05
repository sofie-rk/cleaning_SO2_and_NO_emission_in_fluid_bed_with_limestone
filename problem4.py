from sympy import symbols, Eq, solve

from values import *
y_NO, y_NH3 = symbols('y_NO,y_NH3')


y_NO_in = 400
beta = 2
T = 850+273


def k(A, E, T):
    R = 8.314   # [J/K mol] 
    return A * np.exp(-E/(R*T))

k_ox = k(A1, E1, T)
k_r = k(A2, E2, T)


# Defining equations
eq1 = y_NO_in + t_cyclone*(k_ox*y_NH3 - k_r*10**6/Vm *y_NH3*y_NO) - y_NO

eq2 = beta