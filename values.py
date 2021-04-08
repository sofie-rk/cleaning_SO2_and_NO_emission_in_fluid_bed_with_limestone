from utility_functions import ppmv_to_molm3

# Reactor conditions
T = 850 + 273       # [K]
p = 1.01325*10**5   # [Pa] total pressure
t_bed = 0.5         # [s] residence time in the bed
t_free = 1.0        # [s] residence time in the free board
t_cyclone = 1.5     # [s] residence time in the cyclone


# Particle conditions
S_0 = 2.15 * 10**7  # [m^2 CaO / m^3 CaO]
C_CaO = 0.15        # [m^3 CaO / m^3 solid]
e = 0.63            # [m^3 gas / m^3 bed] bed porosity


# Constants
R = 8.314 # [J/Kmol] gas constant


# Rate constants, catalytic reactions
k_14 = 4.9 * 10**(-6)   # [m^3gas/m^2CaO s]
k2 = 1.6*10**(-5)       # [m^6gas/m^2CaO mol s]
k3 = 6.8*10**(-5)       # [m^6gas/m^2CaO mol s]

# Rate constants, gas phase reactions
A1 = 2.21*10**(14)  # [1/s]
E1 = 317.3*10**3    # [J/mol] activation energy
A2 = 2.45*10**(14)  # [m^3/mol s]
E2 = 244.4 *10**3   # [J/mol] activation energy
E1a = 9.47*10**3    # [J/mol] activation energy correction
E2a = 5.37 * 10**3  # [J/mol] activation energy correction


# From problem 1
d1 = 2*10**(-3)     # [m] particle diameter
r1 = d1/2           # [m] particle radius
C_SO2_1 = ppmv_to_molm3(1200, 850+273) # [mol/m3] concentration of SO2
tau1 = 225*60       # [s] residence time

rho_over_De = tau1*6*C_SO2_1/(r1**2*60*60) # [per hour]

