from values import Vm

# Finding concentration in [mol/m^3] when concentration is
# given in ppmv, using molar volume

def ppmv_to_molm3(conc_ppmv):
    return conc_ppmv*10**(-6)/Vm

def conc_to_ppmv(conc_molm3):
    return conc_molm3*10**6*Vm