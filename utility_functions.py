from values import Vm


def ppmv_to_molm3(conc_ppmv):       # from [ppmv]  
    return conc_ppmv*10**(-6)*Vm    # to [mol/m3]

def conc_to_ppmv(conc_molm3):   # from [mol/m3]
    return conc_molm3*10**6/Vm  # to [ppmv]