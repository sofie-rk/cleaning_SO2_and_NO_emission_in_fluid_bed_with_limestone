

def molar_volume(T):
    R = 8.314   # [J/Kmol]
    p = 1.01325*10**5 # [Pa] total pressure
    return (R*T)/p

def ppmv_to_molm3(conc_ppmv, T):       # from [ppmv]  
    Vm = molar_volume(T)
    return conc_ppmv*10**(-6)/Vm    # to [mol/m3]

def conc_to_ppmv(conc_molm3, T):   # from [mol/m3]
    Vm = molar_volume(T)
    return conc_molm3*10**(6)*Vm  # to [ppmv]