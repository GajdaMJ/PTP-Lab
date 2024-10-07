import numpy as np
import matplotlib.pyplot as plt

x_end = 1 # meter
tend = 4
d = 1e-3
kr = 1.2 #per second

# mew = volumetric flow / cross sectional area
diam = 2.5 # meter
flow = 6.3e-4 #m^3/s
mew = flow / (np.pi*(diam**2)/4)


def model_2_pbr(T, fv1,fv2, C, V=137,n):
    '''Models the behavior of the reaction: Water + Acetic Anhydride -> 2 * Acetic acid in an adiabatic CSTR reactor. \n
    Required Arguments: \n
    T = inlet temperature for the reactor given in unit celsius \n
    C = values [Water, AAH, Acitic Acid, Temperature] \n
    V = volume of the PBR reactor \n
    n = number of tanks in series \n
    '''
    #Chemical constants
    mm_water = 18.01528 # (g/mol)
    rho_water = 0.999842 # (g/ml)
    cw_pure = rho_water/mm_water # (mol/ml)
    mm_AAH = 102.089 # (g/mol)
    rho_AAH = 1.082 # (g/ml)
    caah_pure = rho_AAH/mm_AAH # (mol/ml)

   # Convert flow rates (ml/min to ml/s)
    fv_w_dm3_s = fv1 / 60  # Water flow rate in ml/s
    fv_a_dm3_s = fv2  / 60  # Anhydride flow rate in ml/s

    flow_array = [fv_w_dm3_s, fv_a_dm3_s]
    
    v_pfr = V/n
    #Thermodynamic constants
    params = {
        "C_in_water": (flow_array[0]*cw_pure)/(flow_array[0]+flow_array[1]),
        "C_in_AAH": (flow_array[1]*caah_pure)/(flow_array[0]+flow_array[1]),
        "Inlet temperature": T+273.15,
        "flow": flow_array,
        "V": v_pfr,  # Volume in ml
        "k0": 6e5,          # Reaction rate constant (ml/mol/s)

        # Thermodynamic constants (taken from Asprey et al., 1996)
        "Ea": 45622.34,             # Activation energy (J/mol)
        "R": 8.314,              # Gas constant (J/mol/K)
        "H": -56.6e3,              # Enthalpy change (J/mol)
        "rho": 1,            # Density (g/ml)
        "cp": 4.186             # Heat capacity (J/g/K)
    }

    #initialize the solution matrix
    dcdt = np.zeros(4*n)

    # Getting parameters
    C_in_w = params['C_in_water']
    C_in_AAH = params['C_in_AAH']
    flow = params['flow']
    V = params['V']
    k0 = params['k0']
    Ea = params['Ea']
    R = params['R']
    H = params['H']
    rho = params['rho']
    cp = params['cp']
    inlet_temp = params["Inlet temperature"]

    total_flow = flow[0]+flow[1]
    #solving for each reactor
    for i in range(4*n):
        #for the first reactor it takes the initial conditions
        if i ==0:
            dcdt[i] = (total_flow/V)*(C_in_w - C[0]) - C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3]))
        if i == 1:
            dcdt[i] = (total_flow/V)*(C_in_AAH - C[1])  - C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3]))
        if i == 2:
            dcdt[i] = (total_flow/V)*(0 - C[2]) + 2*C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3])) 
        if i == 3:
            dcdt[i] = (total_flow/V) * (inlet_temp-C[3]) - H/(rho*cp) * C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3]))
        #make a condition for the remain reactors which will not be taking the initial conditions but the outlet of the reactor before
        
        elif np.mod(): #np.mod(divident,divisor) with the output being the remainder
    
        