import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

# Assume reaction is 1st order wrt both components
# Assume isothermal (no exotherm)
# Assume constant density

# Flow parameters
fv_w = 161.10420227050781  # Water Flow Rate (ml/min)
fv_a = 11.524953842163086  # Anhydride Flow Rate (ml/min)

v_cstr = 200  # Volume of CSTR (ml)
# Convert flow rates to consistent units (ml/min to ml/s)
fv_w_dm3_s = fv_w  / 60  # Water flow rate in ml/s
fv_a_dm3_s = fv_a  / 60  # Anhydride flow rate in ml/s

#Inlet Temperature 
T0 = 25 # (celsius)

#Chemical constants
mm_water = 18.01528 # (g/mol)
rho_water = 0.999842 # (g/ml)
cw_pure = rho_water/mm_water

mm_AAH = 102.089 # (g/mol)
rho_AAH = 1.082 # (g/ml)
caah_pure = rho_AAH/mm_AAH

flow_array = [fv_w_dm3_s, fv_a_dm3_s]


#Thermodynamic constants
params = {
    "C_in_water": (flow_array[0]*cw_pure)/(flow_array[0]+flow_array[1]),
    "C_in_AAH": (flow_array[1]*caah_pure)/(flow_array[0]+flow_array[1]),
    "Inlet temperature": T0+273.15,
    "flow": flow_array,
    "V": v_cstr,  # Volume in ml
    "k0": np.exp(16.25)*1000,          # Reaction rate constant (ml/mol/s)

    # Thermodynamic constants (taken from Asprey et al., 1996)
    "Ea": 45622.34,             # Activation energy (J/mol)
    "R": 8.314,              # Gas constant (J/mol/K)
    "H": -56.6e3,              # Enthalpy change (J/mol)
    "rho": 1,            # Density (g/ml)
    "cp": 4.186             # Heat capacity (J/g/K)
}


def der_func(t,C, parameters):
    # Initializing derivative vector
    dcdt = np.zeros(4)

    # Getting parameters
    C_in_w = parameters['C_in_water']
    C_in_AAH = parameters['C_in_AAH']
    flow = parameters['flow']
    V = parameters['V']
    k0 = parameters['k0']
    Ea = parameters['Ea']
    R = parameters['R']
    H = parameters['H']
    rho = parameters['rho']
    cp = parameters['cp']
    inlet_temp = parameters["Inlet temperature"]
    
    reaction_rate = C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3])) 


    total_flow = flow[0]+flow[1]
    #Differential equations
    dcdt[0] =  (flow[0]/V)*(C_in_w) - (total_flow/V)*C[0]    - reaction_rate # reaction_rate # Water Concentration derv
    dcdt[1] =  (flow[1]/V)*(C_in_AAH) - (total_flow/V)*C[1]  - reaction_rate  # Anhydride Concentration derv
    dcdt[2] =  (total_flow/V)*(0 - C[2]) + 2*reaction_rate # 2*reaction_rate # Acetic acid 
    dcdt[3] =  (total_flow/V) * (inlet_temp-C[3]) - H/(rho*cp) * reaction_rate # Temperature part
    return dcdt



tspan = [0,600] # Time in seconds

#For xini, c_water_0, c_AAH_0, C_AA_0, T0(in k)
xini = [cw_pure,0,0,T0+273.15] # Initial Conditions 

sol = scipy.integrate.solve_ivp(der_func, tspan, xini, args=(params,))

#plt.plot(sol.t, sol.y[0], label='Conc. Water')
plt.plot(sol.t, sol.y[1], label='Conc. AAH')
plt.plot(sol.t, sol.y[2], label='Conc. AA')
plt.xlabel('time')
plt.ylabel('Concentration(mol/mL)')
plt.legend()
plt.title('Concentration of various components in a CSTR')
plt.show()

plt.plot(sol.t, sol.y[3])
plt.show()



# check the conversion and graph data due to suspicion of low aah conc

mol_aah = tspan[1]*caah_pure*flow_array[1]
