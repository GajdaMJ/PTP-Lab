import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate


############### concentration at diferent locations ##############
# Assume reaction is 1st order wrt both components
# Assume isothermal (no exotherm)
# Assume constant density

#Reactor
r = 1.25 #(cm)


# Flow parameters
fv_w = 27  # Water Flow Rate (ml/min)
fv_a = 5.2  # Anhydride Flow Rate (ml/min)

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
    dcdv = np.zeros(4)

    # Getting parameters
    flow = parameters['flow']
    k0 = parameters['k0']
    Ea = parameters['Ea']
    R = parameters['R']
    H = parameters['H']
    rho = parameters['rho']
    cp = parameters['cp']
    
    reaction_rate = C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3])) *1e-3


    total_flow = flow[0]+flow[1]
    #Differential equations
    dcdv[0] =  - reaction_rate / total_flow # reaction_rate # Water Concentration derv
    dcdv[1] =  - reaction_rate / total_flow  # Anhydride Concentration derv
    dcdv[2] = 2*reaction_rate / total_flow # 2*reaction_rate # Acetic acid 
    dcdv[3] =  - H/(rho*cp) * reaction_rate # Temperature part
    return dcdv

# (total_flow/V) * (298-C[3]

tspan = [0,221] # Time in seconds

C_in_w = params['C_in_water']
C_in_AAH = params['C_in_AAH']

#For xini, c_water_0, c_AAH_0, C_AA_0, T0(in k)
xini = [C_in_w,C_in_AAH,0,T0+273.15] # Initial Conditions 

sol = scipy.integrate.solve_ivp(der_func, tspan, xini, args=(params,))

x = sol.t/(np.pi*r**2)
plt.plot(x, sol.y[0], label='Conc. Water')
plt.xlabel('Length')
plt.ylabel('Concentration(mol/mL)')
plt.legend()
plt.title('Concentration of water in a Steady State PFR')
plt.show()

plt.plot(x, sol.y[1], label='Conc. AAH')
plt.plot(x, sol.y[2], label='Conc. AA')
plt.xlabel('Length')
plt.ylabel('Concentration(mol/mL)')
plt.legend()
plt.title('Concentration of AAH, and AA in a Steady State PFR')
plt.show()

plt.plot(x, sol.y[3]-273.15, label = 'Temperature')
plt.xlabel('Length')
plt.ylabel('Temperature (K)')
plt.show()



# check the conversion and graph data due to suspicion of low aah conc

#mol_aah = tspan[1]*caah_pure*flow_array[1]


# x_data = list(np.arange(5, 57.5 + 7.5, 7.5))