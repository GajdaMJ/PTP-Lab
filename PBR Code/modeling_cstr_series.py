#modeling the pfr as a cstr system in series with 7 cstr
import numpy as np
import matplotlib.pyplot as plt 
import scipy.integrate


# Flow parameters
fv_w = 161.10420227050781 # Water Flow Rate (ml/min)
fv_a = 11.524953842163086 # Anhydride Flow Rate (ml/min)
v_cstr = 1 # Volume of CSTR (ml)

#Convert flows to units (ml/s)
fv_w_s = fv_w/60
fv_a_s = fv_a/60

#Make an array with the flows
flow_array = [fv_w_s,fv_a_s]

#Thermodynamic constants
params = {
    "flow": flow_array,
    "V": v_cstr,  # Volume in ml
    "k0": 0.02, #1.8725e3,          # Reaction rate constant (ml/mol/s)

    # Thermodynamic constants (taken from Asprey et al., 1996)
    "Ea": 45622.34,             # Activation energy (J/mol)
    "R": 8.314,              # Gas constant (J/mol/K)
    "H": 56.6e3,              # Enthalpy change (J/mol)
    "rho": 1,            # Density (g/ml)
    "cp": 4.186             # Heat capacity (J/g/K)
}



def der_func(t,C, parameters):
    # Initializing derivative vector
    dcdt = np.zeros(4)

    # Getting parameters
    flow = parameters['flow']
    V = parameters['V']
    k0 = parameters['k0']
    Ea = parameters['Ea']
    R = parameters['R']
    H = parameters['H']
    rho = parameters['rho']
    cp = parameters['cp']

    #reaction_rate = C[0]*C[1] * k0  * np.exp(-Ea/R*C[3])
    total_flow = flow[0]+flow[1]

    #Differential equations
    dcdt[0] =  (flow[0]/V)*(2 - C[0])    - k0*C[0]*C[1] # reaction_rate # Water Concentration derv
    dcdt[1] =  (flow[1]/V)*(1 - C[1])    - k0*C[0]*C[1] # reaction_rate # Anhydride Concentration derv
    dcdt[2] =  (total_flow/V)*(0 - C[2]) + 2*k0*C[0]*C[1] # 2*reaction_rate # Acetic acid 
    dcdt[3] =  0#(total_flow/V) * (ini_cond[3]-C[3]) - H/(rho*cp) * reaction_rate # Temperature part
    return dcdt


sol_y = []
sol_t = []

tspan = [0, 60]  # Time in seconds

for n in range(7):
    if n < 1:
        xini = [3, 0, 0, 0]  # Initial Conditions
    else:
        # Use the last values from the previous solution as initial conditions
        xini = [sol_y[-1][i][-1] for i in range(4)]  # Access last time-step of previous solution

    sol = scipy.integrate.solve_ivp(der_func, tspan, xini, args=(params,))
    
    # Append solutions to the lists without reassigning the append result
    sol_y.append(sol.y)  # sol.y is the array of results
    sol_t.append(sol.t)  # sol.t is the time points

print(sol_y)


