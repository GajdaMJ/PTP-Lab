import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

# Flow parameters
fv_w = 161.10420227050781  # Water Flow Rate (ml/min)
fv_a = 11.524953842163086  # Anhydride Flow Rate (ml/min)
v_cstr = 500  # Volume of CSTR (ml)

#Inlet Temperature 
T0 = 25 # (celsius)

# Convert flow rates to consistent units (ml/min to ml/s)
fv_w_dm3_s = fv_w  / 60  # Water flow rate in ml/s
fv_a_dm3_s = fv_a  / 60  # Anhydride flow rate in ml/s

#Chemical constants
mm_water = 18.01528 # (g/mol)
rho_water = 0.999842 # (g/ml)
cw_pure = rho_water/mm_water

mm_AAH = 102.089 # (g/mol)
rho_AAH = 1.082 # (g/ml)
caah_pure = rho_AAH/mm_AAH

flow_array = [fv_w_dm3_s, fv_a_dm3_s]

cw0 = (flow_array[0]*cw_pure)/(flow_array[0]+flow_array[1])
caah0 = (flow_array[0]*caah_pure)/(flow_array[0]+flow_array[1])

# Thermodynamic constants
t_cstr = 25  # CSTR Temperature (°C)
k_1 = 0.02  # Reaction rate constant (dm³/mol/s)
# k_1 is a function of both temperature and

# Define the differential function
def der_func(t, C, flow=flow_array, V=1, k0=k_1):
    dcdt = np.zeros(4)
    
    # Reaction rate (assuming isothermal conditions, no temperature dependency)
    reaction_rate = k0 * C[0] * C[1]
    total_flow = flow[0] + flow[1]

    # Differential equations
    dcdt[0] = (flow[0] / V) * (2 - C[0]) - reaction_rate  # Water concentration derivative
    dcdt[1] = (flow[1] / V) * (1 - C[1]) - reaction_rate  # Anhydride concentration derivative
    dcdt[2] = (total_flow / V) * (0 - C[2]) + 2 * reaction_rate  # Acetic acid concentration derivative
    dcdt[3] = 0  # No temperature change assumed (isothermal)

    return dcdt

# Time span for integration (minutes)
tspan = [0, 60]  # Time in minutes


# Initial conditions: concentrations of Water, Anhydride, Acetic Acid, Temperature
xini = [2.0, 1.0, 0.0, 25]  # Initial concentrations (mol/L) and initial temperature (°C)


# Solve ODE
sol = scipy.integrate.solve_ivp(der_func, tspan, xini)

# Plotting results
plt.plot(sol.t, sol.y[0], label='Water (Ca)')
plt.plot(sol.t, sol.y[1], label='Anhydride (Cb)')
plt.plot(sol.t, sol.y[2], label='Acetic Acid (Cc)')
plt.xlabel('Time (min)')
plt.ylabel('Concentration (mol/mL)')
plt.legend()
plt.title('Concentration of Various Components in a CSTR')
plt.show()
