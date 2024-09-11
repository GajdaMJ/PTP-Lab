import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

# Constants
fv_w = 161.10420227050781  # Water Flow Rate (ml/min)
fv_a = 11.524953842163086  # Anhydride Flow Rate (ml/min)
t_cstr = 25  # CSTR Temperature (°C)
k_1 = 10e5  # Reaction rate constant (dm^3/mol/s)
v_cstr = 500  # Volume of CSTR (ml)

# Convert flow rates to consistent units (ml/min to dm^3/s)
fv_w_dm3_s = fv_w / 1000 / 60  # Water flow rate in dm^3/s
fv_a_dm3_s = fv_a / 1000 / 60  # Anhydride flow rate in dm^3/s

# Initial concentrations (mol/L)
initial_conc_water = 1.0  # Example concentration of Water in mol/L
initial_conc_anhydride = 0.5  # Example concentration of Anhydride in mol/L
initial_conc_acid = 0.0  # Initial concentration of Acetic Acid (product) in mol/L
initial_temp = t_cstr  # Initial temperature in Celsius

# Arrhenius equation parameters
Ea = 50000  # Activation energy in J/mol
R = 8.314  # Universal gas constant in J/mol/K

# Define the differential function
def der_func(t, C, flow=(fv_w_dm3_s, fv_a_dm3_s), V=0.5, k0=k_1, Ea=Ea, R=R, H=-50000, rho=1000, cp=4.18):
    dcdt = np.zeros(4)
    T = C[3] + 273.15  # Convert temperature to Kelvin

    # Reaction rate using Arrhenius equation
    reaction_rate = C[0] * C[1] * k0 * np.exp(-Ea / (R * T))
    total_flow = flow[0] + flow[1]

    # Differential equations
    dcdt[0] = (flow[0] / V) * (initial_conc_water - C[0]) - reaction_rate  # Water concentration derivative
    dcdt[1] = (flow[1] / V) * (initial_conc_anhydride - C[1]) - reaction_rate  # Anhydride concentration derivative
    dcdt[2] = (total_flow / V) * (initial_conc_acid - C[2]) + 2 * reaction_rate  # Acetic acid concentration derivative
    dcdt[3] = (total_flow / V) * (t_cstr - C[3]) - (H / (rho * cp)) * reaction_rate  # Temperature change

    return dcdt

# Time span for integration (seconds)
tspan = [0, 60]  # 1 hour in seconds

# Initial conditions
xini = [initial_conc_water, initial_conc_anhydride, initial_conc_acid, initial_temp]

# Solve ODE
sol = scipy.integrate.solve_ivp(der_func, tspan, xini, t_eval=np.linspace(tspan[0], tspan[1], 100))

# Plotting results
plt.plot(sol.t, sol.y[0], label='Water (Ca)')
plt.plot(sol.t, sol.y[1], label='Anhydride (Cb)')
plt.plot(sol.t, sol.y[2], label='Acetic Acid (Cc)')
plt.xlabel('Time (s)')
plt.ylabel('Concentration (mol/L)')
plt.legend()
plt.title('Concentration of Various Components in a CSTR')
plt.show()
