import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

# Flow parameters
fv_w = 161.10420227050781  # Water Flow Rate (ml/min)
fv_a = 11.524953842163086  # Anhydride Flow Rate (ml/min)
v_cstr = 500  # Volume of CSTR (ml)

# Convert flows to units (ml/s)
fv_w_s = fv_w / 60
fv_a_s = fv_a / 60

# Make an array with the flows
flow_array = [fv_w_s, fv_a_s]

# Thermodynamic constants
params = {
    "flow": flow_array,
    "V": v_cstr / 1000,  # Volume in dm³
    "k0": 1.8725e3,          # Reaction rate constant (ml/mol/s)

    # Thermodynamic constants (taken from Asprey et al., 1996)
    "Ea": 45622.34,             # Activation energy (J/mol)
    "R": 8.314,              # Gas constant (J/mol/K)
    "H": 56.6e3,              # Enthalpy change (J/mol)
    "rho": 1,            # Density (g/ml)
    "cp": 4.186             # Heat capacity (J/g/K)
}

def der_func(t, C, **kwargs):
    # Initializing derivative vector
    dcdt = np.zeros(4)

    # Debug: print keys in kwargs
    print("Available keys in kwargs:", kwargs.keys())

    # Getting parameters
    flow = kwargs.get('flow', [0, 0])  # Use default value to avoid KeyError
    V = kwargs.get('V', 1)
    k0 = kwargs.get('k0', 0)
    Ea = kwargs.get('Ea', 0)
    R = kwargs.get('R', 1)
    H = kwargs.get('H', 0)
    rho = kwargs.get('rho', 1)
    cp = kwargs.get('cp', 1)

    # Reaction rate (assuming temperature in Kelvin)
    reaction_rate = C[0] * C[1] * k0 * np.exp(-Ea / (R * (C[3] + 273.15)))  # Convert C[3] to Kelvin
    total_flow = flow[0] + flow[1]

    # Differential equations
    dcdt[0] = (flow[0] / V) * (2 - C[0]) - reaction_rate  # Water Concentration derivative
    dcdt[1] = (flow[1] / V) * (1 - C[1]) - reaction_rate  # Anhydride Concentration derivative
    dcdt[2] = (total_flow / V) * (0 - C[2]) + 2 * reaction_rate  # Acetic acid concentration derivative
    dcdt[3] = 0  # Temperature change is assumed to be zero (isothermal)

    return dcdt

tspan = [0, 60]  # Time in seconds

# Initial conditions: concentrations of Water, Anhydride, Acetic Acid, Temperature
xini = [2, 1, 0, 25]  # Initial Conditions

# Solve ODE
sol = scipy.integrate.solve_ivp(der_func, tspan, xini, method='RK45', t_eval=np.linspace(tspan[0], tspan[1], 100), kwargs=params)

# Plotting results
plt.plot(sol.t, sol.y[0], label='Water (Ca)')
plt.plot(sol.t, sol.y[1], label='Anhydride (Cb)')
plt.plot(sol.t, sol.y[2], label='Acetic Acid (Cc)')
plt.xlabel('Time (min)')
plt.ylabel('Concentration (mol/L)')
plt.legend()
plt.title('Concentration of Various Components in a CSTR')
plt.show()
