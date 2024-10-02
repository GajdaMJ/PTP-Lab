import numpy as np
import matplotlib.pyplot as plt

# Flow parameters
fv_w = 23.4317 # Water Flow Rate (ml/min)
fv_a = 5.4812  # Anhydride Flow Rate (ml/min)

v_pfr = 337/7  # Volume of CSTR (ml)
# Convert flow rates to consistent units (ml/min to ml/s)
fv_w_dm3_s = fv_w / 60  # Water flow rate in ml/s
fv_a_dm3_s = fv_a / 60  # Anhydride flow rate in ml/s

# Inlet Temperature
T0 = 30
 # (Celsius)

# Chemical constants
mm_water = 18.01528  # (g/mol)
rho_water = 0.999842  # (g/ml)
cw_pure = rho_water / mm_water

mm_AAH = 102.089  # (g/mol)
rho_AAH = 1.082  # (g/ml)
caah_pure = rho_AAH / mm_AAH

flow_array = [fv_w_dm3_s, fv_a_dm3_s]

    #Thermodynamic constants
params = {
    "C_in_water": (flow_array[0]*cw_pure)/(flow_array[0]+flow_array[1]),
    "C_in_AAH": (flow_array[1]*caah_pure)/(flow_array[0]+flow_array[1]),
    "Inlet temperature": T0+273.15,
    "flow": flow_array,
    "V": v_pfr,  # Volume in ml
    "k0": 1e7,          # Reaction rate constant (ml/mol/s) other possible values
    # Thermodynamic constants (taken from Asprey et al., 1996)
    "Ea": 45622.34,             # Activation energy (J/mol)
    "R": 8.314,              # Gas constant (J/mol/K)
    "H": -56.6e4,              # Enthalpy change (J/mol)
    "rho": 1,            # Density (g/ml)
    "cp": 4.186             # Heat capacity (J/g/K)
}

def der_func(t, C, parameters):
    dcdt = np.zeros(4)
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

    reaction_rate = C[0] * C[1] * k0 * np.exp(-Ea / (R * C[3])) * 1e-3

    total_flow = flow[0] + flow[1]
    dcdt[0] = (flow[0] / V) * (C_in_w - C[0]) - reaction_rate  # Water Concentration
    dcdt[1] = (flow[1] / V) * (C_in_AAH - C[1]) - reaction_rate  # Anhydride Concentration
    dcdt[2] = (total_flow / V) * (0 - C[2]) + 2 * reaction_rate  # Acetic acid
    dcdt[3] = (total_flow / V) * (inlet_temp - C[3]) - H / (rho * cp) * reaction_rate  # Temperature
    return dcdt


def der_func_continue(t, C, parameters, c_old):
    dcdt = np.zeros(4)
    flow = parameters['flow']
    V = parameters['V']
    k0 = parameters['k0']
    Ea = parameters['Ea']
    R = parameters['R']
    H = parameters['H']
    rho = parameters['rho']
    cp = parameters['cp']

    reaction_rate = C[0] * C[1] * k0 * np.exp(-Ea / (R * C[3])) * 1e-3

    total_flow = flow[0] + flow[1]
    dcdt[0] = (flow[0] / V) * (c_old[0] - C[0]) - reaction_rate  # Water
    dcdt[1] = (flow[1] / V) * (c_old[1] - C[1]) - reaction_rate  # Anhydride
    dcdt[2] = (total_flow / V) * (c_old[2] - C[2]) + 2 * reaction_rate  # Acetic acid
    dcdt[3] = (total_flow / V) * (c_old[3] - C[3]) - H / (rho * cp) * reaction_rate  # Temperature
    return dcdt


def master_function(fun, tspan, y0, method='rk4', number_of_points=100):
    '''General function to solve system of differential equations.'''
    dt = (tspan[1] - tspan[0]) / number_of_points
    t = np.linspace(tspan[0], tspan[1], number_of_points + 1)
    y = np.zeros((number_of_points + 1, len(y0)))  # len(y0) for each derivative

    # Set initial conditions
    y[0, :] = y0

    if method == 'midpoint':
        for i in range(number_of_points):
            k1 = fun(t[i], y[i, :])
            k2 = fun(t[i] + dt * 0.5, y[i, :] + 0.5 * dt * k1)
            y[i + 1, :] = y[i, :] + dt * k2
    elif method == 'euler':
        for i in range(number_of_points):
            y[i + 1, :] = y[i, :] + dt * fun(t[i], y[i, :])
    elif method == 'rk2':
        for i in range(number_of_points):
            k1 = fun(t[i], y[i, :])
            k2 = fun(t[i] + dt, y[i] + dt * k1)
            y[i + 1, :] = y[i] + dt * 0.5 * (k1 + k2)
    elif method == 'rk4':
        for i in range(number_of_points):
            k1 = fun(t[i], y[i, :])
            k2 = fun(t[i] + dt * 0.5, y[i, :] + 0.5 * dt * k1)
            k3 = fun(t[i] + dt * 0.5, y[i, :] + 0.5 * dt * k2)
            k4 = fun(t[i] + dt, y[i, :] + dt * k3)
            y[i + 1, :] = y[i] + dt * ((1 / 6) * k1 + (1 / 3) * (k2 + k3) + (1 / 6) * k4)
    else:
        raise ValueError('Unknown method specified. Check documentation for supported methods.')
    
    return t, y


tspan = [0, 6000]  # Time in seconds
xini = [cw_pure, 0, 0, T0 + 273.15]  # Initial conditions

# Simulate the first tank
sol_tank1 = master_function(lambda t, C: der_func(t, C, params), tspan, xini, method='rk4', number_of_points=100)

# Initialize list of solutions for tanks 2 to 7
sol_tanks = [sol_tank1]
for tank_num in range(2, 8):
    sol_tank_prev = sol_tanks[-1][1]  # Take the previous tank solution
    sol_tank_current = np.zeros_like(sol_tank_prev)  # Prepare to store the current solution
    
    for i, t in enumerate(sol_tanks[-1][0]):
        c_old = sol_tank_prev[i, :]  # Get previous tank's output at current time step
        sol_tank = master_function(lambda t, C: der_func_continue(t, C, params, c_old), [t, t + (tspan[1] - tspan[0]) / 100], c_old, method='rk4', number_of_points=1)
        sol_tank_current[i, :] = sol_tank[1][-1, :]  # Store current step

    sol_tanks.append((sol_tanks[-1][0], sol_tank_current))  # Append current tank's solution



# # Print the first temperature value for the first reactor
# first_temp_reactor_i = sol_tanks[0][1][0, 3] - 273.15  # Subtract 273 to convert from Kelvin to Celsius
# print(f"The first temperature for the first reactor is: {first_temp_reactor_i:.2f} °C")

# # Print the 20th temperature value for the first reactor
# xth_temp_reactor_i = sol_tanks[0][1][2000, 3] - 273.15  # Subtract 273 to convert from Kelvin to Celsius
# print(f"The xth temperature for the first reactor is: {xth_temp_reactor_i:.2f} °C")

all_temp_reactor_i = sol_tanks[0][1][:, 3] - 273.15  # Subtract 273 to convert from Kelvin to Celsius
print(f"The xth temperature for the first reactor is: {all_temp_reactor_i} °C")



# Create subplots
fig, ax = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

# Plotting concentration of AAH for each tank
for i in range(7):
    ax[0].plot(sol_tanks[i][0], sol_tanks[i][1][:, 1], label=f'Conc. AAH Tank {i + 1}')

ax[0].set_ylabel('Concentration (mol/mL)')
ax[0].set_title('Concentration of AAH in 7 Reactors in Series')
ax[0].legend()  # Add legend to differentiate tanks

# Plotting temperature for each tank
for i in range(7):
    ax[1].plot(sol_tanks[i][0], sol_tanks[i][1][:, 3]-sol_tanks[i][1][0, 3], label=f'Temperature Tank {i + 1}')

ax[1].set_xlabel('Time (s)')
ax[1].set_ylabel('Temperature (K)')
ax[1].set_title('Temperature in 7 Reactors in Series')
ax[1].legend()  # Add legend to differentiate tanks

# Show the plots
plt.show()
