import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from numpy.polynomial import Polynomial

############################################# MODELING THE ACTUAL DATA ######################################################
# Load data for different temperatures
my_data = np.genfromtxt('Data/PFR/25.09.30C.csv', delimiter=';', dtype=None, names=True, encoding='ISO-8859-1')  # Using ISO-8859-1 encoding

def temp_extract(data, x, offset=0):
    rows = data[data['TagName'] == x]
    # Remove invalid rows
    valid_rows = [row for row in rows if row['vValue'] not in ['(null)', None]]
    
    # Parse DateTime for valid rows
    date_times = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_rows]
    vvalues = [float(row['vValue'])+offset for row in valid_rows]
    
    # Calculate elapsed time in minutes
    start_time = date_times[0]
    elapsed_time = [(dt - start_time).total_seconds() / 60 for dt in date_times]

    return elapsed_time, vvalues

# Extract temperature data
elapsed_time_208 =  temp_extract(my_data, 'T208_PV')[0]
temp_208 = temp_extract(my_data, 'T208_PV', offset=-29.1)[1]


#################################################### MODELED DATA ##########################################################


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

    total_flow = flow[0] + flow[1]      #change the flow rate to be the total flow rate going into the next reactor
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
for tank_num in range(1, 8):
    sol_tank_prev = sol_tanks[-1][1]  # Take the previous tank solution>
    sol_tank_current = np.zeros_like(sol_tank_prev)  # Prepare to store the current solution
    
    for i, t in enumerate(sol_tanks[-1][0]):
        c_old = sol_tank_prev[i, :]  # Get previous tank's output at current time step
        sol_tank = master_function(lambda t, C: der_func_continue(t, C, params, c_old), [t, t + (tspan[1] - tspan[0]) / 100], c_old, method='rk4', number_of_points=1)
        sol_tank_current[i, :] = sol_tank[1][-1, :]  # Store current step

    sol_tanks.append((sol_tanks[-1][0], sol_tank_current))  # Append current tank's solution


######################################################## PLOTTING ###############################################################
# Plot the data

# Plot for 208 data
plt.plot(elapsed_time_208, temp_208, '-', label='T208 Raw Data', color='#ff7f0e')  # Orange for Raw Data
plt.plot(elapsed_time_208, sol_tanks[6][1][:, 3]-sol_tanks[6][1][0,3], color = 'red',label=f'Temperature Tank {8}')
plt.title('T208 Data')
plt.xlabel('Elapsed Time (min)')
plt.ylabel('Temperature (°C)')
plt.grid(True)
plt.legend()

plt.show()