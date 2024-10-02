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
elapsed_time_208, temp_208 = temp_extract(my_data, 'T208_PV')
elapsed_time_207, temp_207 = temp_extract(my_data, 'T207_PV')
elapsed_time_206, temp_206 = temp_extract(my_data, 'T206_PV')
elapsed_time_205, temp_205 = temp_extract(my_data, 'T205_PV')
elapsed_time_204, temp_204 = temp_extract(my_data, 'T204_PV')
elapsed_time_203, temp_203 = temp_extract(my_data, 'T203_PV')
elapsed_time_202, temp_202 = temp_extract(my_data, 'T202_PV')
elapsed_time_201, temp_201 = temp_extract(my_data, 'T201_PV')


# Generate evenly spaced time values for interpolation
elapsed_time_interp = np.linspace(min(elapsed_time_208), max(elapsed_time_208), 1000)

# Polynomial fit for the data (degree 5 for a better fit, especially on plateaus)
poly_208 = Polynomial.fit(elapsed_time_208, temp_208, 20)
poly_207 = Polynomial.fit(elapsed_time_207, temp_207, 20)
poly_206 = Polynomial.fit(elapsed_time_206, temp_206, 20)
poly_205 = Polynomial.fit(elapsed_time_205, temp_205, 20)
poly_204 = Polynomial.fit(elapsed_time_204, temp_204, 20)
poly_203 = Polynomial.fit(elapsed_time_203, temp_203, 20)
poly_202 = Polynomial.fit(elapsed_time_202, temp_202, 20)
poly_201 = Polynomial.fit(elapsed_time_201, temp_201, 20)

# Use the polynomial to generate smooth curves
temp_208_smooth = poly_208(elapsed_time_interp)
temp_207_smooth = poly_207(elapsed_time_interp)
temp_206_smooth = poly_206(elapsed_time_interp)
temp_205_smooth = poly_205(elapsed_time_interp)
temp_204_smooth = poly_204(elapsed_time_interp)
temp_203_smooth = poly_203(elapsed_time_interp)
temp_202_smooth = poly_202(elapsed_time_interp)
temp_201_smooth = poly_201(elapsed_time_interp)

#convirting the temps into arrays, so that they the delta T can be plotted
temp_208 = np.array(temp_208)
temp_207 = np.array(temp_207)
temp_206 = np.array(temp_206)
temp_205 = np.array(temp_205)
temp_204 = np.array(temp_204)
temp_203 = np.array(temp_203)
temp_202 = np.array(temp_202)
temp_201 = np.array(temp_201)




#################################################### MODELED DATA ##########################################################
# Flow parameters
fv_w = 23.4317 # Water Flow Rate (ml/min)
fv_a = 5.4812  # Anhydride Flow Rate (ml/min)

v_pfr = 131/7  # Volume of CSTR (ml)
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
    "k0": 1e6, #np.exp(16.25)          # Reaction rate constant (ml/mol/s)
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


######################################################## PLOTTING ###############################################################

# Plot the data
fig, ax = plt.subplots(2, 4, figsize=(30, 10), sharex=True, sharey=True)

# Plot for 208 data
ax[0,0].plot(elapsed_time_208, temp_208 - temp_208[0], 'o', label='T208 Raw Data', color='#1f77b4')  # Blue
ax[0,0].plot(elapsed_time_interp, temp_208_smooth - temp_208[0], '-', color='#ff7f0e', label='T208 Smoothed Curve')  # Orange
ax[0,0].set_title('T208 Data')
ax[0,0].set_xlabel('Elapsed Time (min)')
ax[0,0].set_ylabel('Temperature (°C)')
ax[0,0].grid(True)
ax[0,0].legend()

# Plot for 207 data
ax[0,1].plot(elapsed_time_207, temp_207 - temp_207[0], 'o', label='T207 Raw Data', color='#2ca02c')  # Green
ax[0,1].plot(elapsed_time_interp, temp_207_smooth - temp_207[0], '-', color='#d62728', label='T207 Smoothed Curve')  # Red
ax[0,1].set_title('T207 Data')
ax[0,1].set_xlabel('Elapsed Time (min)')
ax[0,1].grid(True)
ax[0,1].legend()

# Plot for 206 data
ax[0,2].plot(elapsed_time_206, temp_206 - temp_206[0], 'o', label='T206 Raw Data', color='#9467bd')  # Purple
ax[0,2].plot(elapsed_time_interp, temp_206_smooth -temp_206[0], '-', color='#8c564b', label='T206 Smoothed Curve')  # Brown
ax[0,2].set_title('T206 Data')
ax[0,2].set_xlabel('Elapsed Time (min)')
ax[0,2].grid(True)
ax[0,2].legend()

# Plot for 205 data
ax[0,3].plot(elapsed_time_205, temp_205 - temp_205[0], 'o', label='T205 Raw Data', color='#e377c2')  # Pink
ax[0,3].plot(elapsed_time_interp, temp_205_smooth - temp_205[0], '-', color='#7f7f7f', label='T205 Smoothed Curve')  # Gray
ax[0,3].set_title('T205 Data')
ax[0,3].set_xlabel('Elapsed Time (min)')
ax[0,3].grid(True)
ax[0,3].legend()

# Plot for 204 data
ax[1,0].plot(elapsed_time_204, temp_204- temp_204[0], 'o', label='T204 Raw Data', color='#bcbd22')  # Yellow-green
ax[1,0].plot(elapsed_time_interp, temp_204_smooth-temp_204[0], '-', color='#17becf', label='T204 Smoothed Curve')  # Cyan
ax[1,0].set_title('T204 Data')
ax[1,0].set_xlabel('Elapsed Time (min)')
ax[1,0].set_ylabel('Temperature (°C)')
ax[1,0].grid(True)
ax[1,0].legend()

# Plot for 203 data
ax[1,1].plot(elapsed_time_203, temp_203 - temp_203[0], 'o', label='T203 Raw Data', color='#ff9896')  # Light Red
ax[1,1].plot(elapsed_time_interp, temp_203_smooth - temp_203[0], '-', color='#c49c94', label='T203 Smoothed Curve')  # Light Brown
ax[1,1].set_title('T203 Data')
ax[1,1].set_xlabel('Elapsed Time (min)')
ax[1,1].grid(True)
ax[1,1].legend()

# Plot for 202 data
ax[1,2].plot(elapsed_time_202, temp_202 - temp_202[0], 'o', label='T202 Raw Data', color='#aec7e8')  # Light Blue
ax[1,2].plot(elapsed_time_interp, temp_202_smooth - temp_202[0], '-', color='#ffbb78', label='T202 Smoothed Curve')  # Light Orange
ax[1,2].set_title('T202 Data')
ax[1,2].set_xlabel('Elapsed Time (min)')
ax[1,2].grid(True)
ax[1,2].legend()

# Plot for 201 data
ax[1,3].plot(elapsed_time_201, temp_201 - temp_201[0], 'o', label='T201 Raw Data', color='#98df8a')  # Light Green
ax[1,3].plot(elapsed_time_interp, temp_201_smooth - temp_201[0], '-', color='#f7b6d2', label='T201 Smoothed Curve')  # Light Pink
ax[1,3].set_title('T201 Data')
ax[1,3].set_xlabel('Elapsed Time (min)')
ax[1,3].grid(True)
ax[1,3].legend()

# Adjust layout to prevent label overlap and set a global title
fig.suptitle('T200_PV Temperature Data over Time', fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust for title space

# Show the plot
plt.show()
