import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import scipy.integrate
import math
# Assume reaction is 1st order wrt both components
# Assume isothermal (no exotherm)
# Assume constant density

def PBR_model(T,fv1,fv2, V=131, tspan = [0,3600], n=6):
    '''Models the behavior of the reaction: Water + Acetic Anhydride -> 2 * Acetic acid in an adiabatic CSTR reactor. \n
    Required Arguments: \n
    T = inlet temperature for the reactor given in units celsius \n
    fv1 = flow rate of water in units ml/min \n
    fv2 = flow rate of acetic anhydride ml/min \n
    Optional Arguments: \n
    V = volume of the reactor in units ml (default set to 500ml) \n
    tspan = list of evaluation time in units seconds (default set to [0,3600]) \n
    This function was built for the course "Practical Process Technology (6P4X0)" 
    '''
    v_pfr_tank = V/n

    # Convert flow rates (ml/min to ml/s)
    fv_w_dm3_s = fv1 / 60  # Water flow rate in ml/s
    fv_a_dm3_s = fv2  / 60  # Anhydride flow rate in ml/s
    #dont use ml/s use l/s
    #Chemical constants

    #Water
    mm_water = 18.01528 # (g/mol)
    rho_water = 0.999842 # (g/ml)
    cw_pure = rho_water/mm_water # (mol/ml)

    #Acetic acid
    mm_AAH = 102.089 # (g/mol)
    rho_AAH = 1.082 # (g/ml)
    caah_pure = rho_AAH/mm_AAH # (mol/ml)

    flow_array = [fv_w_dm3_s, fv_a_dm3_s]

    # Calculations for glass beads
    V_total = 337 #cm3
    V_beads = 337-V #cm3 should be like 206cm3
    epsilon = (V_total-V_beads)/(V_total) #void fraction
    diameter_bead = 2e-1 # 2mm diameter bu i want it in cm
    A_total = (3*V_beads*diameter_bead)/2
    A_per_tank = A_total/n


    
    params = { # Stores the relevant thermodynamic constants as a dictionary 
        "C_in_water": (flow_array[0]*cw_pure)/(flow_array[0]+flow_array[1]),
        "C_in_AAH": (flow_array[1]*caah_pure)/(flow_array[0]+flow_array[1]),
        "Inlet temperature": T+273.15, # Temp but now in kelvin
        "flow": flow_array,
        "V": v_pfr_tank,  # Volume in ml
        "k0": 6.5e5,#np.exp(16),#7e6,          # Reaction rate constant (ml/mol/s)

        # Thermodynamic constants (taken from Asprey et al., 1996)
        "Ea": 45187.2,             # Activation energy (J/mol)
        "R": 8.314,              # Gas constant (J/mol/K)
        "H": -56.6e3,              # Enthalpy change (J/mol)
        "rho_water": 1,            # Density (g/ml)
        "rho_glass": 2.4,          # Density (g/ml)
        "epsilon" : epsilon, 
        "cp_water": 4.186,             # Heat capacity (J/g/K)
        "cp_glass": 0.84,            #Heat capacity
        "Area_bead_per_tank": A_per_tank, # Area of beads per "tank"
        "U" : 1.2122e-4#0.12122 W/m2*K but we want in W/cm2*K so e-4 
    }

    xini_temp = [cw_pure,0,0,T+273.15, T+273.15] # Initial Conditions 
    xini = np.zeros(5*n)
    for i in range(5*n):
        if np.mod(i,5)==0:
            xini[i] = xini_temp[0]
        elif np.mod(i,5)==1:
            xini[i] = xini_temp[1]
        elif np.mod(i,5)==2:
            xini[i] = xini_temp[2]
        elif np.mod(i,5)==3:
            xini[i] = xini_temp[3]
        elif np.mod(i, 5)==4:
            xini[i] = xini_temp[4]

    sol_me = scipy.integrate.solve_ivp(der_func, tspan, xini, t_eval=np.linspace(tspan[0], tspan[1], 400), args=(params, n)) 
    return sol_me

def der_func(t,C, parameters, n=6):
    '''This function contains the differential equations to solve the reaction A+B->2C in an adiabatic 
    CSTR. \n
    t=time (seconds) \n
    c = Concentration vector like [c_water, c_AAH, c_AA, Temperature]\n
    parameters = dictionary containing thermodynamic constants
    '''
    # Initializing derivative vector
    dcdt = np.zeros(5*n)
    # array of 4 zeros corresponding to [c_water/dt, c_AAH/dt, c_AA/dt, dT/dt]

    # Getting parameters out of our dictionary 
    C_in_w = parameters['C_in_water']
    C_in_AAH = parameters['C_in_AAH']
    flow = parameters['flow']
    V = parameters['V']
    k0 = parameters['k0']
    Ea = parameters['Ea']
    R = parameters['R']
    H = parameters['H']
    rho_water = parameters['rho_water']
    rho_glass = parameters['rho_glass']
    epsilon = parameters['epsilon']
    cp_water = parameters['cp_water']
    cp_glass = parameters['cp_glass']
    inlet_temp = parameters["Inlet temperature"]
    A = parameters["Area_bead_per_tank"]
    U = parameters["U"] 

    total_flow = flow[0]+flow[1]
    
    #Differential equations
    for i in range(5*n):
        if i < 5:  # For the first set of concentrations and temperatures
           dcdt[0] = (total_flow / V) * (C_in_w - C[0]) - C[0] * C[1] * k0 * np.exp(-Ea / (R * C[3]))  # Water Concentration derivative
           dcdt[1] = (total_flow / V) * (C_in_AAH - C[1]) - C[0] * C[1] * k0 * np.exp(-Ea / (R * C[3]))  # Anhydride Concentration derivative
           dcdt[2] = (total_flow / V) * (0 - C[2]) + 2 * C[0] * C[1] * k0 * np.exp(-Ea / (R * C[3]))  # Acetic acid concentration derivative
           dcdt[3] = (total_flow / V) * (inlet_temp - C[3]) - H / (rho_water * cp_water) * C[0] * C[1] * k0 * np.exp(-Ea / (R * C[3])) + (U * A) / (rho_water * cp_water * V) * (C[i+1] - C[3])  # Reactor temperature derivative
           # Glass bead temperature derivative
           dcdt[4] = (U * A) / (rho_glass * cp_glass * V) * (C[3] - C[4])  # Temperature change of glass beads
        else:
            # Loop for additional reactors
            if np.mod(i, 5) == 0:
                dcdt[i] = (total_flow / V) * (C[i - 5] - C[i]) - C[i] * C[i + 1] * k0 * np.exp(-Ea / (R * C[i + 3])) #Water
            elif np.mod(i, 5) == 1:
                dcdt[i] = (total_flow / V) * (C[i - 5] - C[i]) - C[i - 1] * C[i] * k0 * np.exp(-Ea / (R * C[i + 2])) #AAH
            elif np.mod(i, 5) == 2:
                dcdt[i] = (total_flow / V) * (C[i - 5] - C[i]) + 2 * C[i - 2] * C[i - 1] * k0 * np.exp(-Ea / (R * C[i + 1])) # AA
            elif np.mod(i, 5) == 3:
                dcdt[i] = (total_flow / V) * (C[i - 5] - C[i]) - H / (rho_water * cp_water) * C[i - 3] * C[i - 2] * k0 * np.exp(-Ea / (R * C[i])) + (U * A) / (rho_water * cp_water * V) * (C[i+1] - C[i]) 
            elif np.mod(i, 5) == 4:
                dcdt[i] = (U * A) / (rho_glass * cp_glass * V) * (C[i - 1] - C[i])  # Temperature change of glass beads for additional reactors
    return dcdt

def temp_extract(data, x="T200_PV", offset=0):
    '''Function to extract data from csv files\n
    data = data path for your csv file. Give as a string \n
    x = Name of the instrument that you want. Default set to T200_PV (CSTR internal temperature) \n
    offset = linear offset for values. Default set to zero \n
    returns elapsed time and values for your
    '''
    # Extract the flow data to determine the starting time
    flow_rows = data[data['TagName'] == "P120_Flow"]
    valid_flow_rows = [row for row in flow_rows if row['vValue'] not in ['(null)', None]]
    flow_values = [float(row['vValue']) for row in valid_flow_rows]
    flow_dates = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_flow_rows]
    start_time = None
    for i in range(1, len(flow_values)):
        if flow_values[i-1] < 1 and flow_values[i] > 1:     # Loop that checks when the AAH pump is turned on and sets that as the start time
            start_time = flow_dates[i]
            break # Stop the loop once flow starts

    temp_rows = data[data['TagName'] == x]  # Only choose the rows for that particular instrument 
    valid_temp_rows = [row for row in temp_rows if row['vValue'] not in ['(null)', None]] # You want to remove the values when theres null otherwise it does weird things
    
    temp_dates = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_temp_rows] #Converts the weird csv time format to python
    temp_values = [float(row['vValue']) + offset for row in valid_temp_rows]

    # Calculate elapsed time in minutes from the start_time
    elapsed_time = [(dt - start_time).total_seconds() / 60 for dt in temp_dates]

    return elapsed_time, temp_values

def data_extract(data, x, offset=0):
    # Extract the flow data to determine the starting time
    flow_rows = data[data['TagName'] == "P120_Flow"]
    valid_flow_rows = [row for row in flow_rows if row['vValue'] not in ['(null)', None]]
    flow_values = [float(row['vValue']) for row in valid_flow_rows]
    flow_dates = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_flow_rows]

    start_time = None
    for i in range(1, len(flow_values)):
        if flow_values[i-1] < 1 and flow_values[i] > 1:
            start_time = flow_dates[i]
            break
    
    # Extract temperature data starting from the transition point
    temp_rows = data[data['TagName'] == x]
    valid_temp_rows = [row for row in temp_rows if row['vValue'] not in ['(null)', None]]
    
    temp_dates = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_temp_rows]
    temp_values = [float(row['vValue']) + offset for row in valid_temp_rows]

    elapsed_time = [(dt - start_time).total_seconds() / 60 for dt in temp_dates] # Calculate elapsed time in minutes from the start_time

    return elapsed_time, temp_values, (flow_dates[0] - start_time).total_seconds() / 60 


if __name__ == '__main__':
    my_data = np.genfromtxt('Data\PFR\\18.09.40C_again.csv', delimiter=';', dtype=None, names=True, encoding='ISO-8859-1')

    # Extracting all temperature data
    t_values = ['T208_PV','T207_PV','T206_PV','T205_PV','T204_PV','T203_PV','T202_PV','T201_PV','T200_PV']
    results = {}

    for t_value in t_values:
        elap_time, temp_c, offset_time = data_extract(my_data, t_value)
        results[t_value] = {'elapsed_time': elap_time, 'temperature': temp_c, 'offset_time':offset_time}

    # Get AAH Flowrate and Water Flowrate
    elapsed_time_c_aah, aah_flowrate_c_vector, offset_time = data_extract(my_data, x="P120_Flow")
    elapsed_time_c_water, water_flowrate_c_vector, offset_time = data_extract(my_data, x='P100_Flow')

    # Find initial temperature and flowrates
    initial_temperature = np.min(temp_c)
    aah_flowrate_c = np.median(aah_flowrate_c_vector)
    water_flowrate_c = np.median(water_flowrate_c_vector)
    
    n_tanks=16

    # Run PBR model simulation
    sol_me = PBR_model(initial_temperature, water_flowrate_c, aah_flowrate_c, V=131, tspan=[0, 3600], n=n_tanks)

    # Create subplots for each reactor stage
    fig, ax = plt.subplots(2, 4, figsize=(20, 8), sharex=True, sharey=True)
    ax = ax.flatten()

    retention_time = 2 + 2 / 60  # minutes

    for i in range(0, 8):
        # Extract experimental temperature data
        temp_data = np.array(results[t_values[-(i + 1)]]['temperature'])
        elapsed_time = results[t_values[-(i + 1)]]['elapsed_time']
        if i == 0:
            tank = 1
        else:
            tank = math.ceil((i * n_tanks) / (8))+1

        # Plot real temperature data
        ax[i].plot(elapsed_time, temp_data - temp_data[0], color='#ff7f0e', label='Real Data', linewidth=2)

        # Plot model temperature data for the corresponding stage
        ax[i].plot(sol_me.t / 60, sol_me.y[3 + tank * 5, :] - 273.15 - initial_temperature, color='#1f77b4', label='Model Prediction', linewidth=2)

        # Set plot title, labels, and grid
        ax[i].set_title(f'Temperature Probe {i + 1}, Reactor {tank + 1}', fontsize=14, fontweight='bold')
        ax[i].set_xlabel('Elapsed Time (min)', fontsize=12)
        ax[i].set_ylabel('Change in Temperature (°C)', fontsize=12)
        ax[i].set_xlim(0, 20)
        ax[i].grid(True)
        ax[i].legend(fontsize=10)

    # Set global title and adjust layout
    fig.suptitle('Reactor Temperature Data Comparison', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()