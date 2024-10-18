import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import scipy.integrate
# Assume reaction is 1st order wrt both components
# Assume isothermal (no exotherm)
# Assume constant density

def CSTR_model(T1, T2,fv1,fv2, V=500, tspan = [0,3600], t_change=1800):
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
    v_cstr=V

    # Convert flow rates (ml/min to ml/s)
    fv_w_dm3_s = fv1 / 60  # Water flow rate in ml/s
    fv_a_dm3_s = fv2  / 60  # Anhydride flow rate in ml/s

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

    
    params = { # Stores the relevant thermodynamic constants as a dictionary 
        "C_in_water": (flow_array[0]*cw_pure)/(flow_array[0]+flow_array[1]),
        "C_in_AAH": (flow_array[1]*caah_pure)/(flow_array[0]+flow_array[1]),
        "Inlet temperature": T1+273.15, # Temp but now in kelvin
        "flow": flow_array,
        "V": v_cstr,  # Volume in ml
        "k0": 7e6,          # Reaction rate constant (ml/mol/s)

        # Thermodynamic constants (taken from Asprey et al., 1996)
        "Ea": 45622.34,             # Activation energy (J/mol)
        "R": 8.314,              # Gas constant (J/mol/K)
        "H": -56.6e3,              # Enthalpy change (J/mol)
        "rho": 1,            # Density (g/ml)
        "cp": 4.186             # Heat capacity (J/g/K)
    }
    # print(params['C_in_AAH']*params['C_in_water'])
    xini = [cw_pure,0,0,T1+273.15] # Initial Conditions 
    sol_1 =  scipy.integrate.solve_ivp(der_func, [tspan[0], t_change], xini, args=(params,)) 
    
    #Change over
    params["Inlet temperature"] = T2+273.15 #Change temp
    xini = [sol_1.y[0,-1], sol_1.y[1,-1], sol_1.y[2,-1], sol_1.y[3,-1]]   #Take last row as initial conditions for next solution iteration 
    sol_2 = scipy.integrate.solve_ivp(der_func, [t_change,tspan[1]], xini, args=(params,))

    combined_time = np.concatenate((sol_1.t, sol_2.t))  # Combine time points
    combined_y = np.concatenate((sol_1.y, sol_2.y), axis=1)  # Combine solution arrays along axis 1 (columns)

    return combined_time, combined_y

def der_func(t,C, parameters):
    '''This function contains the differential equations to solve the reaction A+B->2C in an adiabatic 
    CSTR. \n
    t=time (seconds) \n
    c = Concentration vector like [c_water, c_AAH, c_AA, Temperature]\n
    parameters = dictionary containing thermodynamic constants
    '''
    # Initializing derivative vector
    dcdt = np.zeros(4)
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
    rho = parameters['rho']
    cp = parameters['cp']
    inlet_temp = parameters["Inlet temperature"]
    
    reaction_rate = C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3])) # reaction rate is repeated so just calculate once

    total_flow = flow[0]+flow[1]
    
    #Differential equations
    dcdt[0] =  (total_flow/V)*(C_in_w - C[0])    - reaction_rate # Water Concentration derv
    dcdt[1] =  (total_flow/V)*(C_in_AAH - C[1])  - reaction_rate  # Anhydride Concentration derv
    dcdt[2] =  (total_flow/V)*(0 - C[2]) + 2*reaction_rate  # Acetic acid 
    dcdt[3] =  (total_flow/V) * (inlet_temp-C[3]) - H/(rho*cp) * reaction_rate # Temperature part
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

def data_extract(data_path):
    '''Extracts the initial conditions for a the reaction \n
    Data_Path = relative path to the csv document'''
    data_numpy = np.genfromtxt(data_path, delimiter=';', dtype=None, names=True, encoding=None) #built in numpy function to extract data

    #Get temperature
    elapsed_time, temp = temp_extract(data_numpy) 

    #Get AAH Flowrate
    elapsed_time_aah, aah_flowrate_vector = temp_extract(data_numpy, x="P120_Flow")

    #Get Water Flowrate
    elapsed_time_water, water_flowrate_vector = temp_extract(data_numpy, x='P100_Flow')

    initial_temperature = np.min(temp) # Minimum temp = ini temp
    aah_flowrate = np.median(aah_flowrate_vector) # better than the average because sometimes we press prime before the experiment starts
    water_flowrate = np.median(water_flowrate_vector) # the signal is also kinda noisy 
    return elapsed_time, temp, initial_temperature, aah_flowrate, water_flowrate


if __name__ == '__main__':
    data_22c = data_extract('Data\Data from trade\CSTR\experiment14.10.csv')

    sol_time, sol_y = CSTR_model(26.9, 29.99, 185.83, 14.89, V=567)

    plt.figure(figsize=(10, 6))  # Larger figure size for better visibility

    # Plot simulation results (CSTR model)
    plt.plot(sol_time/60, sol_y[3, :] - 273.15, label='Predicted (Model)', linestyle='--', color='b', linewidth=2)

    # Plot experimental data
    plt.plot(data_22c[0], data_22c[1], label='Experimental Data', linestyle='-', color='r', linewidth=2)

    # Enhancing the plot aesthetics
    plt.xlabel('Time [minutes]', fontsize=14, weight='bold')
    plt.ylabel('Temperature [°C]', fontsize=14, weight='bold')
    plt.title('CSTR Temperature Profile: Model vs Experimental Data (Step Change)', fontsize=16, weight='bold')
    
    # Adding grid for clarity
    plt.grid(True, which='both', linestyle='--', linewidth=0.7)

    # Adjusting legend for clarity
    plt.legend(fontsize=12, loc='best')

    # Limiting x-axis to the relevant time frame
    plt.xlim(0, 55)

    # Applying tight layout for better spacing of elements
    plt.tight_layout()

    # Displaying the plot
    plt.show()