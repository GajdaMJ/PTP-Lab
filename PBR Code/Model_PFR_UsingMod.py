import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
# Assume reaction is 1st order wrt both components
# Assume isothermal (no exotherm)
# Assume constant density

def master_function(fun,tspan, y0, method='rk4', number_of_points=100):
    '''General function to solve system of differential equations. Does not work on single differential equations. \n
    fun = function 
    y0 = vector of initial conditions
    optional:\n
    method = You can select the method with which your system of differential equations will be evaluated. Default set to second order Runge-Kutta. \n
    Supported methods : midpoint method ('midpoint'), euler method ('euler'), Classical second order Runge-Kutta ('rk2'), classical fourth order Runge-Kutta ('rk4').
    number_of_points = how many steps. Default set to 100. Increasing this reduces error but increases computation time. '''
    dt = (tspan[1] - tspan[0])/number_of_points
    t = np.linspace(tspan[0], tspan[1], number_of_points+1)
    y = np.zeros((number_of_points+1, len(y0))) # len(y0) because you would need an initial condition for each derivative.
    for i in range(len(y0)): #initial conditions as a loop to ensure universability.
        y[0,i] = y0[i]
    if method == 'midpoint':
        for i in range(number_of_points):
            k1 = fun(t[i], y[i,:])
            k2 = fun(t[i] + dt*0.5, y[i,:] + 0.5*dt*k1)
            y[i+1,:] = y[i,:] + dt * k2
    elif method == 'euler':
        for i in range(number_of_points):
            y[i+1,:] = y[i,:] + dt * fun(t[i], y[i,:])
    elif method == 'rk2':
        for i in range(number_of_points):
            k1 = fun(t[i], y[i,:])
            k2 = fun(t[i] + dt, y[i] + dt*k1)
            y[i+1,:] = y[i] + dt*0.5*(k1+k2)
    elif method == 'rk4':
        for i in range(number_of_points):
            k1 = fun(t[i], y[i,:])
            k2 = fun(t[i] + dt*0.5, y[i,:] + 0.5*dt*k1)
            k3 = fun(t[i] + dt*0.5, y[i,:] + 0.5*dt*k2)
            k4 = fun(t[i] +dt, y[i,:] + dt*k3)
            y[i+1,:] = y[i] + dt*((1/6)*k1 + (1/3)*(k2+k3) + (1/6)*k4)
    else:
        return 'Unknown method specified. Check documentation for supported methods' # In case an unknown method is specified
    return t, y

def CSTR_model(T,fv1,fv2, V=123, tspan = [0,3600], n=6):
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

    
    params = { # Stores the relevant thermodynamic constants as a dictionary 
        "C_in_water": (flow_array[0]*cw_pure)/(flow_array[0]+flow_array[1]),
        "C_in_AAH": (flow_array[1]*caah_pure)/(flow_array[0]+flow_array[1]),
        "Inlet temperature": T+273.15, # Temp but now in kelvin
        "flow": flow_array,
        "V": v_pfr_tank,  # Volume in ml
        "k0": 9.092246269480829e+16,          # Reaction rate constant (ml/mol/s)

        # Thermodynamic constants (taken from Asprey et al., 1996)
        "Ea": 104982.28755856157,             # Activation energy (J/mol)
        "R": 8.314,              # Gas constant (J/mol/K)
        "H": -56.6e3,              # Enthalpy change (J/mol)
        "rho": 1,            # Density (g/ml)
        "cp": 4.186             # Heat capacity (J/g/K)
    }
    # print(params['C_in_AAH']*params['C_in_water'])
    xini_temp = [cw_pure,0,0,T+273.15] # Initial Conditions 
    xini = np.zeros(4*n)
    for i in range(4*n):
        if np.mod(i,4)==0:
            xini[i] = xini_temp[0]
        elif np.mod(i,4)==1:
            xini[i] = xini_temp[1]
        elif np.mod(i,4)==2:
            xini[i] = xini_temp[2]
        elif np.mod(i,4)==3:
            xini[i] = xini_temp[3]

    sol_me = master_function(lambda t, C: der_func(t, C, params, n), tspan, xini, method='rk4', number_of_points=300) #Master function is a differential equation solver made for Numerical Methods.
    return sol_me

def der_func(t,C, parameters, n=6):
    '''This function contains the differential equations to solve the reaction A+B->2C in an adiabatic 
    CSTR. \n
    t=time (seconds) \n
    c = Concentration vector like [c_water, c_AAH, c_AA, Temperature]\n
    parameters = dictionary containing thermodynamic constants
    '''
    # Initializing derivative vector
    dcdt = np.zeros(4*n)
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
    for i in range(4*n):
        if i<4:
           dcdt[0] =  (total_flow/V)*(C_in_w - C[0])     - C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3])) # Water Concentration derv
           dcdt[1] =  (total_flow/V)*(C_in_AAH - C[1])   - C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3])) # Anhydride Concentration derv
           dcdt[2] =  (total_flow/V)*(0 - C[2])          + 2*C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3]))  # Acetic acid 
           dcdt[3] =  (total_flow/V) * (inlet_temp-C[3]) - H/(rho*cp) * C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3])) # Temperature part
        else:
            if np.mod(i,4)==0: #np.mod(divident,divisor) with the output being the remainder
                dcdt[i] = (total_flow/V)*(C[i-4] - C[i]) - C[i]*C[i+1] * k0 * np.exp(-Ea/(R*C[i+3]))
            elif np.mod(i,4) == 1:
                dcdt[i] = (total_flow/V)*(C[i-4] - C[i])  - C[i-1]*C[i] * k0 * np.exp(-Ea/(R*C[i+2]))
            elif np.mod(i,4) == 2:
                dcdt[i] = (total_flow/V)*(C[i-4] - C[i]) + 2*C[i-2]*C[i-1] * k0 * np.exp(-Ea/(R*C[i+1]))
            elif np.mod(i,4) == 3:
                dcdt[i] = (total_flow/V) * (C[i-4]-C[i]) - H/(rho*cp) * C[i-3]*C[i-2] * k0 * np.exp(-Ea/(R*C[i]))
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
    data_22c = data_extract('Data\\CSTR\\Runs 16.09\\CSTR 27c.csv')
    sol_me = CSTR_model(20,100,20,V=137,tspan =[0,3600],n=9)

    # plt.plot(sol_me[0], sol_me[1][:, 1], label='Conc. AAH_me')
    # plt.plot(sol_me[0], sol_me[1][:, 2], label='Conc. AA_me')
    plt.plot(sol_me[0]/60, sol_me[1][:, 9*4-1]-273.15, label='think')
    plt.plot(data_22c[0], data_22c[1], label='real')
    plt.xlabel('Time (minutes)')
    plt.xlim(0, np.max(data_22c[1]))
    plt.ylabel('Temperature')
    plt.legend()
    plt.title('Temperature')
    plt.show()
