import numpy as np
import matplotlib.pyplot as plt

x_end = 1 # meter
tend = 4
d = 1e-3
kr = 1.2 #per second

# mew = volumetric flow / cross sectional area
diam = 2.5 # meter
flow = 6.3e-4 #m^3/s
mew = flow / (np.pi*(diam**2)/4)

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

def model_2_pbr(T,fv1,fv2,V=137,tspan =[0,3600],n=6):
    '''Models the behavior of the reaction: Water + Acetic Anhydride -> 2 * Acetic acid in an adiabatic CSTR reactor. \n
    Required Arguments: \n
    T = inlet temperature for the reactor given in units celsius \n
    fv1 = flow rate of water in units ml/min \n
    fv2 = flow rate of acetic anhydride ml/min \n
    Optional Arguments: \n
    V = volume of the reactor in units ml (default set to 500ml) \n
    tspan = list of evaluation time in units seconds (default set to [0,3600]) \n
    n = number of reactors \n
    This function was built for the course "Practical Process Technology (6P4X0)" 
    '''
    #Chemical constants
    mm_water = 18.01528 # (g/mol)
    rho_water = 0.999842 # (g/ml)
    cw_pure = rho_water/mm_water # (mol/ml)
    mm_AAH = 102.089 # (g/mol)
    rho_AAH = 1.082 # (g/ml)
    caah_pure = rho_AAH/mm_AAH # (mol/ml)

   # Convert flow rates (ml/min to ml/s)
    fv_w_dm3_s = fv1 / 60  # Water flow rate in ml/s
    fv_a_dm3_s = fv2  / 60  # Anhydride flow rate in ml/s

    flow_array = [fv_w_dm3_s, fv_a_dm3_s]
    
    v_pfr = V/n
    #Thermodynamic constants
    params = {
        "C_in_water": (flow_array[0]*cw_pure)/(flow_array[0]+flow_array[1]),
        "C_in_AAH": (flow_array[1]*caah_pure)/(flow_array[0]+flow_array[1]),
        "Inlet temperature": T+273.15,
        "flow": flow_array,
        "V": v_pfr,  # Volume in ml
        "k0": 6e5,          # Reaction rate constant (ml/mol/s)

        # Thermodynamic constants (taken from Asprey et al., 1996)
        "Ea": 45622.34,             # Activation energy (J/mol)
        "R": 8.314,              # Gas constant (J/mol/K)
        "H": -56.6e3,              # Enthalpy change (J/mol)
        "rho": 1,            # Density (g/ml)
        "cp": 4.186             # Heat capacity (J/g/K)
    }
    x_cond = [cw_pure,0,0,T+273.15] # Initial Conditions 
    xini = np.zeros(4*n)
    for i in range(4*n):
        if i == 0: 
            xini[i] = x_cond[0]
        if i == 1:
            xini[i] = x_cond[1]
        if i == 2:
            xini[i] = x_cond[2]
        if i == 3:
            xini[i] = x_cond[3]
        elif np.mod(i,4)==0:
            xini[i] = x_cond[0]
        elif np.mod(i,4)==1:
            xini[i] = x_cond[1]
        elif np.mod(i,4)==2:
            xini[i] = x_cond[2]
        elif np.mod(i,4)==3:
            xini[i] = x_cond[3]
            

    sol_me = master_function(lambda t, C: der_func(t, C, params,n), tspan, xini, method='rk4', number_of_points=300) #Master function is a differential equation solver made for Numerical Methods.
    return sol_me


def der_func(t,C, parameters,n):
    '''This function contains the differential equations to solve the reaction A+B->2C in an adiabatic 
    CSTR. \n
    t=time (seconds) \n
    c = Concentration vector like [c_water, c_AAH, c_AA, Temperature]\n
    parameters = dictionary containing thermodynamic constants
    '''
    #initialize the solution matrix
    dcdt = np.zeros(4*n)

    # Getting parameters
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

    total_flow = flow[0]+flow[1]
    #solving for each reactor
    for i in range(4*n):
        #for the first reactor it takes the initial conditions
        if i == 0:
            dcdt[i] = (total_flow/V)*(C_in_w - C[0]) - C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3]))
        if i == 1:
            dcdt[i] = (total_flow/V)*(C_in_AAH - C[1])  - C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3]))
        if i == 2:
            dcdt[i] = (total_flow/V)*(0 - C[2]) + 2*C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3])) 
        if i == 3:
            dcdt[i] = (total_flow/V) * (inlet_temp-C[3]) - H/(rho*cp) * C[0]*C[1] * k0 * np.exp(-Ea/(R*C[3]))
        #make a condition for the remain reactors which will not be taking the initial conditions but the outlet of the reactor before
        
        elif np.mod(i,4)==0: #np.mod(divident,divisor) with the output being the remainder
            dcdt[i] = (total_flow/V)*(C[i-4] - C[i]) - C[i]*C[i+1] * k0 * np.exp(-Ea/(R*C[i+3]))
        
        elif np.mod(i,4) == 1:
            dcdt[i] = (total_flow/V)*(C[i-4] - C[i])  - C[i-1]*C[i] * k0 * np.exp(-Ea/(R*C[i+2]))

        elif np.mod(i,4) == 2:
            dcdt[i] = (total_flow/V)*(C[i-4] - C[i]) + 2*C[i-2]*C[i-1] * k0 * np.exp(-Ea/(R*C[i+1]))
        
        elif np.mod(i,4) == 3:
            dcdt[i] = (total_flow/V) * (C[i-4]-C[i]) - H/(rho*cp) * C[i-3]*C[i-2] * k0 * np.exp(-Ea/(R*C[i]))
    
    return dcdt

sol_me = model_2_pbr(20,100,20,V=137,tspan =[0,3600],n=100)

print(sol_me)

plt.plot(sol_me[0],sol_me[1])
plt.show()