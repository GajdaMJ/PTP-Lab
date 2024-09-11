import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

# Assume reaction is 1st order wrt both components
# Assume isothermal (no exotherm)
# Assume constant density


fv_w = 161.10420227050781 # Water Flow Rate (ml/min)
fv_a = 11.524953842163086 # Anhydride Flow Rate (ml/min)
t_cstr = 25 # CSTR Temperature (c)
k_1 = 0.02 # Reaction Kinetics (dm3/mol/s)
v_cstr = 500 # Volume of CSTR (ml)

flow_array = [fv_w,fv_a]

def der_func(t,C, ini_cond, flow=flow_array, V=1, k0=1, Ea=1, R=1, H=1, rho=1, cp=1):
    dcdt = np.zeros(4)
    reaction_rate = C[0]*C[1] * k0  * np.exp(-Ea/R*C[3])
    total_flow = flow[0]+flow[1]
    dcdt[0] =  (flow[0]/V)*(ini_cond[0] - C[0]) -reaction_rate #Water Concentration derv
    dcdt[1] =  (flow[1]/V)*(ini_cond[1] - C[1]) - reaction_rate#Anhydride Concentration derv
    dcdt[2] =  (total_flow/V)*(ini_cond[2] - C[2]) + 2*reaction_rate#Acetic acid 
    dcdt[3] =  (total_flow/V) * (ini_cond[3]-C[3]) - H/(rho*cp) * reaction_rate#Temperature part
    return dcdt

tspan = [0,60]
xini = [0,0,0,0] 

sol = scipy.integrate.solve_ivp(der_func, tspan, xini, args=(xini,))

plt.plot(sol.t, sol.y[0], label='Ca')
plt.plot(sol.t, sol.y[1], label='Cb')
plt.plot(sol.t, sol.y[2], label='Cc')
plt.xlabel('time')
plt.ylabel('Concentration(mol/L)')
plt.legend()
plt.title('Concentration of various components in a CSTR')
plt.show()
