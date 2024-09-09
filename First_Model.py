import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

# Assume reaction is 1st order wrt both components
# Assume isothermal (no exotherm)
# Assume constant density


fv_w = 1 # Water Flow Rate (ml/min)
fv_a = 1 # Anhydride Flow Rate (ml/min)
t_cstr = 25 # CSTR Temperature (c)
k_1 = 0.02 # Reaction Kinetics (dm3/mol/s)
v_cstr = 200 # Volume of CSTR (ml)

flow = [fv_w,fv_a]
# w + AAH -> 2AA 

# -ra = k_1 * C_A * C_B


def der_func(t,C, phiv = flow, Cain=0.5, Cbin = 0.8, k1=0.1, V = 1):
    dcdt = np.zeros(3)
    dcdt[0] = (phiv[0]/V) * (Cain-C[0]) - k1 * C[0] *C[1] #Water Concentration derv
    dcdt[1] = (phiv[1]/V) * (Cbin- C[1]) - k1 * C[0] *C[1] #Anhydride Concentration derv
    dcdt[2] = ((phiv[0]+phiv[1])/V) * (-C[2]) + 2 * k1 * C[0] *C[1] #Acetic acid 
    return dcdt

tspan = [0,60]
xini = [0,0,0] 

sol = scipy.integrate.solve_ivp(der_func, tspan, xini)

plt.plot(sol.t, sol.y[0], label='Ca')
plt.plot(sol.t, sol.y[1], label='Cb')
plt.plot(sol.t, sol.y[2], label='Cc')
plt.xlabel('time')
plt.ylabel('Concentration(mol/L)')
plt.legend()
plt.title('Concentration of various components in a CSTR')
plt.show()
