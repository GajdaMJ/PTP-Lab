import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import scipy.integrate
from scipy.optimize import curve_fit



cal_cond = [1145,1095,977,735,539]
cal_conc = [1.74,0.87,0.435,0.2175,0.10875] # mol/L

def model_func(x, a, b,d):
    return a * np.exp(b * x) + d

popt, pcov = curve_fit(model_func, cal_cond, cal_conc, p0=(1, -0.001,0))
a, b, d = popt

cal_cond_fit = np.linspace(min(cal_cond), max(cal_cond), 100)
cal_conc_fit = model_func(cal_cond_fit, a, b,d)

plt.scatter(cal_cond, cal_conc, color='red', label='Data Points')
plt.plot(cal_cond_fit, cal_conc_fit, label=f'Fit: {a:.10f} * exp({b:.5f} * x)', color='blue')
plt.xlabel('Conductivity [us/cm]')
plt.ylabel('Concentration [mol/L]')
plt.legend()
plt.title('PBR Conductivity vs Concentration Calibration Curve')
plt.show()

print(f"Fitted equation: y = {a:.10f} * exp({b:.5f}*x + {d:.5f})")

