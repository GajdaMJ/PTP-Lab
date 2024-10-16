import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import scipy.integrate
from scipy.optimize import curve_fit


cstr_data = np.genfromtxt('Data/Data from trade/CSTR/experiment14.10.csv', delimiter=';', dtype=None, names=True, encoding='ISO-8859-1')  # Using ISO-8859-1 encoding

def data_extract(data, x, offset=0):
    rows = data[data['TagName'] == x]
    # Remove invalid rows
    valid_rows = [row for row in rows if row['vValue'] not in ['(null)', None]]
    
    # Parse DateTime for valid rows
    date_times = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_rows]
    vvalues = [float(row['vValue'])+offset for row in valid_rows]
    
    # Calculate elapsed time in minutes
    start_time = date_times[7]
    elapsed_time = [(dt - start_time).total_seconds() / 60 for dt in date_times]

    return elapsed_time, vvalues

######## Extracting experimental conductivity values as steady state #######
t, conductivity = np.array(data_extract(cstr_data, "Q210_PV"))
cool_t = np.array(data_extract(cstr_data, "T400_PV")[1])

cond_35 = np.max(conductivity)
cond_30 = np.mean(conductivity[80:86])
cond_27 = np.mean(conductivity[43:49])

conductivity = [cond_27,cond_30,cond_35]

######## Creating calibration curve for concentration conversion ########
cal_cond = [1695,1636,1594,1530,1429,1274,963,690,523] #conductivity measures
cal_conc = [1.74, 1.566, 1.392, 1.218, 1.044,0.87,0.435,0.2175,0.10875] # known acetic acid concentration [mol/L]

def model_func(x, a, b,d):
    return a * np.exp(b * x) + d

popt, pcov = curve_fit(model_func, cal_cond, cal_conc, p0=(1, -0.001,0))
a, b, d = popt

cal_cond_fit = np.linspace(min(cal_cond), max(cal_cond), 100)
cal_conc_fit = model_func(cal_cond_fit, a, b,d)

plt.scatter(cal_cond, cal_conc, color='red', label='Data Points')
plt.plot(cal_cond_fit, cal_conc_fit, label=f'Fit: {a:.5f} * exp({b:.5f} * x)', color='blue')

plt.xlabel('Conductivity [us/cm]')
plt.ylabel('Concentration [mol/L]')
plt.legend()
plt.title('CSTR Calibration Curve')
plt.show()

print(f"Fitted equation: y = {a:.5f} * exp({b:.5f}*x + {d:.5f})")


######### Converting experimental conductivity into concentration #########
conc_27 = model_func(cond_27, a, b,d)
conc_30 = model_func(cond_30, a, b,d)
conc_35 = model_func(cond_35, a, b,d)

concentrations = [conc_27,conc_30,conc_35] #mol/L


# ##### Inshallah Concentration ####

#water
mm_water = 18.01528 # (g/mol)
rho_water = 999.842 # (g/L)
cw_pure = rho_water/mm_water # (mol/L)

#Acetic acid
mm_AAH = 102.089 # (g/mol)
rho_AAH = 1082 # (g/L)
caah_pure = rho_AAH/mm_AAH # (mol/L)

v = 593.66 #ml
v_w = 1.745/60 #L/s
v_aah = 14/60 #L/s
v_f = v_w + v_aah #L/s
c_w0 = cw_pure * v_w/v_f #mol/L
c_aah0 = caah_pure *v_aah/v_f #mol/L

def k_eq(c_aa):
    return 2*c_aa *v_f /(v*(2*c_w0-c_aa)*(2*c_aah0-c_aa))   

k_val = []

for i in range(0,3):
    k_val.append(k_eq(concentrations[i]))


rep_temp = [1/(27+273.15),1/(30+273.15),1/(35+273.15)]

ln_k = np.log(k_val)

def lin_func(x, a, b):
    return a * x + b

popt1, pcov1 = curve_fit(lin_func, rep_temp, ln_k, p0=(1, 0))

# Extract the optimal parameters a and b
a1, b1 = popt1

# Generate data for plotting the fitted curve
rep_temp_fit = np.linspace(0, max(rep_temp), 100)
ln_k_fit = lin_func(rep_temp_fit, a1, b1)

#plt.scatter(rep_temp, ln_k, color='red', label='Data Points')
plt.plot(rep_temp_fit,ln_k_fit,label=f'Fit: {a1:.5f} * x +({b1:.5f})')
plt.xlabel('1/T [1/K]')
plt.ylabel('ln(k)')
plt.title("Fitted Linear Regression of 1/T vs ln(k)")
plt.xlim(0,np.max(rep_temp))
plt.legend()
plt.show()

print(f"Fitted equation: y = {a1:.5f} * x + ({b1:.5f})")
print(f"k0 = {np.exp(lin_func(0,a1,b1)):.5f}")
print(f"Ea = {(a1*-8.3145)}")

rep_temp_fit_1 = np.linspace(min(rep_temp), max(rep_temp), 100)
ln_k_fit_1 = lin_func(rep_temp_fit_1, a1, b1)

plt.scatter(rep_temp, ln_k, color='red', label='Data Points')
plt.plot(rep_temp_fit_1, ln_k_fit_1, label=f'Fit: {a1:.5f} * x +({b1:.5f})', color='blue')
plt.xlabel('1/T [1/K]')
plt.ylabel('ln(k)')
plt.title("Plot to determine k0 and Ea")
plt.legend()
plt.show()

