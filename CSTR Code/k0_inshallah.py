import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import scipy.integrate
from scipy.optimize import curve_fit
from matplotlib.ticker import ScalarFormatter


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

cond_35 = np.max(conductivity)+50
cond_30 = np.mean(conductivity[80:86])
cond_27 = np.mean(conductivity[43:49])

cond_27_min = np.min(conductivity[43:49])
cond_27_max = np.max(conductivity[43:49])
cond_30_min = np.min(conductivity[80:86])
cond_30_max = np.max(conductivity[80:86])
conductivity = [cond_27,cond_30]


######## Creating calibration curve for concentration conversion ########
cal_cond = [1695,1636,1594,1530,1429,1274,963,690,523] #conductivity measures
#cal_cond = [1695,1636,1586,1531,1419,1139,963,682,506] ##best?
# cal_cond = [1695,1636,1586,1534,1429,1274,966,690,506] 
cal_conc = [1.74, 1.566, 1.392, 1.218, 1.044,0.87,0.435,0.2175,0.10875] # known acetic acid concentration [mol/L]

def model_func(x, a, b,d):
    return a * np.exp(b * x) + d

popt, pcov = curve_fit(model_func, cal_cond, cal_conc, p0=(1, -0.001,0))
a, b, d = popt

cal_cond_fit = np.linspace(min(cal_cond), max(cal_cond), 100)
cal_conc_fit = model_func(cal_cond_fit, a, b,d)

plt.scatter(cal_cond, cal_conc, color='red', label='Data Points')
plt.plot(cal_cond_fit, cal_conc_fit, label=f'Fit: {a:.5f} * exp({b:.5f} * x) + {d:.5f}', color='blue')

plt.xlabel('Conductivity [us/cm]')
plt.ylabel('Concentration [mol/L]')
plt.legend()
plt.title('CSTR Conductivity vs Concentration Calibration Curve')
plt.show()

print(f"Fitted equation: y = {a:.5f} * exp({b:.5f}*x) + {d:.5f}")


predict = []
for i in range (0,9):
    predict.append(model_func(cal_cond[i], a, b,d))
corr_matrix = np.corrcoef(cal_conc, predict)
corr = corr_matrix[0,1]
R_sq = corr**2
 
print(R_sq)

######### Converting experimental conductivity into concentration #########
conc_27 = model_func(cond_27, a, b,d) # mol/L 
conc_30 = model_func(cond_30, a, b,d) # mol/L
conc_35 = model_func(cond_35, a, b,d) # mol/L

#range
conc_27_1 = model_func(cond_27_min,a,b,d)
conc_27_2 = model_func(cond_27_max,a,b,d)
conc_30_1 = model_func(cond_30_min,a,b,d)
conc_30_2 = model_func(cond_30_max,a,b,d)

range_conc = [conc_27_1*1e-3,conc_27_2*1e-3,conc_30_1*1e-3,conc_30_2*1e-3]

concentrations = [conc_27*1e-3,conc_30*1e-3,conc_35*1e-3] #mol/ml
print(f'concentrations = {(concentrations)}')
print(f'concentrations1 = {(range_conc)}')

# ##### Inshallah Concentration ####

#water
mm_water = 18.01528 # (g/mol)
rho_water = .999842 # (g/mL)
cw_pure = rho_water/mm_water # (mol/mL)

#Acetic acid
mm_AAH = 102.089 # (g/mol)
rho_AAH = 1.082 # (g/mL)
caah_pure = rho_AAH/mm_AAH # (mol/mL)

v = 593.66 #mL
v_w = 174.5/60 #mL/s
v_aah = 14/60 #mL/s 
v_f = v_w + v_aah #mL/s
c_w0 = cw_pure * v_w/v_f #mol/mL
c_aah0 = caah_pure *v_aah/v_f #mol/mL

print(c_w0,c_aah0)
def k_eq(c_aa):
    return 2*c_aa *v_f /(v*(2*c_w0-c_aa)*(2*c_aah0-c_aa))   

k_val = []

for i in range(0,2):
    k_val.append(k_eq(concentrations[i]))

k_val_range = []
for i in range(0,4):
    k_val_range.append(k_eq(range_conc[i]))

print(f'kval = {(k_val)}')
rep_temp = [1/(27+273),1/(30+273) ]
rep_temp_range = [1/(27+273),1/(27+273),1/(30+273),1/(30+273) ]


ln_k = np.log(k_val)
ln_k_range = np.log(k_val_range)


print(f'lnkrange = {(ln_k_range)}')

def lin_func(x, a, b):
    return a * x + b

popt1, pcov1 = curve_fit(lin_func, rep_temp, ln_k, p0=(1, 0))

# Extract the optimal parameters a and b
a1, b1 = popt1

# Generate data for plotting the fitted curve
rep_temp_fit = np.linspace(0, max(rep_temp), 100)
ln_k_fit = lin_func(rep_temp_fit, a1, b1)

#####
popt2, pcov2 = curve_fit(lin_func, rep_temp, [ln_k_range[0],ln_k_range[3]], p0=(1, 0))

# Extract the optimal parameters a and b
a2, b2 = popt2

# Generate data for plotting the fitted curve
ln_k_fit_1 = lin_func(rep_temp_fit, a2, b2)

######
popt3, pcov3 = curve_fit(lin_func, rep_temp, [ln_k_range[1],ln_k_range[2]], p0=(1, 0))

# Extract the optimal parameters a and b
a3, b3 = popt3

# Generate data for plotting the fitted curve
ln_k_fit_2 = lin_func(rep_temp_fit, a3, b3)
#####

#plt.scatter(rep_temp, ln_k, color='red', label='Data Points')
plt.plot(rep_temp_fit,ln_k_fit, color = 'red', label=f'Fit: {a1:.5f} * x +({b1:.5f})')
plt.plot(rep_temp_fit,ln_k_fit_1, color = 'blue', linestyle = 'dashed', label=f'Fit: {a2:.5f} * x +({b2:.5f})')
plt.plot(rep_temp_fit,ln_k_fit_1, color = 'green', linestyle = 'dashed', label=f'Fit: {a3:.5f} * x +({b3:.5f})')
# plt.plot(rep_temp_fit, lin_func(rep_temp_fit, 9.825e4/-8.3145, np.log(4.4e14)))
plt.xlabel('1/T [1/K]')
plt.ylabel('ln(k)')
plt.title("Fitted Linear Regression of 1/T vs ln(k)")
plt.xlim(0,np.max(rep_temp))
plt.legend()
plt.show()

print(f"Fitted equation: y = {a1:.5f} * x + ({b1:.5f})")
print(f"k0 = {np.exp(lin_func(0,a1,b1)):.8f}")
print(f"Ea = {(a1*-8.3145)}")

print(f"Fitted equation1: y = {a2:.5f} * x + ({b2:.5f})")
print(f"k0 = {np.exp(lin_func(0,a2,b2)):.8f}")
print(f"Ea = {(a2*-8.3145)}")

print(f"Fitted equation2: y = {a3:.5f} * x + ({b3:.5f})")
print(f"k0 = {np.exp(lin_func(0,a3,b3)):.8f}")
print(f"Ea = {(a3*-8.3145)}")

print(f'k0 avg = {(np.exp(lin_func(0,a1,b1))+np.exp(lin_func(0,a2,b2))+np.exp(lin_func(0,a3,b3)))/3:.8f}')
print(f'Ea avg = {((a1*-8.3145)+(a2*-8.3145)+(a3*-8.3145))/3:.8f}')

rep_temp_fit_1 = np.linspace(min(rep_temp), max(rep_temp), 100)
ln_k_fit_a = lin_func(rep_temp_fit_1, a1, b1)
ln_k_fit_a1 = lin_func(rep_temp_fit_1, a2, b2)
ln_k_fit_a2 = lin_func(rep_temp_fit_1, a3, b3)

plt.scatter(rep_temp, ln_k, color='red', label='Average Data Points')
plt.scatter(rep_temp_range, ln_k_range, color='grey', label='Maximum and Minimum Data Points')
# plt.plot(rep_temp_fit_1, lin_func(rep_temp_fit_1, 9.825e4/-8.3145, np.log(4.4e14)))
plt.plot(rep_temp_fit_1, ln_k_fit_a, label=f'Fit: {a1:.5f} * x +({b1:.5f})', color='red')
plt.plot(rep_temp_fit_1, ln_k_fit_a1, label=f'Fit: {a2:.5f} * x +({b2:.5f})', color='blue',linestyle = 'dashed')
plt.plot(rep_temp_fit_1, ln_k_fit_a2, label=f'Fit: {a3:.5f} * x +({b3:.5f})', color='green',linestyle = 'dashed')
plt.gca().xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.xlabel('1/T [1/K]')
plt.ylabel('ln(k)')
plt.title("Plot to Determine k0 and Ea")
plt.legend()
plt.show()

