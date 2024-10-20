from CSTR_Model import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
data_27c = data_extract('Data\\CSTR_Data\\CSTR 27c.csv')
sol_27 = CSTR_model(data_27c[2], data_27c[4], data_27c[3], V=600)

data_33c = data_extract('Data\CSTR_Data\\42c 25.09.csv')
sol_33 = CSTR_model(data_33c[2], data_33c[4], data_33c[3], V=567)

data_22c = data_extract('Data\\CSTR_Data\\23.09 22c.csv')
sol_22 = CSTR_model(data_22c[2], data_22c[4], data_22c[3], V=567)

data_37c = data_extract('Data\\CSTR_Data\\25.09 37c.csv')
sol_37 = CSTR_model(data_37c[2], data_37c[4], data_37c[3], V=507.11)



fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot temperature vs. time on the left subplot (ax1)
ax1.plot(sol_22.t / 60, sol_22.y[3, :] - 273.15, label='Model prediction 22°C', color='tab:blue', linewidth=2)
ax1.plot(data_22c[0], data_22c[1], 'o', label='Data 22°C', color='tab:blue', markersize=5)

ax1.plot(sol_27.t / 60, sol_27.y[3, :] - 273.15, label='Model prediction 27°C', color='tab:orange', linewidth=2)
ax1.plot(np.array(data_27c[0])-1.33, data_27c[1], 'o', label='Data 27°C', color='tab:orange', markersize=5)

ax1.plot(sol_33.t / 60, sol_33.y[3, :] - 273.15, label='Model prediction 33°C', color='tab:green', linewidth=2)
ax1.plot(np.array(data_33c[0])-2.54, data_33c[1], 'o', label='Data 33°C', color='tab:green', markersize=5)

ax1.plot(sol_37.t / 60, sol_37.y[3, :] - 273.15, label='Model prediction 37°C', color='tab:red', linewidth=2)
ax1.plot(data_37c[0], data_37c[1], 'o', label='Data 37°C', color='tab:red', markersize=5)


# Customize the first plot (temperature vs. time)
ax1.set_xlabel('Time (minutes)', fontsize=14, weight='bold')
ax1.set_ylabel('Temperature (°C)', fontsize=14, weight='bold')
ax1.set_title('Temperature vs. Time', fontsize=16, weight='bold')
ax1.legend(fontsize=12)
ax1.grid(True, linestyle='--', alpha=0.7)
ax1.set_xlim(0,33)

interpolator27 = scipy.interpolate.interp1d(sol_27.t/60, sol_27.y[3, :] - 273.15, kind='linear')
interpolator22 = scipy.interpolate.interp1d(sol_22.t/60, sol_22.y[3, :] - 273.15, kind='linear')
interpolator33 = scipy.interpolate.interp1d(sol_33.t/60, sol_33.y[3, :] - 273.15, kind='linear')
interpolator37 = scipy.interpolate.interp1d(sol_37.t/60, sol_37.y[3, :] - 273.15, kind='linear')

time_22c = np.array(data_22c[0][15:])
temp_22c = np.array(data_22c[1][15:])
sol_interpolated_22c = interpolator22(time_22c)  # Interpolate model at specific times
residuals_22c = temp_22c - sol_interpolated_22c  # Calculate residuals
ax2.plot(time_22c, residuals_22c, label=f'Residual, Probe 22°C', linewidth=2)

# For 27°C data (already done, shown for consistency)
time_27c = np.array(data_27c[0][4:])
temp_27c = np.array(data_27c[1][4:])
sol_interpolated_27c = interpolator27(time_27c)
residuals_27c = temp_27c - sol_interpolated_27c
ax2.plot(time_27c, residuals_27c, label=f'Residual, Probe 27°C', linewidth=2)

# For 33°C data (starting from index 5)
time_33c = np.array(data_33c[0][5:])
temp_33c = np.array(data_33c[1][5:])
sol_interpolated_33c = interpolator33(time_33c)
residuals_33c = temp_33c - sol_interpolated_33c
ax2.plot(time_33c, residuals_33c, label=f'Residual, Probe 33°C', linewidth=2)

# For 37°C data (starting from index 15)
time_37c = np.array(data_37c[0][15:])
temp_37c = np.array(data_37c[1][15:])
sol_interpolated_37c = interpolator37(time_37c)
residuals_37c = temp_37c - sol_interpolated_37c
ax2.plot(time_37c, residuals_37c, label=f'Residual, Probe 37°C', linewidth=2)

delta_T_22c = temp_22c[-1] - temp_22c[0]
delta_T_27c = temp_27c[-1] - temp_27c[0]
delta_T_33c = temp_33c[-1] - temp_33c[0]
delta_T_37c = temp_37c[-1] - temp_37c[0]
ax1.text(34, temp_22c[-1], f'ΔT = {delta_T_22c:.2f}°C', color='tab:blue', fontsize=12)
ax1.text(34, temp_27c[-1], f'ΔT = {delta_T_27c:.2f}°C', color='tab:orange', fontsize=12)
ax1.text(34, temp_33c[-1], f'ΔT = {delta_T_33c:.2f}°C', color='tab:green', fontsize=12)
ax1.text(34, temp_37c[-1], f'ΔT = {delta_T_37c:.2f}°C', color='tab:red', fontsize=12)



# Plot the residuals (difference between real data and model) on the right subplot
# ax2.plot(sol_mereal.t / 60, conc_data-sol_mereal.y[3,:]+273.15, label='Residuals Volume 567 ml', color='tab:blue', linewidth=2)
# ax2.plot(sol_me500.t / 60, conc_data-sol_me500.y[3,:]+273.15, label='Residuals Volume 500 ml', color='tab:orange', linestyle='--', linewidth=2)
# ax2.plot(sol_me600.t / 60, conc_data-sol_me600.y[3,:]+273.15, label='Residuals Volume 600 ml', color='tab:green', linestyle='-.', linewidth=2)

# Customize the second plot (residuals)
ax2.set_xlabel('Time (minutes)', fontsize=14, weight='bold')
ax2.set_ylabel('Residuals (°C)', fontsize=14, weight='bold')
ax2.set_title('Residuals (Model vs. Real Data)', fontsize=16, weight='bold')
ax2.legend(fontsize=12)
ax2.grid(True, linestyle='--', alpha=0.7)
ax2.axhline(0, color='black', linestyle='--', linewidth=1)
ax2.set_xlim(0,33)
# Tight layout for better spacing
plt.tight_layout()

# Show the combined plot
plt.show()