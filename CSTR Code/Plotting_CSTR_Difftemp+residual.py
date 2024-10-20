from Model_3 import *
import matplotlib.pyplot as plt
import numpy as np

data_27c = data_extract('Data\\CSTR\\Runs 16.09\\CSTR 27c.csv')
sol_27 = CSTR_model(data_27c[2], data_27c[4], data_27c[3], V=600)

data_33c = data_extract('Data/CSTR/42c 25.09.csv')
sol_33 = CSTR_model(data_33c[2], data_33c[4], data_33c[3], V=567)

data_22c = data_extract('Data\\CSTR\\23.09 22c.csv')
sol_22 = CSTR_model(data_22c[2], data_22c[4], data_22c[3], V=567)

data_37c = data_extract('Data\\CSTR\\25.09 37c.csv')
sol_37 = CSTR_model(data_37c[2], data_37c[4], data_37c[3], V=507.11)

sol_time_minutes = sol_time / 60 #Sol time in seconds

# Make data np arrays
exp_time = np.array(data_22c[0])
exp_temp = np.array(data_22c[1])

# Interpolate the model predictions to match experimental time points
model_temp_interp = np.interp(exp_time, sol_time_minutes, sol_y[3, :] - 273.15)

# Calculate residuals (experimental - model)
residuals = exp_temp - model_temp_interp

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))  # Create two plots side by side

# original model vs experimental data
ax1.plot(sol_time_minutes, sol_y[3, :] - 273.15, label='Predicted (Model)', linestyle='--', color='b', linewidth=2)
ax1.plot(exp_time, exp_temp, label='Experimental Data', linestyle='-', color='r', linewidth=2)

# Left plot vibes
ax1.set_xlabel('Time [minutes]', fontsize=14, weight='bold')
ax1.set_ylabel('Temperature [°C]', fontsize=14, weight='bold')
ax1.set_title('CSTR Temperature Profile: Model vs Experimental Data', fontsize=16, weight='bold')
ax1.grid(True, which='both', linestyle='--', linewidth=0.7)
ax1.legend(fontsize=12, loc='best')
ax1.set_xlim(0, 55)

# Subtitle
ax1.text(0.5, 0.92, r'$T_1 = 27°C, T_2 = 30°C$', transform=ax1.transAxes, fontsize=12, ha='center')
# Plot on the right: residuals (experimental - model)
ax2.plot(exp_time, residuals, label='Residuals (Exp - Model)', color='g', linestyle='-', linewidth=2)

# Right plot vibes
ax2.set_xlabel('Time [minutes]', fontsize=14, weight='bold')
ax2.set_ylabel('Residuals [°C]', fontsize=14, weight='bold')
ax2.set_title('Residuals of Temperature (Exp - Model)', fontsize=16, weight='bold')
ax2.axhline(0, color='black', linestyle='--', linewidth=1)  # Add a horizontal line at 0 for reference
ax2.grid(True, which='both', linestyle='--', linewidth=0.7)
ax2.set_xlim(0, 55)


plt.tight_layout()
plt.show()