import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.optimize import curve_fit

# Load data for 27°C
my_data = np.genfromtxt('Data\\CSTR\\Runs 16.09\\CSTR 27c.csv', delimiter=';', dtype=None, names=True, encoding=None)

# Extract temperature data
def temp_extract(data):
    t200_pv_rows = data[data['TagName'] == 'T200_PV']
    valid_rows = [row for row in t200_pv_rows if row['vValue'] not in ['(null)', None]]
    
    # Convert date strings to datetime objects
    date_times = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_rows]
    vvalues = [float(row['vValue']) for row in valid_rows]
    
    start_time = date_times[0]
    
    # Calculate elapsed time in minutes
    elapsed_time = [(dt - start_time).total_seconds() / 60 for dt in date_times]
    
    return elapsed_time, vvalues

elapsed_time_27c, temp_27c = temp_extract(my_data)

# Filter data for 4 < t < 32.5 minutes
time_mask = (np.array(elapsed_time_27c) > 4) & (np.array(elapsed_time_27c) < 32.5)
filtered_time = np.array(elapsed_time_27c)[time_mask]
filtered_temp = np.array(temp_27c)[time_mask]

# Define the model function for temperature as a function of time
def temperature_model(t, T0, Tf, k):
    return Tf + (T0 - Tf) * np.exp(-k * t)

# Initial guess for parameters [T0, Tf, k]
initial_guess = [filtered_temp[0], 27 + 273.15, 0.01]

# Fit the model to the filtered data
popt, pcov = curve_fit(temperature_model, filtered_time, filtered_temp, p0=initial_guess)

# Extract fitted parameters
T0_fit, Tf_fit, k_fit = popt
print(f"Fitted T0: {T0_fit:.2f} K, Fitted Tf: {Tf_fit:.2f} K, Fitted k: {k_fit:.5f} min^-1")

# Generate fitted data for plotting
fitted_temps = temperature_model(np.array(filtered_time), *popt)

# Compute dT/dt using numpy's gradient
dT_dt = np.gradient(fitted_temps, filtered_time)

# Plot dT/dt vs T(t)
plt.figure(figsize=(10, 6))
plt.plot(fitted_temps, dT_dt, marker='o', linestyle='-', color='blue')
plt.xlabel('Temperature (K)')
plt.ylabel('dT/dt (K/min)')
plt.title('Rate of Change of Temperature vs Temperature at 27°C')
plt.grid()
plt.show()
