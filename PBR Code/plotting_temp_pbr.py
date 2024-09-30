import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from numpy.polynomial import Polynomial

# Load data for different temperatures
my_data = np.genfromtxt('Data/PFR/25.09.30C.csv', delimiter=';', dtype=None, names=True, encoding='ISO-8859-1')  # Using ISO-8859-1 encoding

def temp_extract(data, x, offset=0):
    rows = data[data['TagName'] == x]
    # Remove invalid rows
    valid_rows = [row for row in rows if row['vValue'] not in ['(null)', None]]
    
    # Parse DateTime for valid rows
    date_times = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_rows]
    vvalues = [float(row['vValue'])+offset for row in valid_rows]
    
    # Calculate elapsed time in minutes
    start_time = date_times[0]
    elapsed_time = [(dt - start_time).total_seconds() / 60 for dt in date_times]

    return elapsed_time, vvalues

# Extract temperature data
elapsed_time_208, temp_208 = temp_extract(my_data, 'T208_PV')
elapsed_time_207, temp_207 = temp_extract(my_data, 'T207_PV')
elapsed_time_206, temp_206 = temp_extract(my_data, 'T206_PV')
elapsed_time_205, temp_205 = temp_extract(my_data, 'T205_PV')
elapsed_time_204, temp_204 = temp_extract(my_data, 'T204_PV')
elapsed_time_203, temp_203 = temp_extract(my_data, 'T203_PV')
elapsed_time_202, temp_202 = temp_extract(my_data, 'T202_PV')
elapsed_time_201, temp_201 = temp_extract(my_data, 'T201_PV')


# Generate evenly spaced time values for interpolation
elapsed_time_interp = np.linspace(min(elapsed_time_208), max(elapsed_time_208), 1000)

# Polynomial fit for the data (degree 5 for a better fit, especially on plateaus)
poly_208 = Polynomial.fit(elapsed_time_208, temp_208, 20)
poly_207 = Polynomial.fit(elapsed_time_207, temp_207, 20)
poly_206 = Polynomial.fit(elapsed_time_206, temp_206, 20)
poly_205 = Polynomial.fit(elapsed_time_205, temp_205, 20)
poly_204 = Polynomial.fit(elapsed_time_204, temp_204, 20)
poly_203 = Polynomial.fit(elapsed_time_203, temp_203, 20)
poly_202 = Polynomial.fit(elapsed_time_202, temp_202, 20)
poly_201 = Polynomial.fit(elapsed_time_201, temp_201, 20)

# Use the polynomial to generate smooth curves
temp_208_smooth = poly_208(elapsed_time_interp)
temp_207_smooth = poly_207(elapsed_time_interp)
temp_206_smooth = poly_206(elapsed_time_interp)
temp_205_smooth = poly_205(elapsed_time_interp)
temp_204_smooth = poly_204(elapsed_time_interp)
temp_203_smooth = poly_203(elapsed_time_interp)
temp_202_smooth = poly_202(elapsed_time_interp)
temp_201_smooth = poly_201(elapsed_time_interp)

# Plot the data
plt.figure(figsize=(12, 8))

# Plot for 208 data
# plt.plot(elapsed_time_208, temp_208, 'o', label='T208 Raw Data', color='#1f77b4')  # Blue
plt.plot(elapsed_time_interp, temp_208_smooth, '-', color='#ff7f0e', label='T208 Smoothed Curve')  # Orange

# Plot for 207 data
# plt.plot(elapsed_time_207, temp_207, 'o', label='T207 Raw Data', color='#2ca02c')  # Green
plt.plot(elapsed_time_interp, temp_207_smooth, '-', color='#d62728', label='T207 Smoothed Curve')  # Red

# Plot for 206 data
# plt.plot(elapsed_time_206, temp_206, 'o', label='T206 Raw Data', color='#9467bd')  # Purple
plt.plot(elapsed_time_interp, temp_206_smooth, '-', color='#8c564b', label='T206 Smoothed Curve')  # Brown

# Plot for 205 data
# plt.plot(elapsed_time_205, temp_205, 'o', label='T205 Raw Data', color='#e377c2')  # Pink
plt.plot(elapsed_time_interp, temp_205_smooth, '-', color='#7f7f7f', label='T205 Smoothed Curve')  # Gray

# Plot for 204 data
# plt.plot(elapsed_time_204, temp_204, 'o', label='T204 Raw Data', color='#bcbd22')  # Yellow-green
plt.plot(elapsed_time_interp, temp_204_smooth, '-', color='#17becf', label='T204 Smoothed Curve')  # Cyan

# Plot for 203 data
# plt.plot(elapsed_time_203, temp_203, 'o', label='T203 Raw Data', color='#ff9896')  # Light Red
plt.plot(elapsed_time_interp, temp_203_smooth, '-', color='#c49c94', label='T203 Smoothed Curve')  # Light Brown

# Plot for 202 data
# plt.plot(elapsed_time_202, temp_202, 'o', label='T202 Raw Data', color='#aec7e8')  # Light Blue
plt.plot(elapsed_time_interp, temp_202_smooth, '-', color='#ffbb78', label='T202 Smoothed Curve')  # Light Orange

# Plot for 201 data
# plt.plot(elapsed_time_201, temp_201, 'o', label='T201 Raw Data', color='#98df8a')  # Light Green
plt.plot(elapsed_time_interp, temp_201_smooth, '-', color='#f7b6d2', label='T201 Smoothed Curve')  # Light Pink

# Labeling the axes
plt.xlabel('Elapsed Time (min)')
plt.ylabel('Temperature (°C)')
plt.title('T200_PV Temperature Data over Time')

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

plt.legend()
plt.grid()
# Show the plot
plt.tight_layout()  # Adjust layout to prevent label overlap
plt.show()
