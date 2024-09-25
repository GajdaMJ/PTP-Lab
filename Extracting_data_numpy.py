import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# Load data for different temperatures
my_data = np.genfromtxt('Data\\CSTR\\Runs 16.09\\CSTR 27c.csv', delimiter=';', dtype=None, names=True, encoding=None)
data_33c = np.genfromtxt('CSTR\\test 18.09.csv', delimiter=';', dtype=None, names=True, encoding=None)
data_22c = np.genfromtxt('Data\\CSTR\\23.09 22c.csv', delimiter=';', dtype=None, names=True, encoding=None)
data_37c = np.genfromtxt('Data\\CSTR\\25.09 37c.csv', delimiter=';', dtype=None, names=True, encoding=None)

# Function to extract temperature data
def temp_extract(data):
    # Extract rows where TagName is 'T200_PV'
    t200_pv_rows = data[data['TagName'] == 'T200_PV']
    
    # Filter out rows where vValue is null or invalid
    valid_rows = [row for row in t200_pv_rows if row['vValue'] not in ['(null)', None]]
    
    # Parse DateTime for valid rows
    date_times = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_rows]
    
    # Extract vValue and convert it to float
    vvalues = [float(row['vValue']) for row in valid_rows]
    
    # Calculate elapsed time in minutes
    start_time = date_times[0]
    elapsed_time = [(dt - start_time).total_seconds() / 60 for dt in date_times]

    return elapsed_time, vvalues

# Extract temperature data
elapsed_time_27c, temp_27c = temp_extract(my_data)
elapsed_time_33c, temp_33c = temp_extract(data_33c)
elapsed_time_22c, temp_22c = temp_extract(data_22c)
elapsed_time_37c, temp_37c = temp_extract(data_37c)

# Print the temperature range for 27°C and 33°C data
print("Temperature range for 27°C data:", np.max(temp_27c) - np.min(temp_27c))
print("Temperature range for 33°C data:", np.max(temp_33c) - np.min(temp_33c))

# Plot the data
plt.figure(figsize=(12, 8))

# Plot for 27°C data
plt.plot(elapsed_time_27c, temp_27c, label='27°C Data', marker='o', linestyle='-', color='blue')

# Plot for 33°C data
plt.plot(elapsed_time_33c, temp_33c, label='33°C Data', marker='o', linestyle='-', color='orange')

# Plot for 22°C data
plt.plot(elapsed_time_22c, temp_22c, label='22°C Data', marker='o', linestyle='-', color='green')
print(np.max(temp_22c))
# Plot for 37°C data
plt.plot(elapsed_time_37c, temp_37c, label='37°C Data', marker='o', linestyle='-', color='red')

# Labeling the axes
plt.xlabel('Elapsed Time (min)')
plt.ylabel('Temperature (°C)')  # Adjust based on the correct unit for T200_PV
plt.title('T200_PV Temperature Data over Time')

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

# Add a legend
plt.legend()

# Show grid
plt.grid()

# Show the plot
plt.tight_layout()  # Adjust layout to prevent label overlap
plt.show()