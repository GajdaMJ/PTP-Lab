import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# Load data for different temperatures
my_data = np.genfromtxt('Data\\CSTR\\Runs 16.09\\CSTR 27c.csv', delimiter=';', dtype=None, names=True, encoding=None)
data_33c = np.genfromtxt('CSTR\\test 18.09.csv', delimiter=';', dtype=None, names=True, encoding=None)
data_22c = np.genfromtxt('Data\\CSTR\\23.09 22c.csv', delimiter=';', dtype=None, names=True, encoding=None)
data_37c = np.genfromtxt('Data\\CSTR\\25.09 37c.csv', delimiter=';', dtype=None, names=True, encoding=None)


def temp_extract(data, x="T200_PV", offset=0):
    # Extract the flow data to determine the starting time
    flow_rows = data[data['TagName'] == "P120_Flow"]
    valid_flow_rows = [row for row in flow_rows if row['vValue'] not in ['(null)', None]]
    flow_values = [float(row['vValue']) for row in valid_flow_rows]
    flow_dates = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_flow_rows]

    # Find the first point where flow transitions from < 1 to > 1
    start_time = None
    for i in range(1, len(flow_values)):
        if flow_values[i-1] < 1 and flow_values[i] > 1:
            start_time = flow_dates[i]
            break

    if start_time is None:
        raise ValueError("No flow transition from < 1 to > 1 found in the data.")

    # Extract temperature data starting from the transition point
    temp_rows = data[data['TagName'] == x]
    valid_temp_rows = [row for row in temp_rows if row['vValue'] not in ['(null)', None]]
    
    temp_dates = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_temp_rows]
    temp_values = [float(row['vValue']) + offset for row in valid_temp_rows]

    # Calculate elapsed time in minutes from the start_time
    elapsed_time = [(dt - start_time).total_seconds() / 60 for dt in temp_dates]

    return elapsed_time, temp_values

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
plt.xlim(0,35)
# Show the plot
plt.tight_layout()  # Adjust layout to prevent label overlap
plt.show()