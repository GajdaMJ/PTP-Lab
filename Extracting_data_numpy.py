import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime


my_data = np.genfromtxt('Data\CSTR\Runs 16.09\CSTR 27c.csv', delimiter=';', dtype=None, names=True, encoding=None)
data_33c =  np.genfromtxt('CSTR\\test 18.09.csv', delimiter=';', dtype=None, names=True, encoding=None)
data_22c= np.genfromtxt('Data\CSTR\\23.09 22c.csv', delimiter=';', dtype=None, names=True, encoding=None)
data_37c= np.genfromtxt('Data\CSTR\\25.09 37c.csv', delimiter=';', dtype=None, names=True, encoding=None)

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



# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(temp_extract(my_data)[0], temp_extract(my_data)[1], marker='o', linestyle='-', color='b')
plt.plot(temp_extract(data_22c)[0], temp_extract(data_22c)[1])
plt.plot(temp_extract(data_37c)[0], temp_extract(data_37c)[1])
plt.plot(temp_extract(data_33c)[0], temp_extract(data_33c)[1])

# Labeling the axes
plt.xlabel('Elapsed Time (min)')
plt.ylabel('vValue (°C)')  # Adjust based on the correct unit for T200_PV
plt.title('T200_PV vValue over Time')

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)
plt.legend()
# Show the plot
plt.tight_layout()  # Adjust layout to prevent label overlap
plt.show()