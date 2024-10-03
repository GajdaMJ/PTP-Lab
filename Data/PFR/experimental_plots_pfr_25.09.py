import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime


data_22c=  np.genfromtxt('25.09.22C(att.55.conductivityweird).csv', delimiter=';', dtype=None, names=True, encoding=None)
data_30c= np.genfromtxt('25.09.30C.csv', delimiter=';', dtype=None, names=True, encoding=None)
data_33c= np.genfromtxt('25.09.33C.csv', delimiter=';', dtype=None, names=True, encoding=None)

def temp_extract(data,x, offset=0):
    # Extract rows where TagName is 'T200_PV'
    t200_pv_rows = data[data['TagName'] == x]
     
    # Filter out rows where vValue is null or invalid
    valid_rows = [row for row in t200_pv_rows if row['vValue'] not in ['(null)', None]]
    
    # Parse DateTime for valid rows
    date_times = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_rows]

    # Extract vValue and convert it to float, calibration added for the thermocouples
    vvalues = [float(row['vValue']) + offset for row in valid_rows]

    
    # Calculate elapsed time in minutes
    start_time = date_times[0]
    elapsed_time = [(dt - start_time).total_seconds() / 60 for dt in date_times]

    return elapsed_time, vvalues

t_values = ['T208_PV','T207_PV','T206_PV','T205_PV','T204_PV','T203_PV','T202_PV','T201_PV','T200_PV']
t_label = ["t9", "t8", 't7', 't6', 't5', 't4', 't3', 't2', 't1']

#########################
# Plot the data

t_22_0 = [21.7,21.7,22.7,22.7,21.7,22.6,23,22.8,23]
plt.figure(figsize=(6, 8))
for i in range (0,8):
    plt.plot(temp_extract(data_22c,t_values[i])[0], temp_extract(data_22c,t_values[i], offset=-t_22_0[i]+22)[1], label = t_label[i])

# Labeling the axes
plt.xlabel('Elapsed Time (min)')
plt.ylabel('Temperature (°C)')  # Adjust based on the correct unit for T200_PV
plt.title('Temperatures at different locations over Time (t0 = 22)')
plt.legend()

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

# Show the plot
plt.tight_layout()  # Adjust layout to prevent label overlap
plt.show()

##################################
# Plot the data 30C

t_30_0 = [29.1,29.1,30.2,30.1,29.2,30.1,30.3,30.4,30.4]
plt.figure(figsize=(6, 8))
for i in range (0,8):
    plt.plot(temp_extract(data_30c,t_values[i])[0], temp_extract(data_30c,t_values[i], offset=-t_30_0[i]+30)[1], label = t_label[i])

# Labeling the axes
plt.xlabel('Elapsed Time (min)')
plt.ylabel('Temperature (°C)')  # Adjust based on the correct unit for T200_PV
plt.title('Temperatures at different locations over Time (t0 = 30)')
plt.legend()

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

# Show the plot
plt.tight_layout()  # Adjust layout to prevent label overlap
plt.show()

###########################################

# Plot the data 33C
t_33_0 = [31.7,31.7,32.9,32.9,31.7,32.9,33.4,33.4,33.2]
plt.figure(figsize=(6, 8))
for i in range (0,8):
    plt.plot(temp_extract(data_33c,t_values[i])[0], temp_extract(data_33c,t_values[i], offset=-t_33_0[i]+33)[1], label = t_label[i])

# Labeling the axes
plt.xlabel('Elapsed Time (min)')
plt.ylabel('Temperature (°C)')  # Adjust based on the correct unit for T200_PV
plt.title('Temperatures at different locations over Time (t0 = 33)')
plt.legend()

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

# Show the plot
plt.tight_layout()  # Adjust layout to prevent label overlap
plt.show()



##########################

# Plot COMPARISSON
# plt.figure(figsize=(6, 8))
# plt.plot(temp_extract(data_22c,'T208_PV')[0], temp_extract(data_22c,'T208_PV', offset=-21.7+22)[1], label = "t0 = 22")
# plt.plot(temp_extract(data_30c,'T208_PV')[0], temp_extract(data_30c,'T208_PV', offset=-29.1+30)[1], label = "t0 = 30")
# plt.plot(temp_extract(data_33c,'T208_PV')[0], temp_extract(data_33c,'T208_PV', offset=-31.7+33)[1], label = "t0 = 33")


# # Labeling the axes
# plt.xlabel('Elapsed Time (min)')
# plt.ylabel('vValue (°C)')  # Adjust based on the correct unit for T200_PV
# plt.title('Last Reactor Temperatures of Different Initial Temperatures Over Times')
# plt.legend()

# # Rotate x-axis labels for better readability
# plt.xticks(rotation=45)

# # Show the plot
# plt.tight_layout()  # Adjust layout to prevent label overlap
# plt.show()

#############################
# Plot comparison diff

plt.figure(figsize=(6, 8))
plt.plot(temp_extract(data_22c,'T208_PV')[0], temp_extract(data_22c,'T208_PV', offset=-21.7)[1], label = "t0 = 22")
plt.plot(temp_extract(data_30c,'T208_PV')[0], temp_extract(data_30c,'T208_PV', offset=-29.1)[1], label = "t0 = 30")
plt.plot(temp_extract(data_33c,'T208_PV')[0], temp_extract(data_33c,'T208_PV', offset=-31.7)[1], label = "t0 = 33")


# Labeling the axes
plt.xlabel('Elapsed Time (min)')
plt.ylabel('Change in Temperature (°C)')  # Adjust based on the correct unit for T200_PV
plt.title('Change in Last Reactor Temperatures of Different Initial Temperatures Over Time')
plt.legend()

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

# Show the plot
plt.tight_layout()  # Adjust layout to prevent label overlap
plt.show()