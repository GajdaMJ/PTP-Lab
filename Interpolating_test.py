import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.optimize import curve_fit

# Load data for 27°C
my_data = np.genfromtxt('Data\CSTR\\25.09 37c.csv', delimiter=';', dtype=None, names=True, encoding=None)


def temp_extract(data, x="T200_PV", offset=0):
    '''Function to extract data from csv files
    data = data path for your csv file. Give as a string 
    x = Name of the instrument that you want. Default set to T200_PV (CSTR internal temperature) 
    offset = linear offset for values. Default set to zero 
    returns elapsed time and values for your
    '''
    # Extract the flow data to determine the starting time
    flow_rows = data[data['TagName'] == "P120_Flow"]
    valid_flow_rows = [row for row in flow_rows if row['vValue'] not in ['(null)', None]]
    
    # Check if we have valid flow rows
    if not valid_flow_rows:
        raise ValueError("No valid flow rows found in the data.")
    
    flow_values = [float(row['vValue']) for row in valid_flow_rows]
    flow_dates = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_flow_rows]
    

    start_time = None
    for i in range(1, len(flow_values)):
        if flow_values[i-1] < 1 and flow_values[i] > 1:     # Loop that checks when the AAH pump is turned on and sets that as the start time
            start_time = flow_dates[i]
            break # Stop the loop once flow starts

    # Handle the case where no start time is found
    if start_time is None:
        start_time = flow_dates[0]  # Use the first flow date if no flow condition is met

    # Extract temperature data for the specified tag
    temp_rows = data[data['TagName'] == x]  # Only choose the rows for that particular instrument 
    valid_temp_rows = [row for row in temp_rows if row['vValue'] not in ['(null)', None]] # Remove null values
    
    temp_dates = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_temp_rows]
    temp_values = [float(row['vValue']) + offset for row in valid_temp_rows]

    # Calculate elapsed time in minutes from the start_time
    elapsed_time = [(dt - start_time).total_seconds() / 60 for dt in temp_dates]

    return elapsed_time, temp_values



time_temp, temp = temp_extract(my_data)
time_conductivity, conductivity = temp_extract(my_data, x="Q210_PV")




plt.plot(time_temp, conductivity)
plt.show()