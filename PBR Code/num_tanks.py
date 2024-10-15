import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime


tracer_data = np.genfromtxt('Data\\PFR\\09.10.TRACER.csv', delimiter=';', dtype=None, names=True, encoding=None)
tracer_data = tracer_data[7:] #adjust for when started

def data_extract(data, x="T200_PV", offset=0):
    '''Function to extract data from csv files
    data = data path for your csv file. Give as a string 
    x = Name of the instrument that you want. Default set to T200_PV (CSTR internal temperature) 
    offset = linear offset for values. Default set to zero 
    returns elapsed time and values for your
    '''

    # Extract temperature data for the specified tag
    temp_rows = data[data['TagName'] == x]  # Only choose the rows for that particular instrument 
    valid_temp_rows = [row for row in temp_rows if row['vValue'] not in ['(null)', None]] # Remove null values
    
    temp_dates = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_temp_rows]
    temp_values = [float(row['vValue']) + offset for row in valid_temp_rows]

    # Calculate elapsed time in minutes from the start_time
    elapsed_time = [(dt - temp_dates[0]).total_seconds() / 60 for dt in temp_dates]

    return elapsed_time, temp_values

tracer_time, conductivity = data_extract(tracer_data, "QT210_PV")
conductivity= conductivity-np.min(conductivity)


plt.plot(tracer_time,conductivity)
plt.show()

tracer_time = np.array(tracer_time)
conductivity = np.array(conductivity)

tau = np.trapz(tracer_time*conductivity,tracer_time)/(np.trapz(conductivity,tracer_time))
sigma = np.trapz((tracer_time-tau)**2 *conductivity,tracer_time)/(np.trapz(conductivity,tracer_time))
print(f"N={tau**2 / sigma}")
e_t = conductivity/(np.trapz(conductivity, tracer_time))
print(np.trapz(e_t, tracer_time))
mrt = np.trapz(tracer_time * e_t, tracer_time)
var = np.trapz(((tracer_time - mrt)**2)*e_t, tracer_time)

print(f'number of tanks={var}')
num_tanks = mrt**2 / var

print(num_tanks)