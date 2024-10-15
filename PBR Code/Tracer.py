import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.optimize import fsolve
from numpy.polynomial import Polynomial

############################################# MODELING THE ACTUAL DATA ######################################################
# Load data for different temperatures
tracer_data = np.genfromtxt('Data/Calibrations/09.10.TRACER.csv', delimiter=';', dtype=None, names=True, encoding='ISO-8859-1')  # Using ISO-8859-1 encoding

def data_extract(data, x, offset=0):
    rows = data[data['TagName'] == x]
    # Remove invalid rows
    valid_rows = [row for row in rows if row['vValue'] not in ['(null)', None]]
    
    # Parse DateTime for valid rows
    date_times = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_rows]
    vvalues = [float(row['vValue'])+offset for row in valid_rows]
    
    # Calculate elapsed time in minutes
    start_time = date_times[7]
    elapsed_time = [(dt - start_time).total_seconds() / 60 for dt in date_times]

    return elapsed_time, vvalues


time_tracer = data_extract(tracer_data,"QT210_PV")[0]
conductivity = data_extract(tracer_data, "QT210_PV")[1]
max = np.max(conductivity)
min = np.min(conductivity)
max_index = np.argmax(conductivity)
t_max = time_tracer[max_index]
baseline = np.full_like(conductivity, min)

print(time_tracer)

print(t_max,max)

h_max = lambda x: (max-min)/2

plt.plot(time_tracer, conductivity)
plt.plot(time_tracer, baseline, label = "baseline", color = "black")
plt.plot(t_max, max, color = 'red', marker = "o")
plt.xlim(0)
plt.title('Conductivity')
plt.xlabel('Elapsed Time (min)')
plt.ylabel('us/cm')
plt.grid(True)
plt.legend()

plt.show()

