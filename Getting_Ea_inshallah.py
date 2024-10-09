import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.optimize import curve_fit
import pandas as pd


# Load data for 27°C
my_data = np.genfromtxt('Data\\CSTR\\Runs 16.09\\CSTR 27c.csv', delimiter=';', dtype=None, names=True, encoding=None)

def temp_extract(data, x="T200_PV", offset=0):
    '''Function to extract data from csv files\n
    data = data path for your csv file. Give as a string \n
    x = Name of the instrument that you want. Default set to T200_PV (CSTR internal temperature) \n
    offset = linear offset for values. Default set to zero \n
    returns elapsed time and values for your
    '''
    # Extract the flow data to determine the starting time
    flow_rows = data[data['TagName'] == "P120_Flow"]
    valid_flow_rows = [row for row in flow_rows if row['vValue'] not in ['(null)', None]]
    flow_values = [float(row['vValue']) for row in valid_flow_rows]
    flow_dates = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_flow_rows]
    start_time = None
    for i in range(1, len(flow_values)):
        if flow_values[i-1] < 1 and flow_values[i] > 1:     # Loop that checks when the AAH pump is turned on and sets that as the start time
            start_time = flow_dates[i]
            break # Stop the loop once flow starts

    temp_rows = data[data['TagName'] == x]  # Only choose the rows for that particular instrument 
    valid_temp_rows = [row for row in temp_rows if row['vValue'] not in ['(null)', None]] # Remove null values
    
    temp_dates = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_temp_rows]
    temp_values = [float(row['vValue']) + offset for row in valid_temp_rows]

    # Calculate elapsed time in minutes from the start_time
    elapsed_time = [(dt - start_time).total_seconds() /60 for dt in temp_dates]

    return elapsed_time, temp_values

# Extract temperature and conductivity
time_temp, temp = temp_extract(my_data)
time_conductivity, conductivity = temp_extract(my_data, x="Q210_PV")

# Conductivity normalization
conductivity_min = np.min(conductivity)
conductivity_max = np.max(conductivity)

window_size = 5
temp = temp[4:]
temp_rolling_avg = pd.Series(temp).rolling(window=window_size, min_periods=1).mean().to_numpy()
x_plot = 1 / (np.array(temp_rolling_avg) + 273.15)
time_conductivity = time_conductivity[4:]
conductivity = conductivity[4:]

X = ((conductivity - conductivity_min) / (conductivity_max - conductivity_min)) ** 2
# Calculate the rolling average with a window size of your choice (e.g., 5)
window_size = 7
X_rolling_avg = pd.Series(X).rolling(window=window_size, min_periods=1).mean().to_numpy()
dXdt_gradient = np.gradient(X_rolling_avg, time_conductivity)

# Define flow rates
V = 473  # ml
elapsed_time_aah, aah_flowrate_vector = temp_extract(my_data, x="P120_Flow")
elapsed_time_water, water_flowrate_vector = temp_extract(my_data, x='P100_Flow')

initial_temperature = np.min(temp)
aah_flowrate = np.median(aah_flowrate_vector)
water_flowrate = np.median(water_flowrate_vector)

total_flowrate = aah_flowrate + water_flowrate
print(total_flowrate)
X = np.clip(X, 1e-10, None)  # Avoid log of negative values or zero
k =  1 / (1-X) * (X * (total_flowrate / V) + dXdt_gradient) 

# Linear fit function
def linear_fit(x, a, b):
    return a * x + b

# Fit the data
params, covariance = curve_fit(linear_fit, x_plot[8:-25], np.log(k[8:-25]))

# Extract parameters
a, b = params
print(f"Fitted parameters: a = {a}, b = {b}")

# Prepare fitted line data for plotting
x_fit = np.linspace(np.min(x_plot[8:-35]), np.max(x_plot[8:-35]), 100)
y_fit = linear_fit(x_fit, a, b)
print(conductivity_max)

print(f"Ea={a*-8.3145}")
print(f"k={np.exp(b)}")

##Plotting
plt.plot(x_plot[8:-35], np.log(k[8:-35]), 'bo', label='Data')
plt.plot(x_fit, y_fit, 'r-', label='Fitted Line: y = {:.2f}x + {:.2f}'.format(a, b))
plt.xlabel('1/(Temperature + 273.15)')
plt.ylabel('log(k)')
plt.title('Linear Fit of log(k) vs 1/(Temperature + 273.15)')
plt.legend()
plt.show()
