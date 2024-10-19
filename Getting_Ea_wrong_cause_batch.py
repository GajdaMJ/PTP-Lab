import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.optimize import curve_fit

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
    valid_temp_rows = [row for row in temp_rows if row['vValue'] not in ['(null)', None]] # You want to remove the values when theres null otherwise it does weird things
    
    temp_dates = [datetime.strptime(row['DateTime'].split('.')[0], '%Y-%m-%d %H:%M:%S') for row in valid_temp_rows] #Converts the weird csv time format to python
    temp_values = [float(row['vValue']) + offset for row in valid_temp_rows]

    # Calculate elapsed time in minutes from the start_time
    elapsed_time = [(dt - start_time).total_seconds()  for dt in temp_dates]

    return elapsed_time, temp_values

time_temp, temp = temp_extract(my_data)
time_conductivity, conductivity = temp_extract(my_data, x="Q210_PV")

conductivity_min = np.min(conductivity)
conductivity_max = np.max(conductivity)

x = 1 / (np.array(temp) + 273.15)
# Clip conductivity_fraction to avoid negative values
conductivity_fraction = ((conductivity_max-conductivity_min)**2 - (conductivity-conductivity_min)**2) / ((conductivity_max-conductivity_min)**2)
conductivity_fraction = np.clip(conductivity_fraction, 1e-10, None)  # Avoid log of negative values or zero

# Clip time_conductivity to avoid division by zero or negative values
time_conductivity = np.clip(time_conductivity, 1e-10, None)

# Logarithmic calculation with clipped values
y = np.log(-1 / np.array(time_conductivity) * np.log(conductivity_fraction))

# Linear function for fitting
def linear_func(x, a, b):
    return a * x + b

# Select the range for x and y (excluding the indices 10 to -45 as in your plot)
x_fit = x[14:-5]
y_fit = y[14:-5]

# Fit the linear model to the data
popt, pcov = curve_fit(linear_func, x_fit, y_fit)

# Extract the parameters a and b
a, b = popt


# Print the results
print(f"Slope (a): {a}")
print(f"Intercept (b): {b}")
print(f"Ea={a*-8.3145}")
print(f"k={np.exp(b)}")
# Plot the original data and the fitted line
plt.plot(x_fit, y_fit, 'bo', label='Data')
plt.plot(x_fit, linear_func(x_fit, a, b), 'r-', label=f'Fit: y={a:.2f}x + {b:.2f}')
plt.xlabel("1 / (Temp + 273.15)")
plt.ylabel("log(-1 / t * log(conductivity_fraction))")
plt.legend()
plt.show()