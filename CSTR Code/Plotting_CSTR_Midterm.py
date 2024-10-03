import numpy as np
import matplotlib.pyplot as plt
from Model_3 import *

# Extracting and modeling data for different temperatures
data_27c = data_extract('Data\\CSTR\\Runs 16.09\\CSTR 27c.csv')
sol_27 = CSTR_model(data_27c[2], data_27c[4], data_27c[3], V=567)

data_33c = data_extract('CSTR\\test 18.09.csv')
sol_33 = CSTR_model(data_33c[2], data_33c[4], data_33c[3], V=567)

data_22c = data_extract('Data\\CSTR\\23.09 22c.csv')
sol_22 = CSTR_model(data_22c[2], data_22c[4], data_22c[3], V=567)

data_37c = data_extract('Data\\CSTR\\25.09 37c.csv')
sol_37 = CSTR_model(data_37c[2], data_37c[4], data_37c[3], V=567)

# Create a 2x2 plot layout
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# Define shared x-limits
x_lim = (0, 35)

# Plotting each solution in the respective subplot
# Top-left: 27°C
axs[0, 0].plot(sol_27[0] / 60, sol_27[1][:, 3] - 273.15, label='Model 27°C')
axs[0, 0].plot(data_27c[0], data_27c[1], label='Real 27°C')
axs[0, 0].set_xlabel('Time (minutes)')
axs[0, 0].set_ylabel('Temperature (°C)')
axs[0, 0].set_title('Temperature at 27°C')
axs[0, 0].legend()
axs[0, 0].grid()
axs[0, 0].set_xlim(x_lim)  # Set shared x-limits

# Top-right: 33°C
axs[0, 1].plot(sol_33[0] / 60, sol_33[1][:, 3] - 273.15, label='Model 33°C')
axs[0, 1].plot(data_33c[0], data_33c[1], label='Real 33°C')
axs[0, 1].set_xlabel('Time (minutes)')
axs[0, 1].set_ylabel('Temperature (°C)')
axs[0, 1].set_title('Temperature at 33°C')
axs[0, 1].legend()
axs[0, 1].grid()
axs[0, 1].set_xlim(x_lim)  # Set shared x-limits

# Bottom-left: 22°C
axs[1, 0].plot(sol_22[0] / 60, sol_22[1][:, 3] - 273.15, label='Model 22°C')
axs[1, 0].plot(data_22c[0], data_22c[1], label='Real 22°C')
axs[1, 0].set_xlabel('Time (minutes)')
axs[1, 0].set_ylabel('Temperature (°C)')
axs[1, 0].set_title('Temperature at 22°C')
axs[1, 0].legend()
axs[1, 0].grid()
axs[1, 0].set_xlim(x_lim)  # Set shared x-limits

# Bottom-right: 37°C
axs[1, 1].plot(sol_37[0] / 60, sol_37[1][:, 3] - 273.15, label='Model 37°C')
axs[1, 1].plot(data_37c[0], data_37c[1] , label='Real 37°C')
axs[1, 1].set_xlabel('Time (minutes)')
axs[1, 1].set_ylabel('Temperature (°C)')
axs[1, 1].set_title('Temperature at 37°C')
axs[1, 1].legend()
axs[1, 1].grid()
axs[1, 1].set_xlim(x_lim)  # Set shared x-limits

# Adjust layout and show plot
plt.tight_layout()
plt.show()
