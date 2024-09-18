import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
my_data = np.genfromtxt('Data/test_11.09_cstr_AAHFlow.csv', delimiter=';', dtype=None, names=True, encoding=None)



date_times = [datetime.strptime(dt.split('.')[0], '%Y-%m-%d %H:%M:%S') for dt in my_data['DateTime']]

# Extract the 'vValue' column as y-axis data
vvalues = my_data['vValue']

# Plotting the data
plt.figure(figsize=(10, 6))
plt.plot(date_times, vvalues, marker='o', linestyle='-', color='b')

# Labeling the axes
plt.xlabel('DateTime')
plt.ylabel('vValue (ml/min)')
plt.title('vValue over Time')

# Rotate x-axis labels for better readability
plt.xticks(rotation=45)

# Show the plot
plt.tight_layout()  # Adjust layout to prevent label overlap
plt.show()