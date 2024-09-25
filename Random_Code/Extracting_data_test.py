# import numpy as np
# import matplotlib.pyplot as plt
# from datetime import datetime
# my_data = np.genfromtxt('Data/CSTR/Runs 16.09/CSTR 27c.csv', delimiter=';', dtype=None, names=True, encoding=None)

# print(my_data)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV data
df = pd.read_csv("CSTR/CSTR_27c.csv", sep=";")

# Extract columns
tagname = df['TagName'].values
time = df['DateTime'].values
value = df['vValue'].values

# Create a DataFrame with tagname, time, and value
data_1 = {'tagname': tagname, 'time': time, 'value': value}
data_frame = pd.DataFrame(data=data_1)

# Convert time column to datetime for better handling
data_frame['time'] = pd.to_datetime(data_frame['time'])

# Grouping by specific TagNames
group_T400PV = data_frame[data_frame['tagname'] == 'T400_PV']
group_T200PV = data_frame[data_frame['tagname'] == 'T200_PV']
group_T151PV = data_frame[data_frame['tagname'] == 'T151_PV']
group_T150PV = data_frame[data_frame['tagname'] == 'T150_PV']
group_Q210PV = data_frame[data_frame['tagname'] == 'Q210_PV']
group_P120Flow = data_frame[data_frame['tagname'] == 'P120_Flow']
group_P100Flow = data_frame[data_frame['tagname'] == 'P100_Flow']

# Example: Check number of groups (i.e., records per tag)
print(f"T400_PV has {group_T400PV.shape[0]} records")
print(f"T200_PV has {group_T200PV.shape[0]} records")

# Print a sample of the data for each group
print(group_T400PV)
# print(group_T200PV.head())
