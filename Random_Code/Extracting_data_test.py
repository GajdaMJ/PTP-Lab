# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt

# # Load the CSV data
# df = pd.read_csv("Data/CSTR/Runs 16.09/CSTR 27c.csv", sep=";")

# # Extract columns
# tagname = df['TagName'].values
# time = df['DateTime'].values
# value = df['vValue'].values

# # Create a DataFrame with tagname, time, and value
# data_1 = {'tagname': tagname, 'time': time, 'value': value}
# data_frame = pd.DataFrame(data=data_1)

# # Convert time column to datetime for better handling
# data_frame['time'] = pd.to_datetime(data_frame['time'])

# # Grouping by specific TagNames
# group_T400PV = data_frame[data_frame['tagname'] == 'T400_PV']
# group_T200PV = data_frame[data_frame['tagname'] == 'T200_PV']
# group_T151PV = data_frame[data_frame['tagname'] == 'T151_PV']
# group_T150PV = data_frame[data_frame['tagname'] == 'T150_PV']
# group_Q210PV = data_frame[data_frame['tagname'] == 'Q210_PV']
# group_P120Flow = data_frame[data_frame['tagname'] == 'P120_Flow']
# group_P100Flow = data_frame[data_frame['tagname'] == 'P100_Flow']

# # Example: Check number of groups (i.e., records per tag)
# print(f"T400_PV has {group_T400PV.shape[0]} records")
# print(f"T200_PV has {group_T200PV.shape[0]} records")

# # Print a sample of the data for each group
# print(group_T400PV)
# # print(group_T200PV.head())

# plt.plot(group_T200PV['time'], group_T200PV['value'])
# plt.show()


######################################################################




import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV data
df = pd.read_csv("Data/CSTR/Runs 16.09/CSTR 27c.csv", sep=";")

# Extract columns
tagname = df['TagName'].values
time = df['DateTime'].values
value = df['vValue'].values

# Create a DataFrame with tagname, time, and value
data_1 = {'tagname': tagname, 'time': time, 'value': value}
data_frame = pd.DataFrame(data=data_1)

# Convert time column to datetime for better handling
data_frame['time'] = pd.to_datetime(data_frame['time'])

# Grouping by specific TagNames (e.g., T200_PV)
group_T400PV = data_frame[data_frame['tagname'] == 'T400_PV'].copy()
group_T200PV = data_frame[data_frame['tagname'] == 'T200_PV'].copy()
group_T151PV = data_frame[data_frame['tagname'] == 'T151_PV'].copy()
group_T150PV = data_frame[data_frame['tagname'] == 'T150_PV'].copy()
group_Q210PV = data_frame[data_frame['tagname'] == 'Q210_PV'].copy()
group_P120Flow = data_frame[data_frame['tagname'] == 'P120_Flow'].copy()
group_P100Flow = data_frame[data_frame['tagname'] == 'P100_Flow'].copy()

# Convert time to minutes relative to the first time value for each group
group_T400PV['elapsed_minutes'] = (group_T400PV['time'] - group_T400PV['time'].iloc[0]).dt.total_seconds() / 60
group_T200PV['elapsed_minutes'] = (group_T200PV['time'] - group_T200PV['time'].iloc[0]).dt.total_seconds() / 60
group_T151PV['elapsed_minutes'] = (group_T151PV['time'] - group_T151PV['time'].iloc[0]).dt.total_seconds() / 60
group_T150PV['elapsed_minutes'] = (group_T150PV['time'] - group_T150PV['time'].iloc[0]).dt.total_seconds() / 60
group_Q210PV['elapsed_minutes'] = (group_Q210PV['time'] - group_Q210PV['time'].iloc[0]).dt.total_seconds() / 60
group_P120Flow['elapsed_minutes'] = (group_P120Flow['time'] - group_P120Flow['time'].iloc[0]).dt.total_seconds() / 60
group_P100Flow['elapsed_minutes'] = (group_P100Flow['time'] - group_P100Flow['time'].iloc[0]).dt.total_seconds() / 60

# List of the grouped dataframes to plot
variables = [group_T200PV, group_Q210PV, group_P120Flow, group_P100Flow]

# Plotting function
def plotting_variables(var):
   fig_var, ax = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
   ax = ax.flatten()
   ax_num = 0
   for df in var:  # Iterate over the dataframes directly
       ax[ax_num].plot(df['elapsed_minutes'], df['value'])
       ax[ax_num].set_xlabel('Time (minutes)')
       ax[ax_num].set_ylabel('Value')
       ax[ax_num].set_title(df['tagname'].iloc[0])  # Use tagname for the title
       ax_num += 1
   fig_var.suptitle('Variables over Time', fontweight='bold')
   plt.tight_layout()

# Plotting the function
plotting_variables(variables)
plt.show()

