import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV data
df = pd.read_csv("Data/PFR/PFR.18.09/PFR.25.18.09.csv", sep=";")

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
#temperatures 
T400 = data_frame[data_frame['tagname'] == 'T400_PV'].copy()
T208 = data_frame[data_frame['tagname'] == 'T208_PV'].copy()
T207 = data_frame[data_frame['tagname'] == 'T207_PV'].copy()
T206 = data_frame[data_frame['tagname'] == 'T206_PV'].copy()
T205 = data_frame[data_frame['tagname'] == 'T205_PV'].copy()
T204 = data_frame[data_frame['tagname'] == 'T204_PV'].copy()
T203 = data_frame[data_frame['tagname'] == 'T203_PV'].copy()
T202 = data_frame[data_frame['tagname'] == 'T202_PV'].copy()
T201 = data_frame[data_frame['tagname'] == 'T201_PV'].copy()
T200 = data_frame[data_frame['tagname'] == 'T200_PV'].copy()
#conductivity and flows 
Q210 = data_frame[data_frame['tagname'] == 'QT210_PV'].copy()
P120 = data_frame[data_frame['tagname'] == 'P120_PV'].copy()
P100= data_frame[data_frame['tagname'] == 'P100_PV'].copy()

# Convert time to minutes relative to the first time value for each group

T400['elapsed_minutes'] = (T400['time'] - T400['time'].iloc[0]).dt.total_seconds() / 60
T208['elapsed_minutes'] = (T208['time'] - T208['time'].iloc[0]).dt.total_seconds() / 60
T207['elapsed_minutes'] = (T207['time'] - T207['time'].iloc[0]).dt.total_seconds() / 60
T206['elapsed_minutes'] = (T206['time'] - T206['time'].iloc[0]).dt.total_seconds() / 60
T205['elapsed_minutes'] = (T205['time'] - T205['time'].iloc[0]).dt.total_seconds() / 60
T204['elapsed_minutes'] = (T204['time'] - T204['time'].iloc[0]).dt.total_seconds() / 60
T203['elapsed_minutes'] = (T203['time'] - T203['time'].iloc[0]).dt.total_seconds() / 60
T202['elapsed_minutes'] = (T202['time'] - T202['time'].iloc[0]).dt.total_seconds() / 60
T201['elapsed_minutes'] = (T201['time'] - T201['time'].iloc[0]).dt.total_seconds() / 60
T200['elapsed_minutes'] = (T200['time'] - T200['time'].iloc[0]).dt.total_seconds() / 60
Q210['elapsed_minutes'] = (Q210['time'] - Q210['time'].iloc[0]).dt.total_seconds() / 60
P120['elapsed_minutes'] = (P120['time'] - P120['time'].iloc[0]).dt.total_seconds() / 60
P100['elapsed_minutes'] = (P100['time'] - P100['time'].iloc[0]).dt.total_seconds() / 60

# # Plotting function
# def plotting_variables(var):
#    fig_var, ax = plt.subplots(7,6, figsize=(10, 8), sharex=True)
#    ax = ax.flatten()
#    ax_num = 0
#    for df in var:  # Iterate over the dataframes directly
#        ax[ax_num].plot(df['elapsed_minutes'], df['value'])
#        ax[ax_num].set_xlabel('Time (minutes)')
#        ax[ax_num].set_ylabel('Value')
#        ax[ax_num].set_title(df['tagname'].iloc[0])  # Use tagname for the title
#        ax_num += 1
#    fig_var.suptitle('Variables over Time', fontweight='bold')
#    plt.tight_layout()

# # Plotting the function
# plotting_variables(groups)
# plt.show()

fig_ax, ax = plt.subplots(1,2, figsize = (10,8))
ax = ax.flatten()
ax[0].plot(T208['elapsed_minutes'],T208['value'])
ax[1].plot(Q210['elapsed_minutes'], Q210['value'])
ax[0].set_title('T208')
ax[1].set_title('Q210')
plt.show()






