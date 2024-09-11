import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt

df_aah = pd.read_csv('Data/test_11.09_cstr_AAHFlow.csv', sep=";")

time_aah = df_aah ['DateTime'].values
aahflow = df_aah ['vValue'].values

#making the data frame with time as index
data_1 = {'aahflow':aahflow}
data_frame = pd.DataFrame(data=data_1, index=time_aah)

#the data frame with time and values in columns
# data_1 = {'aahflow':aahflow}
# data_frame = pd.DataFrame(data=data_1, index=time_aah)

df_temp = pd.read_csv('Data/test_11.09_cstr_Internal_Temp.csv', sep=";")

time_temp = df_temp ['DateTime'].values
temp = df_temp ['vValue'].values

df_w = pd.read_csv('Data/test_11.09_cstr_WaterFlow.csv', sep=";")

time_w = df_w ['DateTime'].values
wflow = df_w ['vValue'].values


df_conductivity1 = pd.read_csv('Data/CSTR_Conductivity_11.09.txt', sep="\t ")

print(df_conductivity1)

#for importing conductivity from a txt file
import pandas as pd

# Define a list to hold the parsed data
data = []

# Open the file and read it line by line
with open('CSTR/CSTR_Conductivity_11.09.txt', 'r') as file:
    lines = file.readlines()

# Process the lines
for i in range(0, len(lines), 2):
    line = lines[i].strip()
    
    # Split the line by tabs to extract relevant data
    parts = line.split('\t')
    
    if len(parts) >= 7:  # Ensure the line contains valid data
        date = parts[0]
        time = parts[1]
        conductivity = parts[2]
        cond_units = parts[3]
        temp = parts[4]
        temp_units = parts[5]
        reftem = parts[6].split('=')[1].strip()  # Extract the reference temperature
        channel = parts[6].split(' ')[-1]  # Extract channel info (CH1)
        
        # Append the data to the list
        data.append([date, time, conductivity, cond_units, temp, temp_units, reftem, channel])

# Define the DataFrame columns
columns = ['date', 'time', 'conductivity', 'cond_units', 'temp', 'temp_units', 'reftem', 'channel']

# Create a DataFrame from the data list
df = pd.DataFrame(data, columns=columns)

# Display the DataFrame
print(df)

