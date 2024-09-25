import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("Data/PFR/PFR.16.09.25C/temperatures_1.csv", sep = ";")
time = df['DateTime'].values
aflow = df['vValue'].values

#making the data frame with time as index
data_1 = {'aflow':aflow}
data_frame = pd.DataFrame(data=data_1, index=time)

#the data frame with time and values in columnsw
# data_1 = {'aflow':aflow}
# data_frame = pd.DataFrame(data=data_1, index=time)

x_data = list(np.arange(5, 57.5 + 7.5, 7.5))
y_data = (26.2, 26.7, 26.6, 25.7, 27.2, 27.5, 27.6, 27)

print(data_frame)
# vl8 = aflow[-8:]
# time_l = time[-6:]

plt.plot(x_data, y_data)
plt.show()

# print(vl8) 