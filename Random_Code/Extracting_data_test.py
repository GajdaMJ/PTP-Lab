import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
my_data = np.genfromtxt('Data\CSTR\Runs 18.09\cstr 33c.csv', delimiter=';', dtype=None, names=True, encoding=None)

print(my_data)



import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("CSTR/CSTR_27c.csv", sep = ";")
tagname = df['TagName'].values
time = df['DateTime'].values
value = df['vValue'].values


#making the data frame with time as index
data_1 = {'tagname' :tagname, 'time': time,'value':value}
data_frame = pd.DataFrame(data=data_1, index=tagname)

#the data frame with time and values in columns
# data_1 = {'aflow':aflow}
# data_frame = pd.DataFrame(data=data_1, index=time)

print(data_frame)

