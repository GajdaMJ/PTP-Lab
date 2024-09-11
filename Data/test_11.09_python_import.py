import numpy as np 
import csv 

file = open('test_11.09_AAHFlow.csv')
type(file)

csvreader = csv.reader(file)

header = []
header = next(csvreader)
header