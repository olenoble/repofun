import os 
import pandas as pd
import numpy as np
import scipy as sp
from pandas.io.data import DataReader
import datetime as dt
import pylab as plt

os.chdir('/home/olivier/Code/Python/tf_rp')

start_date = dt.datetime.strptime('1995-01-01', '%Y-%m-%d')
end_date = dt.datetime.strptime('2015-11-01', '%Y-%m-%d')

ibm = DataReader('IBM', 'yahoo', start_date, end_date)
aapl = DataReader('AAPL', 'yahoo', start_date, end_date)
msft = DataReader('msft', 'yahoo', start_date, end_date)
intel = DataReader('INTC', 'yahoo', start_date, end_date)
hpq = DataReader('HPQ', 'yahoo', start_date, end_date)

x = ibm['Adj Close']
y = x.tolist()
decay_factor = 63
past_wgt = 1.0

for i in range(1, len(y)):
	y[i] = y[i] + past_wgt * np.exp(-1.0 / decay_factor) * y[i-1]
	wgt = (1 - np.exp(-(i+1.0)/decay_factor)) / (1 - np.exp(-1.0/decay_factor))
	y[i] = y[i] / wgt
	past_wgt = wgt


