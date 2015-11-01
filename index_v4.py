import bt
import pandas as pd
import numpy as np
import scipy as sp
import os
import datetime as dt

import pylab as plt

os.chdir('/home/olivier/Code/Python/bttest/')
ts_data = pd.read_csv('alltimesseries.csv')

ts_data['Date'] = [dt.datetime.strptime(x,'%Y-%m-%d') for x in ts_data['Date']]
ts_data.set_index('Date', inplace = True)

s = bt.Strategy('s1', [bt.algos.RunMonthly(),
                       bt.algos.SelectAll(),
                       bt.algos.WeighEqually(),
                       bt.algos.Rebalance()])


test = bt.Backtest(s, ts_data['CL1'])
res = bt.run(test)

res.plot()
plt.show()




s2 = bt.Strategy('s2', [bt.algos.RunWeekly(),
                        bt.algos.SelectAll(),
                        bt.algos.WeighInvVol(),
                        bt.algos.Rebalance()])

# now let's test it with the same data set. We will also compare it with our first backtest.
test2 = bt.Backtest(s2, ts_data)
res2 = bt.run(test, test2)

res2.plot()
plt.show()

res2.plot_security_weights()
