import os 
import pandas as pd
import numpy as np
import scipy as sp
from scipy.optimize import minimize
from pandas.io.data import DataReader
import datetime as dt
import pylab as plt

os.chdir('/home/zoe/Code/Python/rp_test')

def rp_fun(w):
	return -1.0 * np.sum(np.log(w))

### inputs
vol_ref = 0.25
rho_ref = 0.95
vol_vol = 0.5
step = 1.0 / 252

l_ts = 252 * 10
n_undl = 15

volt = 0.05
trend = 0.00


### set up basic parameters
cor_mat = np.array([[rho_ref] * n_undl] * n_undl) + np.diag([1 - rho_ref] * n_undl)
cor_chol = np.linalg.cholesky(cor_mat)
roll_i = range(21, l_ts - 21, 21)
idx_lvl = 100
idxhist = pd.DataFrame()

### run strat
range_i = range(0, len(roll_i)-1) 
for i in range_i:

	print "Roll " + str(i) + "/" + str(len(range_i))
	rolltemp = roll_i[i]

	voltemp = np.diag(vol_ref * np.exp(-0.5 * vol_vol ** 2  + vol_vol * np.random.normal(size = n_undl)))
	covmat = np.dot(voltemp, np.dot(cor_mat, voltemp))
	cons = {'type': 'ineq', 'fun': lambda x, c, v:  -np.dot(x, np.dot(c,x)) + v, 'args': (covmat, volt * volt,)}
	res = minimize(rp_fun, [0.001] * n_undl, bounds = [(0, 10000)] * n_undl, method='SLSQP', constraints=cons, options = {'maxiter': 500})

	if not res['success']:
		print res['message']
		
	wgt = res['x']

	dw = np.random.normal(size = 21 * n_undl); dw = (dw - np.mean(dw)) / np.std(dw)
	dw = dw.reshape(21, n_undl)
	dw = np.dot(cor_chol,dw.T).T
	
	voltemp_post = vol_ref * np.exp(-0.5 * vol_vol ** 2  + vol_vol * np.random.normal(size = n_undl))
	
	tstemp = pd.DataFrame(-0.5 * (voltemp_post ** 2) * step + trend * step + voltemp_post * np.sqrt(step) * dw)
	tstemp = np.exp(tstemp.apply(np.cumsum, 0)) 

	#index_ret = s1_idx.apply(lambda x: x / x[min(x.index)] - 1, 0) * wgt
	#index_new = idx_lvl * (1 + np.cumsum(index_ret.apply(sum, 1)))

	index_new = idx_lvl * np.cumprod(1 + (tstemp.apply(lambda x: x / x.shift() - 1, 0) * wgt).fillna(0).apply(sum, 1))
	idx_lvl = index_new.values[-1]
	idxhist = idxhist.append(pd.DataFrame(index_new[:-1], columns = ['test']), ignore_index=True)

print ""
realvol = (np.std(idxhist.apply(lambda x: x / x.shift() - 1))*np.sqrt(1/step)).values[0]
print "Realized Vol: " + str(np.round(realvol * 10000) / 100) + "%"
print "Adjustment (* volvol): " + str(np.round(np.log(realvol / volt) / (vol_vol ** 2) * 100) / 100)


idxhist.plot(); plt.show()
