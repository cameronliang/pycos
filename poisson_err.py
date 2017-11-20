########################################################################################
#
# poisson_err.py 	(c) Cameron J. Liang
#						University of Chicago
#     				    jwliang@oddjob.uchicago.edu
#     				    Cameron.Liang@gmail.com
#
########################################################################################

"""
A Script for calculating Possion confidence level and/or limits. 
Run this code to produce a poission confidence table if it does not exist. 
"""

########################################################################################

import numpy as np
from scipy.optimize import brentq # Find root of function given interval [a,b]
from math import factorial

########################################################################################

def poisson_func(lam,n):
	return ((lam ** n) * np.exp(-lam)) / factorial(n)

########################################################################################

def limits_calc(count, guess,cl, option = ''):
	interval_begin = guess -0.5*guess
	interval_end = guess + 50	
	if count >= 0.:
		if option == 'upper':
			g = lambda lam: np.sum(poisson_func(lam, n) for n in range(0,count+1)) - 0.5*(1-cl)
		elif option == 'lower':
			g = lambda lam: np.sum(poisson_func(lam, n) for n in range(count)) - 0.5*(1+cl)

		x = np.linspace(0,125,10)
		for i in range(len(x)):
			if g(x[i]) > 0:
				interval_begin = x[i]
				break
		for i in range(len(x)):
			if g(x[i]) < 0:
				interval_end = x[i]
				break				
		return brentq(g,interval_begin,interval_end)
	else:
		return 0.

########################################################################################

if __name__ == '__main__':
	
	f = open('./possion_err.dat','w')
	for i in range(136):
		print i 
		if i == 0:
			up = limits_calc(i,i,0.6827,option='upper')
			f.write(('%d\t%f\t%f\n') %(i,up,0))
		else:
			up = limits_calc(i,i,0.6827,option='upper')
			low = limits_calc(i,i,0.6827,option='lower')
			f.write(('%d\t%f\t%f\n') %(i,up,low))

	f.close()

########################################################################################	