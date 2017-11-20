"""
Wavelength Dependent Fit. 
"""
import numpy as np
import matplotlib.pyplot as pl
import os,re,sys

from utilities import Read_FitsPoints, Make_AllFitsPoints_Plots

	
if __name__ == '__main__':
	input_path = sys.argv[1]
	Make_AllFitsPoints_Plots(input_path)