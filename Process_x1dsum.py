########################################################################################
#
# x1dsum.py  		(c) Cameron J. Liang
#						University of Chicago
#     				    jwliang@oddjob.uchicago.edu
#     				    Cameron.Liang@gmail.com
#
########################################################################################

"""
Script to produce x1dsum ascii files from raw x1dsum.fits
i.e., Produce ascii spectra with corrected low-count possion error

INPUT: 
(1) Full path to QSO directory with x1dum.fits files, qso_directory = '.../.../qso_name_dir'
	-- required all raw data fits file in /qso_name_dir/raw
(2) possion asymmetric error tables in directory = './poisson_err.dat'

OUTPUT:
(a) 'filelist' = a list of all ascii spectra 
(b) 'G130M_filelist' = a list of ascii spectra of G130M grating only 
(c) 'G160M_filelist' = a list of ascii spectra of G160M grating only 
(d) all ascii spectra with format 
	(wave, flux, error, photon_counts, upper_flux_error, lower_flux_error)
"""

########################################################################################

import os
import re
import sys
import operator
import numpy as np
import pyfits as pf
from scipy.interpolate import interp1d

########################################################################################

if len(sys.argv) != 2:
	print 'python Process_x1dsum.py qso_name'
	print 'code exited'
	sys.exit()

filepath = sys.argv[1]
filenames = np.array((os.listdir(filepath + '/raw')))

#output_path = filepath + '/intermediate'
#if not os.path.isdir(output_path):
#	os.mkdir(output_path)

# ------- read in possion error table ------- #
n,dn_plus,dn_minus = (np.loadtxt('./data/poisson_err.dat', 
					  unpack = True, usecols = [0,1,2]))
err_up   = interp1d(n, dn_plus,  kind='linear')
err_down = interp1d(n, dn_minus, kind='linear')


# ----- Write a text file containing list of x1sum filenames, list_f ---- #
list_f 	   = open(filepath + '/filelist','w')
list_G130f = open(filepath + '/G130M_filelist','w')
list_G160f = open(filepath + '/G160M_filelist','w')


# ---- Work on each x1dsum 1D spectrum ------ #
for filename in filenames:
	if re.search('_x1dsum.fits',filename):
		ascii_filename =  filename[:9] + '.x1dsum'
		print ascii_filename		
		
		# Reading info from fits file 
		fitsfile = filepath + '/raw/' + filename
		hdu_list = pf.open(fitsfile)
		wave = hdu_list[1].data['WAVELENGTH'].flatten()
		flux = hdu_list[1].data['FLUX'].flatten()
		error = hdu_list[1].data['ERROR'].flatten()
		net_count_rate = hdu_list[1].data['NET']
		exposure_time = hdu_list[1].data['EXPTIME']
		gcount_rate =  hdu_list[1].data['GROSS']
		
		# Get the net count of the Object 
		array_size_0 = len(net_count_rate[0])
		array_size_1 = len(net_count_rate[1])
		net_count = np.array([np.zeros(array_size_0),np.zeros(array_size_1)])
		net_count[0] = net_count_rate[0] * exposure_time[0]
		net_count[1] = net_count_rate[1] * exposure_time[1] 
		net_count = net_count.flatten()
		gcount = gcount_rate.flatten() * exposure_time[0]

		# Calculate the error from the net counts
		df_minus = np.zeros(len(net_count))
		df_plus  = np.zeros(len(net_count))

		for i in xrange(len(net_count)):			
			netcount = int(net_count[i] + 0.5)
			if net_count[i] > 0 and net_count[i] <= 135:
				if netcount > 0:
					dcount_plus  = err_up(netcount) - netcount
					dcount_minus = netcount - err_down(netcount)
					df_minus[i] = (dcount_minus/ net_count[i]) * flux[i]
					df_plus[i] = (dcount_plus/ net_count[i]) * flux[i]
				
				elif netcount == 0: #if float number < 0.5, it can be rounded to 0. 
					dcount_plus  = err_up(netcount) - netcount
					df_plus[i] = (dcount_plus/ net_count[i]) * flux[i]
					df_minus[i] = 0.

			elif net_count[i] == 0: # if net_count is identically zero == no data
				df_minus[i] = 0.
				df_plus[i] = 0.

			elif net_count[i] < 0:
				gcount_int = int(gcount[i] + 0.5)
				if gcount_int > 0:
					dcount_plus = err_up(gcount_int) - gcount_int
					dcount_minus = gcount_int - err_down(gcount_int)
					df_minus[i] = (dcount_minus/ net_count[i]) * flux[i]
					df_plus[i] = (dcount_plus/ net_count[i]) * flux[i]				
				elif gcount_int == 0: #if float number < 0.5, it can be rounded to 0. 
					dcount_plus = err_up(gcount_int) - gcount_int
					df_plus[i] = (dcount_plus/ net_count[i]) * flux[i]				
					df_minus[i] = 0.

			elif net_count[i] > 135: # If >135, assume gaussian error (at 135, --> 8% frac. error)
				dcount_plus  = np.sqrt(net_count[i])	
				dcount_minus = np.sqrt(net_count[i])	

				df_minus[i] = (dcount_minus/ net_count[i]) * flux[i]
				df_plus[i] = (dcount_plus/ net_count[i]) * flux[i]
			else: 
				df_minus[i] = -1.
				df_plus[i] = -1.

		wave,flux,error,net_count,df_plus,df_minus = [list(x) for x in zip(*sorted(zip(wave,flux,error,net_count,df_plus,df_minus),key=operator.itemgetter(0)))]

		# Produce individual ascii spectrum
		print 'Writing... ', filepath + '/' + ascii_filename
		f = open(filepath + '/' + ascii_filename,'w')
		for i in xrange(len(wave)):
			if flux[i] == 0 and df_plus[i]==[0] and df_minus[i] == 0:
				continue
			else:
				f.write("%f %.32f %.32f %.32f %.32f %.32f\n" % (wave[i],flux[i],error[i],net_count[i],df_plus[i],df_minus[i]))
		f.close()

		# Write a file with a list of ascii spectra files
		list_f.write("%s\n" % ascii_filename)
		if max(wave) > 1500 and min(wave) >= 1300:
			list_G160f.write("%s\n" % ascii_filename)
		else:
			list_G130f.write("%s\n" % ascii_filename)

list_f.close()
list_G130f.close()
list_G160f.close()
print 'Done.'

########################################################################################