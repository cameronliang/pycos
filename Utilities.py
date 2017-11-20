########################################################################################
#
# Utilities.py  	(c) Cameron J. Liang
#						University of Chicago
#     				    jwliang@oddjob.uchicago.edu
#     				    Cameron.Liang@gmail.com
#
########################################################################################

"""
Utitilies used throughout PyCos
"""

import sys
import os
import warnings
import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

# Speed of light
c = 299792.458

########################################################################################
#
# Plotting/Align Utitilies
#
########################################################################################


def zoom_region(central_wavelength,dwave,dflux_window_up,dflux_window_down):
	''' Define zoom in window in spectrum '''
	wave_window_left = central_wavelength - dwave
	wave_window_right = central_wavelength + dwave
	flux_window_up = -1e-14 + dflux_window_up
	flux_window_down = 2e-14 + dflux_window_down
	return [wave_window_left,wave_window_right],[flux_window_up, flux_window_down]

########################################################################################

def Select_Data(spec,wave_edges):
	"""
	Select the wavelength range to give a straight line 
	or gaussian fitting, or apparent optical depth. 
	or equivalent width. 
	"""
	if len(wave_edges) < 2: 
		raise ValueError('must be at least two bin edges!')

	wave,flux,error,dfp,dfm= spec
	inds = np.where( (wave > wave_edges[0]) & (wave <= wave_edges[1]) )[0]

	return wave[inds], flux[inds],error[inds], dfp[inds], dfm[inds]

########################################################################################

def Print_LineInfo(gauss_params,logN,transition_name):
	amplitude,centroid_wave,sigma_width = gauss_params
	transition_dict = ReadTransitionData()
	rest_wave = transition_dict[transition_name].wave
	fitted_redshift = round((centroid_wave/rest_wave) - 1., 5)
	sigma_v = (sigma_width/rest_wave)*c

	# print function for Gaussian fit.
	print '\n'
	print 'Gaussian Fit(amp,centroid,sigma)[A] = ', amplitude,centroid_wave,sigma_width
	print 'Redshift = ', fitted_redshift
	print 'Sigma[km/s] = ', sigma_v
	print 'FWHM[km/s] = ', sigma_v*2*np.sqrt(2*np.log(2))
	print 'logN = ', logN
	print '\n'	

########################################################################################

def PrintIonsNames():
	printLine()
	print 'High ions:'
	print 'NeVIIIa 780.324     |', 'NVb     1242.804    |'
	print 'NeVIIIb 770.409     |', 'SiIVa   1393.76018  |'
	print 'OVIa    1031.9261   |', 'SiIVb   1402.77291  |'
	print 'OVIb    1037.6167   |', 'CIVa    1548.2049   |'
	print 'NVa     1238.821    |', 'CIVb    1550.77845  |'
	printLine()
	print 'Low ions:'
	print 'CIII    977.0201    |', 'SIIc    1259.519    |'
	print 'FeIIa   1144.9379   |', 'SiIIc   1260.4221   |'
	print 'SiIIa   1190.4158   |', 'OI      1302.1685   |'
	print 'SiIIb   1193.2897   |', 'SiIId   1304.3702   |'
	print 'NIa     1199.5496   |', 'CII     1334.5323   |'
	print 'NIb     1200.2233   |', 'CII*    1335.7077   |' 
	print 'NIc     1200.7098   |', 'SiIIe   1526.70698  |'
	print 'SiIII   1206.5000   |', 'FeIIb   1608.45085  |'
	print 'SIIa    1250.584    |', 'AlII    1670.7886   |'
	print 'SIIb    1253.811    |', 'NiII	1703.4119   |'
	printLine()
	print 'Warning: be careful for the choice of MW lines for the use of obsolute calibration.'

########################################################################################

def UserInputs():
	transition_dict = ReadTransitionData()
	transition_name = raw_input('Transition Name: ')
	redshift 		= 0.0
	line_region = (1+redshift) * transition_dict[transition_name].wave
	return transition_name,redshift,line_region

########################################################################################

def write_shiftspec(spec, params, spec_path):
	output_fname = raw_input('Output File name:\n')
	output_fname = spec_path + '/' +output_fname

	wave,flux,error,dfp,dfm = spec
	
	pivot_wave = np.median(wave)
	delta_wave = np.polyval(params,wave-pivot_wave)
	
	rectified_wave = wave + delta_wave

	final_spec_file = open(output_fname,'w')
	final_spec_file.write('# Wave\tflux\tmean_df\tdf_up\tdf_down\n')
	for i in xrange(len(rectified_wave)):
		(final_spec_file.write('%f\t%.32f\t%.32f\t%.32f\t%.32f\n' % (
			rectified_wave[i],
			flux[i],
			error[i],
			dfp[i],
			dfm[i])))
	final_spec_file.close()
	print 'Written %s\n' % output_fname
	return None

########################################################################################

def plot_dwave_fit(rest_wave,dlambda,pivot_wave):
	fit_order = determine_fit_order(dlambda)
	fit_params = np.polyfit(rest_wave-pivot_wave,dlambda,fit_order)
	fit_dwave = np.polyval(fit_params,rest_wave-pivot_wave)

	# Visually show the recorded data and best fit
	pl.figure(2,figsize=(10,6))
	pl.plot(rest_wave,dlambda,'mo',ms = 10,label=r'$\rm Recorded shifts$')
	pl.plot(rest_wave,fit_dwave,lw=2,label=r'$\rm Best Fit$')
	pl.legend(loc='best')
	pl.xlabel(r'$\rm Wavelength \lbrack \AA \rbrack$')
	pl.ylabel(r'$\Delta \lambda \rm \lbrack \AA \rbrack$')
	pl.show()	

def print_dwave_fit(transitions,transtion_rest_wave,obs_centroid,dlambda):
	print 'Transition\tRest wave\tObs. Wave \tdwave'
	for i in xrange(len(dlambda)):
		print '%s\t%f\t%f\t%f' % (transitions[i], transtion_rest_wave[i],obs_centroid[i],dlambda[i])
	return 

########################################################################################

def determine_fit_order(dlambda):
	if len(dlambda) == 1:
		print 'Only one data point available - constant offset is used.'
		return 0
	else:
		print 'Choosing fitting orders...'
		print '0 = constant offset; 1 = linear ...'
		fit_order = int(raw_input("Fitting order: "))
		return fit_order


########################################################################################
#
# Computations
#
########################################################################################

def Gaussian_function(x, norm, mean, sigma):
	return norm * np.exp(-(x-mean)**2./(2.* sigma**2.))

def straight_line(p1,p2,spec):
	"""
	Given two points make a straight line. 
	Output: an interpolated linear function 
	"""
	x1,y1 = p1; x2,y2 = p2
	
	m = (y2-y1) / (x2-x1)
	x = np.arange(x1-0.01, x2+0.01, 0.01) 

	y = m *(x-x1) + y1
	f = np.vectorize(interp1d(x,y,kind = 'linear'))
	wave,flux,error = Select_Data(spec,[x1,x2])
	return wave,f(wave)

########################################################################################

def Fit_Gaussian(spec,estimated_cont_level):
	wave,flux,error,dfp,dfm= spec

	try:
		popt,pcov = curve_fit(Gaussian_function, wave, flux-estimated_cont_level, 
							p0 = np.array([-estimated_cont_level,np.mean(wave),0.2]), sigma = error)
	except RuntimeError:
		print("Error - Gaussian Fit Failed")
		return None

	amplitude,centroid_wave,sigma_width = popt
	return amplitude,centroid_wave,sigma_width

########################################################################################

def ComputeAppColumn(spec,transition_name):
	"""
	Estimated of Apparent column density
	"""
	c  = 2.99792458e10    # speed of light
	me = 9.10938291e-28   # electron mass
	e  = 4.80320425e-10   # electron charge
	mH = 1.67e-24         # proton mass
	ang_cm = 1.e-8 		  # Convert Angstrom to cm

	transition_dict = ReadTransitionData()
	rest_wave = transition_dict[transition_name].wave * ang_cm
	osc_f = transition_dict[transition_name].osc_f

	wave,flux,error,dfp,dfm = spec

	dwave = (wave[1] - wave[0])
	temp_sum = np.sum(np.log(1/flux))*dwave * ang_cm
	N = (me*c**2 / (np.pi*osc_f*(rest_wave*e)**2)) * temp_sum

	return np.log10(N)

def chi_squared_calc(wave1, flux1, sig1, wave2, flux2, sig2, wa,wb):
	wave1 = np.array(wave1); flux1 = np.array(flux1)
	wave2 = np.array(wave2); flux2 = np.array(flux2)
	
	start_ind1 = np.where(abs(wave1 - wa) < 0.1)[0][0]
	end_ind1 = np.where(abs(wave1 - wb) < 0.1)[0][0]

	start_ind2 = np.where(abs(wave2 - wa) < 0.1)[0][0]

	num_of_shift = 100
	dpix = np.zeros(num_of_shift); chi2 = np.zeros(num_of_shift)
	
	tempflux1 = flux1[start_ind1:end_ind1]
	tempwave1 = wave1[start_ind1:end_ind1]
	tempsig1 = sig1[start_ind1:end_ind1]

	delta = end_ind1 - start_ind1

	for i in xrange(-num_of_shift/2,num_of_shift/2):
		dpix[i] =  -1*i	
		start  = start_ind2 + i
		end = start_ind2 + delta + i

		tempwave2 = wave2[start:end]
		tempflux2 = flux2[start:end]
		tempsig2 = sig2[start:end]

		#Errors are correlated?
		sigma2 = tempsig1**2 + tempsig2**2
		chi2[i] = np.sum(((tempflux1- tempflux2)**2) / sigma2) / (end_ind1- start_ind1)
		#chi2[i] = np.log10(np.sum(((tempflux1- tempflux2)**2)) / (end_ind1- start_ind1))
	dpix, chi2 = zip(*sorted(zip(dpix, chi2)))
	return dpix, chi2

########################################################################################

def break_wave_calc(wave):
	for i in xrange(len(wave)):
		if wave[i+1] - wave[i] > 5.0:
			break_wavelength = wave[i] + 1.0
			break
	return break_wavelength

########################################################################################


########################################################################################
#
# Testing/Other Utilities
#
########################################################################################


def rounded_obs_wave(obs_wave,rest_wave):
	"""
	Return rest wavelength if obs_wave is within one resolution element
	"""
	dv = 15 # km/s
	obs_wave_left = (dv*c/rest_wave) - rest_wave
	obs_wave_right = (dv*c/rest_wave) + rest_wave

	if obs_wave_left <= obs_wave < obs_wave_right:
		return rest_wave
	else:
		return obs_wave

########################################################################################

def Read_FitsPoints(input_path,fname):
	"""
	Read a list of points of (wave, dwave) from 
	output of align.py
	"""
	wave = []; dwave = []
	break_ind = 0 	

	filename_fullpath = input_path + fname
	with open(filename_fullpath) as f:
		lines = f.readlines()
		for i in xrange(1,len(lines)):
			if re.search('a',lines[i]):
				break_ind = i -1

	wave = np.loadtxt(filename_fullpath,usecols=[0])
	dwave = np.loadtxt(filename_fullpath,usecols=[1])
	
	# Section b
	wave_b = wave[:break_ind]; dwave_b = dwave[:break_ind]
	
	# Section a
	wave_a = wave[break_ind:]; dwave_a = dwave[break_ind:]
	
	pl.plot(wave_a,dwave_a,'o',label='Segment a')
	pl.plot(wave_b,dwave_b,'o',label='Segment b')
	pl.legend(loc='best')
	pl.ylim([-0.1,0.1])
	pl.xlabel(r'Wavelength $\AA$')
	pl.ylabel(r'$\Delta \lambda \AA$')
	pl.savefig(input_path + 'plots/' + fname + '.png')

	pl.clf()

	return [wave_a,dwave_a], [wave_b,dwave_b]

########################################################################################

def Make_AllFitsPoints_Plots(input_path):

	if not os.path.isdir(input_path + '/plots'):
		os.makedirs(input_path + '/plots')

	filelists = os.listdir(input_path)
	for i in xrange(len(filelists)):
		if re.search('fits_',filelists[i]):
			Wave_dwave(input_path, filelists[i])		
	return 

########################################################################################

def Create_Final_NativeWave(spec_cube,combine_grating):
	"""
	Create a wavelength array that uses the original 
	wavelength arrays as the final resolution.
	It uses the wavelength from smaller values first. 
	then switch to another wavelength array when data 
	runs out. It also fill the gaps using a constant dlambda 
	where dlambda = median(resolution)

	Parameters
	---------------------------------------------------------------------------
	spec_cube: array
		Multi-dimensional array of spectral data 

	Returns
	---------------------------------------------------------------------------	
	final_wave_bucket: 1D array
		Final wavelength array based on the native resolution. 

	see also
	---------------------------------------------------------------------------
	Coadd_interp and Coadd_func
	"""
	Num_file = len(spec_cube)

	# Order the input wavelength arrays from small to large
	starting_wave = np.zeros(Num_file)
	inds_orders = np.arange(0,Num_file)
	for n in xrange(Num_file):
		starting_wave[n] = spec_cube[n][0][0]
	new_starting_wave, new_inds_orders = zip(*sorted(zip(starting_wave,inds_orders)))

	# Create new arrays of wavelengths. 
	wave_bucket = [spec_cube[int(new_inds_orders[0])][0]]
	for n in xrange(Num_file):
		wave_bucket.append(spec_cube[int(new_inds_orders[n])][0])

	# Stitch together the wavelength into one
	final_wave_bucket = [wave_bucket[0]]
	for n in xrange(1,Num_file+1):

		# Find the max wavelength from the last spec
		max_wave = np.max(wave_bucket[n-1]) 

		# Find all the inds in the current that is bigger than the last. 
		inds = np.where(wave_bucket[n] > max_wave)[0] 

		# Append those pixels from the last. 
		final_wave_bucket.append(wave_bucket[n][inds])

	final_wave_bucket = np.array(np.hstack(final_wave_bucket))

	# If combine G130M and G160M grating only:
	# This is required because the two graitings give a bimodal distribution 
	# of resolution. The median res is ill-defined. 
	if combine_grating == True:
		return final_wave_bucket
	else:
		dwave = final_wave_bucket[1:] - final_wave_bucket[:-1]
		inds = np.where(dwave>np.median(dwave))[0]
		if len(inds)>0:
			# Fill all the data gaps (if any)
			for i in xrange(len(inds)):
				wave_fillgap = np.arange(final_wave_bucket[inds[i]]+np.median(dwave), 
										 final_wave_bucket[inds[i]+1],np.median(dwave))
				final_wave_bucket = np.append(final_wave_bucket,wave_fillgap)
			return sorted(np.hstack(final_wave_bucket))
		else:
			print 'No filling between wavelength gaps...'
			return final_wave_bucket
	
########################################################################################

def Create_Final_NativeWave_test(spec_cube,combine_grating):
	"""
	A test to speed up and fixes the Create_Final_NativeWave function

	Parameters
	---------------------------------------------------------------------------
	spec_cube: array
		Multi-dimensional array of spectral data 

	Returns
	---------------------------------------------------------------------------	
	final_wave_bucket: 1D array
		Final wavelength array based on the native resolution. 

	see also
	---------------------------------------------------------------------------
	Coadd_interp and Coadd_func
	"""
	Num_file = len(spec_cube)

	# Order the input wavelength arrays from small to large
	starting_wave = np.zeros(Num_file)
	inds_orders = np.arange(0,Num_file)
	for n in xrange(Num_file):
		starting_wave[n] = spec_cube[n][0][0]
	new_starting_wave, new_inds_orders = zip(*sorted(zip(starting_wave,inds_orders)))

	# Create new arrays of wavelengths. 
	wave_bucket = [spec_cube[int(new_inds_orders[0])][0]]
	for n in xrange(Num_file):
		wave_bucket.append(spec_cube[int(new_inds_orders[n])][0])

	# Stitch together the wavelength into one
	final_wave_bucket = [wave_bucket[0]]
	for n in xrange(1,Num_file+1):

		# Find the max wavelength from the last spec
		max_wave = np.max(wave_bucket[n-1]) 

		# Find all the inds in the current that is bigger than the last. 
		inds = np.where(wave_bucket[n] > max_wave)[0] 

		# Append those pixels from the last. 
		final_wave_bucket.append(wave_bucket[n][inds])

	final_wave_bucket = np.array(np.hstack(final_wave_bucket))

	#return final_wave_bucket

	# Now fill all the gaps in the wavelength array 
	dwave = final_wave_bucket[1:] - final_wave_bucket[:-1]
	inds = np.where(dwave>np.median(dwave))[0]
	if len(inds)>0:
		# Fill all the data gaps (if any)
		for i in xrange(len(inds)):
			wave_fillgap = np.arange(final_wave_bucket[inds[i]]+np.median(dwave), 
									 final_wave_bucket[inds[i]+1],np.median(dwave))
			final_wave_bucket = np.append(final_wave_bucket,wave_fillgap)
		return sorted(np.hstack(final_wave_bucket))
	else:
		print 'No filling between wavelength gaps...'
		return final_wave_bucket


########################################################################################

def Create_Constant_WavelengthArray(spec_cube,final_wave_start,final_wave_end):
	""" Does not work  -- patten noise is large """
	dwave = np.zeros(len(spec_cube))
	for n in xrange(len(spec_cube)):
		temp_final_wave  = spec_cube[n][0] # Take one of the spectrum use its resolution
		dwave[n] = np.median(temp_final_wave[1:] - temp_final_wave[:-1])
	dwave = np.max(dwave)
	final_wave = np.arange(final_wave_start,final_wave_end,dwave)
	print 'Since input dv = 0 -> median resolution (constant) dwave = %f angstrom is used.' % dwave
	return final_wave
	

########################################################################################

def read_fits_spec(fname):
	"""
	Read a fits spectrum written from standard
	pyfits where data are located in the 1st extention HDU
	"""
	import pyfits as pf

	hdu = pf.open(fname)
	wave 	 = hdu[1].data.wavelength
	flux 	 = hdu[1].data.flux
	df 		 = hdu[1].data.df
	df_plus  = hdu[1].data.df_plus
	df_minus = hdu[1].data.df_minus
	return wave,flux,df,df_plus,df_minus

########################################################################################

def printLine():

	print('-------------------------------------------------------------------------------------')

	return

########################################################################################

def printShortLine():

	print('--------------------------')

	return	

def RelativeCalibrationKeyMap():
	print '\n'
	printLine()
	print '?	Show keys map (What is shown here)'
	printLine()
	print 'WINDOW CONTROL KEYS:'
	print '}		shift to right with 5.0 Angstrom'
	print '{		shift to left with 5.0 Angstrom'
	print ']		shift to right with 0.5 Angstrom'
	print '[		shift to left with 0.5 Angstrom'
	print 'shift +/-	Zoom in/out by 0.5 Angstrom'
	print '+/-		Zoom in/out by 0.5 Angstrom'
	print 'T/t		Zoom top by 2e-15'
	print 'B/b		Zoom bottom by 2e-15'
	print 'U/u		Zoom top by 5e-14'
	print 'M/m		Zoom bottom by 5e-14'
	print 'r		replot'
	printLine()
	print 'SPEC CONTROL KEYS:'
	print 'left	shift spec -1 pixel'
	print 'right	shift spec +1 pixel'
	print 'up	scale up'
	print 'down	scale down'

	printLine()
	print 'FITTING SPEC KEYS:'
	print 'enter	create points based on current offset'
	print 'x	Compute chi-squared based on current window'
	print 'D	Delete most recent fitted point in segment b'
	print 'd	Delete most recent fitted point in segment a'
	print '0 to 3:	Assign polynomial fit of nth order in segment a'
	print 'shift + 0 to 3:	Assign polynomial fit of nth order in segment b'
	print 'w	Write a rectified spectrum based on polynomial order fits.'
	printLine()

	print 'Hints: If No lines to fit at all. Enter one point with zero offset'
	printLine()	

########################################################################################
#
# Polynomial Fit Utilities
#
########################################################################################

def best_poly_fit(x_data, y_data, plot_fig=False):
	"""

	Make a best fit polynomial model based of reduced chi squared

	-- this function only fits one segment 
	This function is currently not used in the RelativeCalibration.py
	Returns

		model parameters
	"""
	central_x = np.mean(x_data)
	models = []
	for i in xrange(3):
		with warnings.catch_warnings():
			warnings.filterwarnings('error')
			try: 
				models.append(np.polyfit(x_data - central_x, y_data,i +1))
			except np.RankWarning:
				print "not enought data for %d order" % i
	reduced_chi_squared = []
	if len(x_data) == 2:
		reduced_chi_squared.append(0.0)
	else:
		reduced_chi_squared.append(np.sum((np.polyval(models[0],x_data - central_x) - y_data)**2) / (len(y_data) - 1. - 1.))
		reduced_chi_squared.append(np.sum((np.polyval(models[1],x_data - central_x) - y_data)**2) / (len(y_data) - 2. - 1.))
		reduced_chi_squared.append(np.sum((np.polyval(models[2],x_data - central_x) - y_data)**2) / (len(y_data) - 3. - 1.))

	reduced_chi_squared = np.array(reduced_chi_squared)

	# If there are only 3 data points or less, use 1st order
	if len(x_data) <= 3:
		model = models[0]
	else:
		min_ind = np.argmin(reduced_chi_squared)
		model = models[min_ind]

	if plot_fig:
		fig1 = pl.figure(1)
		pl.plot(x_data,y_data,'o')

		model_x = np.linspace(np.min(x_data),np.max(x_data), 200)
		model_y = np.polyval(model,model_x - central_x)
		return model_x,model_y # return arrays of values to plot. 
	else:
		# parameters array, and central wavelength
		return np.array(model),central_x

########################################################################################

def shift_spec(rectified_wave,scalings,fname):
	"""
	Write out the wavelength corrected spectra with flux 
	scaling. 

	Parameters
	---------------------------------------------------------------------------
	rectified_wave: array
		wavelength array that was corrected from a delta lambda (of some order) shift
	mean_scaling: float
		flux scaling to match

	"""
	scale_sega,scale_segb,break_wavelength = scalings
	wave,flux,error,counts,dfp,dfm = np.loadtxt(fname, unpack = True)

	print len(rectified_wave), len(wave), len(flux)

	inds_a = np.where(wave>=break_wavelength)[0]
	inds_b = np.where(wave<break_wavelength)[0]

	flux  = np.hstack((flux[inds_b]*scale_segb,  flux[inds_a]*scale_sega))
	error = np.hstack((error[inds_b]*scale_segb, error[inds_a]*scale_sega))
	dfp   = np.hstack((dfp[inds_b]*scale_segb, 	dfp[inds_a]*scale_sega))
	dfm   = np.hstack((dfm[inds_b]*scale_segb,  dfm[inds_a]*scale_sega))

	# The exposure names have to be 9 characters long for next line to work properly.
	final_spec_file = open(fname[:-16] +'rect_' + fname[-16:],'w')
	for i in xrange(len(rectified_wave)):
		(final_spec_file.write('%f %.32f %.32f %f %.32f %.32f\n' % (
			rectified_wave[i],
			flux[i],
			error[i],
			counts[i],
			dfp[i],
			dfm[i])))
	final_spec_file.close()

########################################################################################

def WriteFitPoints(fname, central_wave_segb, shift_segb, central_wave_sega, shift_sega):
	f = open(fname[:-16] +'fits_' + fname[-16:],'w')

	f.write('# Segment b\n')
	for i in xrange(len(central_wave_segb)):
		f.write('%f %f\n'  % (central_wave_segb[i], shift_segb[i]))

	f.write('# Segment a\n')
	for i in xrange(len(central_wave_sega)):
		f.write('%f %f\n'  % (central_wave_sega[i], shift_sega[i]))

	f.close()
	return None

########################################################################################

def Perform_Fits(spec2, segment_a,segment_b, break_wavelength):
	"""
	Perform a polynomail fit of order n to correct wavelength array 
	based on dlamvs vs lambda point from segment a and b

	Parameters
	---------------------------------------------------------------------------
	spec2: str
		full path to the spectrum 	
	segment_a: array
		segment_a = central_wave_sega,shift_sega, fit_order_a, scale_sega

	Returns
	---------------------------------------------------------------------------

	"""

	# Read in spectrum to be corrected 
	wave2,flux2,error2 = np.loadtxt(spec2,unpack=True,usecols = [0,1,2])
	
	# Unpack segment; len(central_wave_sega) = number of points/lines fitted = len(shift_sega)
	central_wave_sega,shift_sega, fit_order_a, scale_sega = segment_a
	central_wave_segb,shift_segb, fit_order_b, scale_segb = segment_b

	# Get pivot for the fit. 
	pivot_wave_b = np.mean(central_wave_segb)
	pivot_wave_a = np.mean(central_wave_sega)

	# Get the fitted parameters
	params_segb = np.polyfit(central_wave_segb-pivot_wave_b, shift_segb, fit_order_b)
	params_sega = np.polyfit(central_wave_sega-pivot_wave_a, shift_sega, fit_order_a)

	# Find all the wavelength for segment
	best_fit_wave_segb = [w2 for w2 in wave2 if w2 <= break_wavelength]	
	best_fit_wave_sega = [w2 for w2 in wave2 if w2 > break_wavelength]

	# Get the delta lambda for all the wavelength in each segment 
	best_fit_shift_segb = np.polyval(params_segb, best_fit_wave_segb-pivot_wave_b)
	best_fit_shift_sega = np.polyval(params_sega, best_fit_wave_sega-pivot_wave_a)

	# Add the offset to the original wavelength 
	wave_segb = best_fit_wave_segb + best_fit_shift_segb
	wave_sega = best_fit_wave_sega + best_fit_shift_sega
	new_wave = np.hstack((wave_segb,wave_sega))

	# Write rectified spec
	scalings = np.mean(scale_sega),np.mean(scale_segb),break_wavelength
	shift_spec(new_wave,scalings,spec2)
	#shiftspec(spec2,params_sega,params_segb,pivot_wave_a,pivot_wave_b,0,break_wavelength)

	# Write out points and fit 
	WriteFitPoints(spec2, central_wave_segb, shift_segb, central_wave_sega, shift_sega)
	print 'Written rectified spectrum %s\n'	% spec2[:-16] +'rect_' + spec2[-16:]

	segment_a_fit = best_fit_wave_sega,best_fit_shift_sega,
	segment_b_fit = best_fit_wave_segb,best_fit_shift_segb
	return segment_a_fit, segment_b_fit



def Get_Best_Fit_curve(spec2,segment,break_wavelength,segment_opt = ''):
	wave2,flux2,error2 = np.loadtxt(spec2,unpack=True,usecols = [0,1,2])
	central_wave_seg,shift_seg, fit_order, scale_seg = segment
	pivot_wave = np.mean(central_wave_seg)
	params_seg = np.polyfit(central_wave_seg-pivot_wave, shift_seg, fit_order)
	if segment_opt == 'a':
		best_fit_wave_seg = [w2 for w2 in wave2 if w2 > break_wavelength]
	elif segment_opt == 'b':
		best_fit_wave_seg = [w2 for w2 in wave2 if w2 <= break_wavelength]
	best_fit_shift_seg = np.polyval(params_seg, best_fit_wave_seg-pivot_wave)
	return best_fit_wave_seg, best_fit_shift_seg

########################################################################################

def ChooseFitOrders(fit_order_a,fit_order_b,central_wave_sega,central_wave_segb):

	if len(central_wave_segb) == 1:
		fit_order_b = 0
	elif fit_order_b >= (len(central_wave_segb) - 1):
		print 'Entered order is greater than n-1 points'
		print 'Used linear order instead'
		fit_order_b = 1

	if len(central_wave_sega) == 1:
		fit_order_a = 0;
	elif fit_order_a >= (len(central_wave_sega) - 1):
		print 'Entered order is greater than n-1 points'
		print 'Used linear order instead'
		fit_order_a = 1

	print 'Seg A Fit order %d' % fit_order_a
	print 'Seg B Fit order %d' % fit_order_b

	return fit_order_a,fit_order_b

########################################################################################

class Transition:
    def __init__(self,name,wave,osc_f,gamma,mass):
        self.name  = name
        self.wave  = wave
        self.osc_f = osc_f
        self.gamma = gamma
        self.mass  = mass

########################################################################################

def ReadTransitionData():
	"""
	Dictionary for list of properties of a transition
	"""
	data_file = './data/atomic_data.dat'
	amu = 1.66053892e-24   # 1 atomic mass in grams
	
	name  = np.loadtxt(data_file, dtype=str, usecols=[0])
	wave  = np.loadtxt(data_file, usecols=[1])
	osc_f = np.loadtxt(data_file, usecols=[2])
	gamma = np.loadtxt(data_file, usecols=[3])
	mass  = np.loadtxt(data_file, usecols=[4]) * amu

	transition_dict = {}
	for i in xrange(len(name)):
		transition_dict[str(name[i])] = Transition(name[i],wave[i],osc_f[i],gamma[i],mass[i])
	return transition_dict


########################################################################################


if __name__ == '__main__':
	"""	Test"""
	transition_dict = ReadTransitionData()
	print transition_dict['CIVa'].wave

########################################################################################
