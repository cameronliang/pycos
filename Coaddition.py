########################################################################################
#
# Coaddition.py  	(c) Cameron J. Liang
#						University of Chicago
#     				    jwliang@oddjob.uchicago.edu
#     				    Cameron.Liang@gmail.com
#
########################################################################################

"""
Generic coaddition code
"""

import os
import re
import sys
import time
import numpy as np
from scipy.interpolate import interp1d
import pylab as pl
import Utilities

########################################################################################

# Constants
c = 299792.458	# speed of light in km/s 

########################################################################################

def ReadFileList(path_to_filelist, filelist):
	"""
	Reading in list of file names 

	filelist is an ascii filelist
	that contains a list of ascii spectra file,
	each with wave, flux, error
	"""
	filelist  = path_to_filelist + '/' + filelist
	files = np.loadtxt(filelist,dtype = 'str')
	
	print filelist
	if len(np.atleast_1d(files)) == 1:
		spec_list = [path_to_filelist  + '/' + str(files)]

	elif len(np.atleast_1d(files)) > 1:
		spec_list = []
		print 'List of files to coadd:'
		for f in files:			
			spec_list.append(path_to_filelist + '/' + 'rect_' + f)
			print 'rect_' + f
	spec_list = np.array(spec_list)
	return spec_list

########################################################################################

def ReadCOSx1dsumSpectrum(filename):
	"""
	filename with full path
	Purporse is to have other variation 
	of files and differnet way of reading in. 
	"""
	wave,flux,dfp,dfm = np.loadtxt(filename,unpack=True,usecols=[0,1,4,5])
	return np.array([wave,flux,dfp,dfm])

########################################################################################

def ReadSpecCube_x1dsum(spec_list):
	"""
	Produce a 3D cube, by reading in the all the spectra. 
	spec_list = 1D array with full path to spectra for coadding


	OUTPUT: specs[spec_index][0=wave, 1=flux, 2= dfp, 3=dfm][pixel number]

	eg.
	wave_array of fist spec = specs[0][0]
	flux_array of fist spec = specs[0][1]
	"""
	specs = []
	for f in spec_list:
		specs.append(ReadCOSx1dsumSpectrum(f))
	return specs

########################################################################################

def ReadSpecCube_Normal(spec_list):
	"""
	Reads all other ascii spectrum (except x1dsum)
	assumed: spec = [wave,flux,error,dfp,dfm]
	"""
	def Read_NormalSpec(filename):
		# collect [wave,flux,dfp,dfm] only
		# so that the specs have same dimension as output from ReadCOSx1dsumSpectrum
		return np.loadtxt(filename,unpack=True,usecols=[0,1,3,4])

	specs = []
	for f in spec_list:
		specs.append(Read_NormalSpec(f))
	return specs

########################################################################################

def ComputeWeights(spec_cube):
	"""
	Compute a weight for each spectrum
	Definiton of weight used: weight = (Median S/N) ** 2
	
	Parameters
	---------------------------------------------------------------------------	
	spec_cube: array
		Multi-dimensional array of spectral data 

	Returns
	---------------------------------------------------------------------------	
	weights: array
		Normalized weights based on squared of the signal-to-noise for weighting of coadd. 
		length = number of spectra

	"""
	Num_file = len(spec_cube)
	weights = np.zeros(Num_file)
	for n in xrange(Num_file):
		# Remove all negative values from SNR array
		signal = spec_cube[n][1] # flux 
		noise = spec_cube[n][2]  # error 
		SNR = np.clip(signal/noise, 0, np.max(signal/noise))
		median_SNR = np.median(np.trim_zeros(np.sort(SNR)))
		#median_SNR = np.median(1/noise) 
		weights[n] = median_SNR**2

	# Normalize 1D weights
	weights = weights / np.sum(weights)#np.linalg.norm(weights)
	return weights

########################################################################################

def ComputeAverages(array, weights, opt=''):
	"""Weighted Average"""
	axis_number = 0
	if opt == 'simple':
		return np.mean(array,axis = axis_number)
	elif opt == 'simple_weighted':
		return np.average(array, weights = weights,axis = axis_number)
	elif opt == 'other_way_of_weighting':
		return None
	else:
		print 'option not reached.'
		print 'code exited'
		exit()

########################################################################################

def ComputeWeightedError(weighted_mean, flux_array, error_array, weights):
	"""Unbiased Sample Weighted Standard Deviation"""
	V1 = np.sum(weights); V2 = np.sum(weights**2)
	biased_weighted_error = np.sqrt(np.sum(weights*(flux_array - weighted_mean)**2))
	biased_correction_factor = np.sqrt(V1  - V2/V1)

	return biased_weighted_error / biased_correction_factor

########################################################################################

def ComputeSimpleWeightedError(error_array, weights):
	"""
	Simple Weighted Standard Deviation 
	(assumed uncorrelated points - not entirely correct)
	"""
	top = np.sqrt(np.sum((weights*error_array)**2))
	bottom = np.sum(weights)
	return top/bottom

########################################################################################

def find_nearest_neighbor(array,value):
	"""
	Find the two indices of array of the nearest neightbor 
	of the input value
	"""
	ind = (np.abs(array-value)).argmin()

	if ind == 0 or ind == len(array)-1:
		return np.array([ind,ind])
	else:
		if array[ind] <= value:
			ind2 = ind + 1
		else:
			ind2 = ind - 1
	return np.array([ind,ind2])

########################################################################################

def Find_2points(inds1,inds2,wave,flux,error):
	#wave = np.array(wave); flux = np.array(flux); error = np.array(error);
	p1 = np.array([wave[inds1],flux[inds1],error[inds1]])
	p2 = np.array([wave[inds2],flux[inds2],error[inds2]])
	return p1,p2

########################################################################################

def ComputeError_Interp(x1,y1,dy1,x2,y2,dy2,xp):
	"""
	Flux error from the interpolation coadd scheme
	Compute the propogated error from m and b
	(yp = m xp + b)
	"""
	#x1,y1,dy1 = p1; x2,y2,dy2 = p2

	x_pivot = 0.5*(x1+x2)
	x1 = x1 - x_pivot; x2 = x2 - x_pivot
	xp = xp - x_pivot

	if x2 == x1:
		return np.nan
	else:
		dm_sqr = (dy1**2 + dy2**2)/(x2-x1)**2
		db_sqr = ((1.+x1/(x2-x1))*dy1)**2 + ((x1/(x2-x1))*dy2)**2
		return np.sqrt((xp**2)*dm_sqr + db_sqr)


########################################################################################

def Coadd_interp(spec_cube,weights,dv,weighting_option = ''):
	"""
	Coaddition function based on interpolation scheme. A given pixel is splited between
	the neightboring 2 pixels.

	Parameters
	---------------------------------------------------------------------------	
	spec_cube: array
		Multi-dimensional array of spectral data
		spec_cube[i][0] = wave
		spec_cube[i][1] = flux
		spec_cube[i][2] = dfp
		spec_cube[i][3] = dfm
	weights: array
		Normalized weights based on squared of the signal-to-noise for weighting of coadd.
		length = number of spectra
	dv: float
		Resolution element ; [dv] = km/s; If dv = 0 -> choose native resolution
		see also: Create_Final_NativeWave from Utilities
	weighting_option: str
		Not used in this function -- added only for uniformity from the other Coadd function

	Returns
	---------------------------------------------------------------------------
	spec: array
		final spectrum file with columns: final_wave, final_flux,final_df,
		final_dfp, final_dfm
	
	See also
	---------------------------------------------------------------------------
	Coadd_func 
		A different coadd scheme with re-grouping
	"""
	
	Num_file = len(spec_cube) 	# Number of files/spectra
	interp_order = 'linear'     # Interpolation order

	# Determine the start/end of the wavelength array from min/max wavelength of all spectra
	temp_wave_start = np.zeros(Num_file); temp_wave_end = np.zeros(Num_file)
	for n in xrange(Num_file):
		temp_wave_start[n] = np.min(spec_cube[n][0])
		temp_wave_end[n]   = np.max(spec_cube[n][0])
	final_wave_start = np.min(temp_wave_start)
	final_wave_end   = np.max(temp_wave_end)

	# Choose resolution scheme and define final wavelength array
	if dv == 0:
		# Choose native resolution based on the original wavelength array
		final_wave = Utilities.Create_Final_NativeWave(spec_cube,combine_grating)
		TOTAL_NUMBER_PIXEL = len(final_wave)
	else:
		# Constant velocity resolution (dv)
		TOTAL_NUMBER_PIXEL = np.int(np.log10(final_wave_end/final_wave_start) 
									/ np.log10(1 + dv/c) + 0.5)
		array_index = np.arange(0,TOTAL_NUMBER_PIXEL,1)
		final_wave  = final_wave_start * ((1 + dv/c)**array_index) # one can work out this formula

	final_flux  = np.zeros((Num_file,TOTAL_NUMBER_PIXEL))
	final_dfp 	= np.zeros((Num_file,TOTAL_NUMBER_PIXEL))
	final_dfm 	= np.zeros((Num_file,TOTAL_NUMBER_PIXEL))
	final_weights = np.zeros((Num_file,TOTAL_NUMBER_PIXEL))

	# Interpolate each spectrum to get flux/error at new wavelength
	for n in xrange(Num_file):
		# interpolate functions wave,flux,error in individual spectrum 
		flux_interp = interp1d(spec_cube[n][0],spec_cube[n][1], kind = interp_order,
								bounds_error=False, fill_value=np.nan)	
		
		# Evaluate interpolated flux and error at new wavelengths. 
		final_flux[n] = flux_interp(final_wave) 

		# Calculate propagated error from linear fit
		temp_min_ind = np.searchsorted(final_wave, np.min(spec_cube[n][0]) )
		temp_max_ind = np.searchsorted(final_wave, np.max(spec_cube[n][0]) )
		temp_wave = final_wave[temp_min_ind:temp_max_ind]
		inds2 = np.searchsorted(np.sort(spec_cube[n][0]),temp_wave) 
		inds1 = inds2-1; inds1[0] = 0

		p1,p2 = Find_2points(inds1,inds2,spec_cube[n][0],spec_cube[n][1],spec_cube[n][2])
		final_dfp[n][temp_min_ind:temp_max_ind] = ComputeError_Interp(p1[0],p1[1],p1[2],
																	p2[0],p2[1],p2[2],temp_wave)

		p1,p2 = Find_2points(inds1,inds2,spec_cube[n][0],spec_cube[n][1],spec_cube[n][3])
		final_dfm[n][temp_min_ind:temp_max_ind] = ComputeError_Interp(p1[0],p1[1],p1[2],
																	p2[0],p2[1],p2[2],temp_wave)

		# Create the same weights per spectrum
		final_weights[n] = weights[n]*np.ones_like(final_flux[n])

	# mask out nans	
	final_flux = np.ma.masked_array(final_flux, np.isnan(final_flux))
	final_dfp  = np.ma.masked_array(final_dfp, np.isnan(final_dfp))
	final_dfm  = np.ma.masked_array(final_dfm, np.isnan(final_dfm))


	# mask out zeros (from the array initialization)
	final_flux = np.ma.masked_array(final_flux, final_flux==0)
	final_dfp  = np.ma.masked_array(final_dfp, final_dfp==0)
	final_dfm  = np.ma.masked_array(final_dfm, final_dfp==0)

	# Mask the same elements where the spec is masked. 
	final_weights  = np.ma.masked_array(final_weights, np.ma.getmask(final_flux))
	final_weights  = np.ma.masked_array(final_weights, np.ma.getmask(final_dfp))
	final_weights  = np.ma.masked_array(final_weights, np.ma.getmask(final_dfm))
	
	# Re-normalize weights after masking
	weight_normalization = np.sum(final_weights,axis=0)
	final_weights = final_weights / weight_normalization

	# Average over spectra (over n); 
	final_flux 	= np.ma.average(final_flux, weights=final_weights, axis=0)

	# Same as ComputeSimpleWeightedError function 
	final_dfp = np.sqrt(np.ma.sum( (final_weights*final_dfp)**2, axis = 0))
	final_dfm = np.sqrt(np.ma.sum( (final_weights*final_dfm)**2, axis = 0))
	final_df 	= 0.5*(final_dfp+final_dfm)

	return np.array([final_wave, final_flux,final_df, final_dfp, final_dfm])

########################################################################################

def Coadd_func(spec_cube, weights, dv, weighting_option = ''):
	"""
	Coaddition function based on rebinning scheme.

	Parameters
	---------------------------------------------------------------------------	
	spec_cube: array
		Multi-dimensional array of spectral data
		spec_cube[i][0] = wave
		spec_cube[i][1] = flux
		spec_cube[i][2] = dfp
		spec_cube[i][3] = dfm
	weights: array
		Normalized weights based on squared of the signal-to-noise for weighting of coadd. 
		length = number of spectra
	dv: float
		Resolution element ; [dv] = km/s; If dv = 0 -> choose native resolution 
		see also: Create_Final_NativeWave from Utilities
	weighting_option: str
		Not used in this function -- added only for uniformity from the other Coadd function

	Returns
	---------------------------------------------------------------------------
	spec: array
		final spectrum file with columns: final_wave, final_flux,final_df,
		final_dfp, final_dfm
	
	See also
	---------------------------------------------------------------------------
	Coadd_interp 
		A different coadd scheme with interpolation between flux
	"""
	Num_file = len(spec_cube) 	# Number of files/spectra

	# Determine the start/end of the wavelength array from min/max wavelength of all spectra
	temp_wave_start = np.zeros(Num_file); temp_wave_end = np.zeros(Num_file)
	for i in xrange(Num_file):
		temp_wave_start[i] = np.min(spec_cube[i][0])
		temp_wave_end[i] = np.max(spec_cube[i][0])
	final_wave_start = np.min(temp_wave_start)
	final_wave_end = np.max(temp_wave_end)

	if dv == 0:
		final_wave = Utilities.Create_Final_NativeWave(spec_cube,combine_grating)
		TOTAL_NUMBER_PIXEL = len(final_wave)
	else:
		TOTAL_NUMBER_PIXEL = np.int(np.log10(final_wave_end/final_wave_start)
											 / np.log10(1 + dv/c) + 0.5)
		array_index = np.arange(0,TOTAL_NUMBER_PIXEL,1)
		final_wave  = final_wave_start * ((1 + dv/c)**array_index) # one can work out this formula

	final_flux  = np.zeros(TOTAL_NUMBER_PIXEL)
	final_dfp 	= np.zeros(TOTAL_NUMBER_PIXEL)
	final_dfm 	= np.zeros(TOTAL_NUMBER_PIXEL)

	# Copy each weight of spec into length of spec
	all_weights = []
	all_wave = []; all_flux = []
	all_dfp = []; all_dfm = []
	for i in xrange(Num_file):
		all_weights.append(np.ones(len(spec_cube[i][1]))*weights[i])
		all_wave.append(spec_cube[i][0]); all_flux.append(spec_cube[i][1]); 
		all_dfp.append(spec_cube[i][2]); all_dfm.append(spec_cube[i][3]); 

	all_wave = np.hstack(all_wave); all_flux = np.hstack(all_flux)	
	all_dfp = np.hstack(all_dfp); all_dfm = np.hstack(all_dfm)
	
	all_weights = np.hstack(all_weights)
	
	# Sort everything by wavelength 
	all_wave,all_flux,all_dfp,all_dfm,all_weights = np.array(zip(*sorted(zip(all_wave,all_flux,all_dfp,all_dfm ,all_weights))))


	# Compute weighted average; i.e stack
	for j in xrange(1,TOTAL_NUMBER_PIXEL):
		index = np.where((all_wave >= final_wave[j-1]) 
						& (all_wave < final_wave[j]))[0]

		if len(index) == 0:
			# if no data in this bin, set final = Nan
			final_flux[j] = float('NaN')
			final_dfp[j] = float('NaN')
			final_dfm[j] = float('NaN')
		else:
			final_flux[j] = ComputeAverages(all_flux[index], all_weights[index], 
											opt=weighting_option)
			final_dfp[j] = ComputeSimpleWeightedError(all_dfp[index],all_weights[index])
			final_dfm[j] = ComputeSimpleWeightedError(all_dfm[index],all_weights[index])

	good_index = np.where(np.logical_not(np.isnan(final_dfp)))[0]
	
	final_wave = final_wave[good_index]; final_flux = final_flux[good_index]
	final_dfp = final_dfp[good_index]; final_dfm = final_dfm[good_index]
	final_df = 0.5*(final_dfp+final_dfm)

	return np.array([final_wave, final_flux,final_df, final_dfp, final_dfm])

########################################################################################

def WriteSpectrum(wave,flux,error,file_option=''):
	if file_option == 'ascii':

		final_output = path_to_spec + '/' + output_filename + '.spec'
		f = open(final_output,'w')
		f.write('# wavelength\tflux\terror\n')
		for i in xrange(len(wave)):
			f.write('%f\t%.32f\t%.32f\n' % (wave[i], flux[i],error[i]))
		f.close()
	
	elif file_option =='fits':
		import pyfits 
		col1 = pyfits.Column(name='wavelength',format='E',array=wavelength)
		col2 = pyfits.Column(name='flux',format='E',array=flux)
		cols = pyfits.ColDefs([col1,col2])
		tbhdu = pyfits.BinTableHDU.from_columns(cols)
		tbhdu.writeto(fname + '.fits')

	print 'Written: ', final_output
	
	return None

########################################################################################

def WriteSpectrumAsymmetricError(path_to_spec, output_filename, data_array ,file_option=''):
	wave, flux, df, dfp, dfm = data_array
	if file_option == 'ascii':
		final_output = path_to_spec + '/' + output_filename + '.spec'
		f = open(final_output,'w')
		f.write('# wavelength\tflux\tdf_mean\tdf_plus\tdf_minus\n')
		for i in xrange(len(wave)):
			f.write('%f\t%.32f\t%.32f\t%.32f\t%.32f\n' % (wave[i], flux[i],df[i], dfp[i], dfm[i]))
		f.close()
		print 'Written: ', final_output

	elif file_option =='fits':
		import pyfits
		final_output = path_to_spec + '/' + output_filename
		col1 = pyfits.Column(name='wavelength',format='E',array=wave)
		col2 = pyfits.Column(name='flux',format='E',array=flux)
		col3 = pyfits.Column(name='df',format='E',array=df)
		col4 = pyfits.Column(name='df_plus',format='E',array=dfp)
		col5 = pyfits.Column(name='df_minus',format='E',array=dfm)
		cols = pyfits.ColDefs([col1,col2,col3,col4,col5])
		tbhdu = pyfits.BinTableHDU.from_columns(cols)
		tbhdu.writeto(final_output + '.fits')

	return None

########################################################################################

def main(path_to_spec,spec_filelist,output_filename,dv):
	"""Carry out interpolation based on various inputs"""
	
	# Read in file list of spectra to co-add
	# Read in the cube of spectra
	if isinstance(spec_filelist,list):
		spec_cube = ReadSpecCube_Normal(spec_filelist)
	else:
		spec_list = ReadFileList(path_to_spec, spec_filelist)
		spec_cube = ReadSpecCube_x1dsum(spec_list)

	# Compute weights of weighted-average of coaddition
	weights = ComputeWeights(spec_cube)

	# Coadd the spectra
	var = raw_input(" (a) Rebinning coadd or (b) Interpolated coadd:\n")
	if var == 'a':
		data_array = Coadd_func(spec_cube,weights,dv,weighting_option = 'simple_weighted')
	elif var == 'b':
		data_array = Coadd_interp(spec_cube,weights,dv,weighting_option = 'simple_weighted')
	else:
		print '%s Not valid option; exit code.' % var
		exit()

	# Write to file
	WriteSpectrumAsymmetricError(path_to_spec, output_filename, data_array,file_option = 'ascii')

def main_interp(path_to_spec,spec_filelist,output_filename,dv):
	"""Carry out interpolation based on various inputs"""
	
	# Read in file list of spectra to co-add
	# Read in the cube of spectra
	if isinstance(spec_filelist,list):
		spec_cube = ReadSpecCube_Normal(spec_filelist)
	else:
		spec_list = ReadFileList(path_to_spec, spec_filelist)
		spec_cube = ReadSpecCube_x1dsum(spec_list)

	# Compute weights of weighted-average of coaddition
	weights = ComputeWeights(spec_cube)

	# Coadd the spectra
	data_array = Coadd_interp(spec_cube,weights,dv,weighting_option = 'simple_weighted')

	# Write to file
	WriteSpectrumAsymmetricError(path_to_spec, output_filename, data_array,file_option = 'ascii')


########################################################################################

def interactive_main():
	"""Options to give filelists to coadd. """
	if len(sys.argv) == 1:
		path_to_spec 	= raw_input('Full path to spectrum:\n')

		opt = raw_input('(a) G130M, (b) G160M, (c) manually enter file name\n')
		if opt == 'a':
			spec_filelist = 'G130M_filelist'
		elif opt == 'b':
			spec_filelist = 'G160M_filelist'
		elif opt == 'c':
			spec_filelist 	= raw_input('File containing lists of exposures:\n')
		else:
			print '%s is not an option; exit program.\n' % opt
			exit()
		output_filename = raw_input('Output file name (w/o extension):\n')

		opt = raw_input('(a) Native Resolution, (b) 7.5 km/s, (c) 15 km/s\n')
		if opt == 'a':
			dv = 0. # When dv=0 is chosen, native resolution is used. 
		elif opt == 'b':
			dv = 7.5
		elif opt == 'c':
			dv = 15.0
		else:
			print '%s is not an option; exit program.\n' % opt
			exit()
	elif len(sys.argv) == 5:
		path_to_spec 	= sys.argv[1]			# full Path to spectrum 
		spec_filelist 	= sys.argv[2]			# List of files containing the spectra for coadd
		output_filename = sys.argv[3]			# Final output spectrum filename
		dv 				= float(sys.argv[4]) 	# Resolution of spectrum, i.e binning
		main_interp(path_to_spec,spec_filelist,output_filename,dv)
		exit()
	else:
		print '\n'
		Utilities.printLine()
		print 'Run code as below:'
		print 'python Coaddition.py path_to_spectrum spec_list output_filename(w/o extension) resolution_in_kms\n'
		print 'or simply:'
		print 'python Coaddition.py'
		Utilities.printLine()
		print '\n'
		exit()

	main(path_to_spec,spec_filelist,output_filename,dv)
	return 

########################################################################################

def manual_coadd():
	"""
	Manually enter full path of the file name to coadd. 
	"""
	spec_list = []
	keep_going = True
	while keep_going:
		var = raw_input('Full Path to spectrum with file name, or\n'
						'Press (d) when done:\n')
		if var =='d':
			if len(spec_list) == 0:
				print 'Array of filenames is empty; exit code'
				exit()

			path_to_spec = raw_input('\nOutput path:\n')
			output_filename = raw_input('\nOutput filename (without extension):\n')
			dv = float(raw_input('\nResolution dv(km/s) = \n'
								'(Choose 0 km/s for native resolution)\n'))

			main(path_to_spec,spec_list,output_filename,dv)
			exit()
		else:
			spec_list.append(var)

########################################################################################

if __name__ == '__main__':
	ComputeError_Interp = np.vectorize(ComputeError_Interp)
	combine_grating = False	
	if len(sys.argv) > 2:
		interactive_main()
		exit()
	var = raw_input('(a) Read from file lists\n'
					'(b) Combine two gratings?\n')
	if var == 'a':
		interactive_main()
	elif var == 'b':
		combine_grating = True
		manual_coadd()

########################################################################################
