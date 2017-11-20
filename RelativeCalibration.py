########################################################################################
#
# align.py  		(c) Cameron J. Liang
#						University of Chicago
#     				    jwliang@oddjob.uchicago.edu
#     				    Cameron.Liang@gmail.com
#
########################################################################################

"""
This module corrects wavelength-dependent offset between various xd1sum 1D QSO 
spectra to regisiter offsets (shifts) against a reference spectra (also a x1dsum 
1D spectra). All spectra are read from a file list based on the output of x1dsum.py 
These are G130M, G160M graiting and the combined list. With the offsets, users can
choose to a polynomial of degree n (n = 1-3) fit to correct the wavelength array. 
New spectrum will be written with a 'rect_' prefix to the original name. 
"""

########################################################################################

import sys
import os
import shutil
import numpy as np
import matplotlib.pyplot as pl

import Utilities

# Make sure not to use any latex related to speed up rendering
pl.rc('font', family='Bitstream Vera Sans')
pl.rc('text', usetex=False) 

#For interactive key_press_event
pl.switch_backend('tkAgg') 

########################################################################################

def readfile(path_to_filelist, filelist):
	"""
	Reading in list of file names and append the input path to
	each individual file name. 

	Parameters
	---------------------------------------------------------------------------
	path_to_filelist: str
		FUll path to the ascii file and the associated ascii spectra 

	filelist: str
		Full path to the file which contains a list of names of ascii spectra
		produced by x1dsum.py; Each ascii contains at least (wave, flux, error)

	Returns
	---------------------------------------------------------------------------
	files_with_path: array of str
		Full path to each indivial ascii spectra
	
	see also
	---------------------------------------------------------------------------
	x1dsum.py
	"""

	files = np.loadtxt(filelist,dtype = 'str')
	files_with_path = []
	for f in files:
		files_with_path.append(path_to_filelist + '/'+f)
	return files_with_path

########################################################################################

def zoom_region(central_wavelength,dwave,dflux_window_up,dflux_window_down):
	"""
	Define zoom in window in spectrum 
	OUTPUT goes to arguments for pl.xlim, pl.ylim in plottting
	"""
	wave_window_left = central_wavelength - dwave
	wave_window_right = central_wavelength + dwave
	flux_window_up = -1e-14 + dflux_window_up
	flux_window_down = 2e-14 + dflux_window_down
	return [wave_window_left,wave_window_right],[flux_window_up, flux_window_down]

########################################################################################

def wave_res(spec1,line_region):
	"""
	Calculate the resolution in wavelength of a pixel
	element at a given wavelength. (take median of all the resolution elements)
	"""
	wave1,flux1,error1 = np.loadtxt(spec1,unpack=True,usecols = [0,1,2])
	if not line_region:
		pass
	else:
		waves = wave1[(wave1 >= line_region - 0.1) & (wave1 < line_region + 0.1)]

	# Raw resolution of COS at a given wavelength
	delta_wave_per_pix = np.median(wave1[1:] - wave1[:-1]) 
	return delta_wave_per_pix

########################################################################################

def plotting(spec1, spec2, line_region, dwave = 10, 
			 dflux_window_up = 0.0, dflux_window_down=0.0, replot = True):
	"""
	Plots a references spectrum 'spec1' and another to be aligned 'spec2'; Go 
	through the alignment process to output new spec2 with rectified wavelength 
	array based on a polynomial fit correction for the wavelength-dependent offsets.

	Parameters
	---------------------------------------------------------------------------
	spec1: str
		Full path to filename of the reference spectrum
	spec2: str
		Full path to filename of the to-be-aligned spectrum
	line_region: float
		A wavelength region for the plotting window (to start with)
	dwave: float
		plotting window from the centroid of the absorption line
		By default it assumes +/- 10 Angstrom
	dflux_window_up: float
		plotting window in the flux direction; mainly used for initialzation in the 
		argument
	dflux_window_down: float
		plotting window in the flux direction; mainly used for initialzation in the 
		argument
	replot: bool
		By default it makes plotting window based on the intial window parameters. 
		It can be reset to be 'True' after zooming of window of any kind. 

	Returns
	---------------------------------------------------------------------------
	None
	"""

	wave1,flux1,error1 = np.loadtxt(spec1,unpack=True,usecols = [0,1,2])
	wave2,flux2,error2 = np.loadtxt(spec2,unpack=True,usecols = [0,1,2])

	if not line_region:
		both_wave = np.array(wave1,wave2).flatten()
		line_region = np.min(both_wave)
	else:
		waves = wave1[(wave1 >= line_region - 0.1) & (wave1 < line_region + 0.1)]

	fig1 = pl.figure(figsize=(10, 20))

	# segment b points to fit
	central_wave_segb = []; shift_segb = []; scale_segb = []
	
	# segment b best fit 
	break_wavelength = Utilities.break_wave_calc(wave2)
	best_fit_wave_segb = [w2 for w2 in wave2 if w2 <= break_wavelength]
	best_fit_shift_segb= np.zeros(len(best_fit_wave_segb))

	# segment a points to fit
	central_wave_sega = []; shift_sega = []; scale_sega = []
	
	# segment a best fit 
	best_fit_wave_sega = [w2 for w2 in wave2 if w2 > break_wavelength]
	best_fit_shift_sega= np.zeros(len(best_fit_wave_sega))

	############################################################
	# Top Panel: Offset as a function of wavelength
	############################################################
	ax1 = fig1.add_subplot(311)
	dfits_segb, = pl.plot(central_wave_segb,shift_segb, 'bo',label= 'segment b')
	best_fit_line_b, = pl.plot(best_fit_wave_segb, best_fit_shift_segb, c = 'b')

	dfits_sega, = pl.plot(central_wave_sega,shift_sega, 'go',label= 'segment a')
	best_fit_line_a, = pl.plot(best_fit_wave_sega, best_fit_shift_sega, c = 'g')
	ax1.set_xlabel(r'$\lambda$ ($\AA$)')
	ax1.set_ylabel(r'$\Delta \lambda$')
	pl.xlim([min(wave2),max(wave2)])
	pl.legend(loc='best')

	############################################################
	# Middle Panel: Chi Squared as a function of shift in pixels
	############################################################
	ax2 = fig1.add_subplot(312)
	dpix = np.zeros(100); chi2 = np.zeros(100)
	my_chi_waveshift = [0,0]
	my_chi2 = [0,10]
	chi2_line, = pl.step(dpix,chi2)
	my_shift_vline, = pl.plot(my_chi_waveshift,my_chi2,linewidth = 2.0)
	pl.xlim([-60,60])
	ax2.set_xlabel(r'$\Delta$pix')
	ax2.set_ylabel(r'$\chi^2$')

	############################################################
	# Bottom Panel: Spectrum plot
	############################################################
	ax3 = fig1.add_subplot(313)
	ax3.clear()
	pl.step(wave1,flux1,color = 'k', label = str(spec1[-16:-7])) # Reference spec
	line, = pl.step(wave2,flux2,color = 'r', label=str(spec2[-16:-7]))
	pl.legend(loc=1)
	ax3.set_xlabel(r'$\lambda$ ($\AA$)')
	ax3.set_ylabel('Flux')

	if replot:
		pass
	else:
		pl.xlim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		pl.ylim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])

	delta_wave_per_pix = wave_res(spec2,line_region)


	# Initialize dummy variables
	dx = 0.0; 			   	 scale_factor = 1.0
	dwave = 10; 		   	 new_line_region = 0.0
	dflux_window_up = 0.0; 	 dflux_window_down = 0.0

	small_shift = 0.5; 	 	 big_shift = 5.0 		# Shift Spectrum left/right
	zoom = 0.1; 	 		 big_zoom = 0.5 		# Wavelength zoom in/out
	flux_zoom = 2e-15; 		 big_flux_zoom = 5e-14	# Flux zoom in/out
	
	# Shift and Scaling records
	shift_record = []; scale_record = []
	
	########################################################################################

	def shift_spec(event):
		"""
		Interactive click/plotting Event

		Write out rectified spectra from a polynomial fit correction of the 
		wavelength-dependent offset from side effects of the function. 

		Parameters
		---------------------------------------------------------------------------
		event: obj
			Mouse clicks or button press when focus on the plotting window

		Returns
		---------------------------------------------------------------------------
		None
		"""		
		global dx; global scale_factor
		global dwave; global line_region
		global dflux_window_up; global dflux_window_down
		global dpix; global chi2

		global fit_order_b; global fit_order_a
	
		best_fit_wave_segb = []; best_fit_shift_segb = []
		best_fit_wave_sega = []; best_fit_shift_sega = []
		
		#######################################################################
		#														 			  #
		#	WINDOW Control 					 					 			  #
		#																	  #
		#######################################################################
		if event.key == '}':
			line_region += big_shift
			pl.xlim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == '{':
			line_region -= big_shift
			pl.xlim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == ']':
			line_region += small_shift
			pl.xlim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == '[':
			line_region -= small_shift
			pl.xlim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == '-':
			dwave += zoom
			pl.xlim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == '=':
			dwave -= zoom
			pl.xlim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == '_':
			dwave += big_zoom
			pl.xlim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == '+':
			dwave -= big_zoom
			pl.xlim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key =='b':
			dflux_window_up += flux_zoom
			pl.ylim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])
		
		elif event.key =='B':
			dflux_window_up -= flux_zoom
			pl.ylim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])
		
		elif event.key =='t':
			dflux_window_down -= flux_zoom
			pl.ylim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])			
		
		elif event.key =='T':
			dflux_window_down += flux_zoom
			pl.ylim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])
		
		elif event.key =='m':
			dflux_window_up += big_flux_zoom
			pl.ylim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])
		
		elif event.key =='M':
			dflux_window_up -= big_flux_zoom
			pl.ylim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])
		
		elif event.key =='u':
			dflux_window_down -= big_flux_zoom
			pl.ylim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])			
		
		elif event.key =='U':
			dflux_window_down += big_flux_zoom
			pl.ylim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])

		elif event.key =='r':
			dwave = 4; dflux_window_up = 0.0; dflux_window_down = 0.0
			pl.xlim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
			pl.ylim(zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])
		
		#######################################################################
		#														 			  #
		#	Fitting Control 					 					 		  #
		#																	  #
		#######################################################################
		elif event.key == 'enter':
			# Register a aligned point delta lambda based on the given offset. 
			print '\n'
			print 'central Wavelength = %f' % line_region
			print 'Shift = %f' % np.sum(shift_record)
			print 'Scale = %f\n' % np.prod(scale_record)
			
			if line_region >= break_wavelength:
				central_wave_sega.append(line_region)	
				shift_sega.append(np.sum(shift_record))
				scale_sega.append(np.prod(scale_record))
			else:
				central_wave_segb.append(line_region)
				shift_segb.append(np.sum(shift_record))
				scale_segb.append(np.prod(scale_record))
		
		elif event.key == 'x':
			# Compute chi^2 based on data currently displayed on plotting window
			wa = line_region - dwave
			wb = line_region + dwave
			print 'Range for chi^2 = (%f, %f)' % (wa,wb)
			temp_chi2_output = Utilities.chi_squared_calc(wave1,flux1,error1,wave2,flux2,error2,wa,wb)
			dpix = np.array(temp_chi2_output[0])
			chi2 = np.array(temp_chi2_output[1])

		elif event.key == '4' or event.key == '5' or event.key == '6' or event.key == '7':
			# Choose an polynomial fit of order n and write out rectified 
			# spectrum based on the correction to wavelength array
			fit_order_b = int(event.key); fit_order_a = int(event.key)
			if len(central_wave_sega) == 0 or len(central_wave_sega) == 0:
				print '\n'
				print 'Both Segment a and b needs to be have at least 1 point to make fits.\n'
				pass
			else:
				fit_order_a,fit_order_b = Utilities.ChooseFitOrders(fit_order_a,fit_order_b,
											central_wave_sega,central_wave_segb)
			segment_a = central_wave_sega,shift_sega, fit_order_a, scale_sega
			segment_b = central_wave_segb,shift_segb, fit_order_b, scale_segb
			segment_a_fit, segment_b_fit = Utilities.Perform_Fits(spec2,
											segment_a,segment_b, break_wavelength)
			# For Displaying fits on plot. 
			best_fit_wave_sega,best_fit_shift_sega = segment_a_fit
			best_fit_wave_segb,best_fit_shift_segb = segment_b_fit

		elif event.key == ')' or event.key == '!' or event.key == '@' or event.key == '#':
			if len(central_wave_segb) == 0:
				print 'Segment a needs to be have at least 1 point to make fits.\n'
				pass
			else:	
				if event.key == ')': 
					if len(central_wave_segb) >= 1:
						fit_order_b = 0
					else:
						print 'Require at least 1 point to fit 0 order.'
				elif event.key == '!': 
					if len(central_wave_segb) >= 2:
						fit_order_b = 1
					else:
						print 'Require at least 2 point to fit 1 order.'
				elif event.key == '@': 
					if len(central_wave_segb) >= 3:
						fit_order_b = 2
					else:
						print 'Require at least 3 point to fit 2 order.'
				elif event.key == '#': 
					if len(central_wave_segb) >= 4:
						fit_order_b = 3
					else:
						print 'Require at least 4 point to fit 3 order.'

			try:
				fit_order_b
			except NameError:
				print 'Segment a order not assigned.'
			else:
				print 'Segment a Polynomial order = ', fit_order_b
				segment_b = central_wave_segb,shift_segb, fit_order_b, scale_segb
				best_fit_wave_segb,best_fit_shift_segb = Utilities.Get_Best_Fit_curve(spec2,
															segment_b,break_wavelength,'b')

		elif event.key == '0' or event.key == '1' or event.key == '2' or event.key == '3':
			if len(central_wave_sega) == 0:
				print 'Segment b needs to be have at least 1 point to make fits.\n'
				pass
			else:	
				if event.key == '0': 
					if len(central_wave_sega) >= 1:
						fit_order_a = 0
					else:
						print 'Require at least 1 point to fit 0 order.'
				elif event.key == '1': 
					if len(central_wave_sega) >= 2:
						fit_order_a = 1
					else:
						print 'Require at least 2 point to fit 1 order.'
				elif event.key == '2': 
					if len(central_wave_sega) >= 3:
						fit_order_a = 2
					else:
						print 'Require at least 3 point to fit 2 order.'
				elif event.key == '3': 
					if len(central_wave_sega) >= 4:
						fit_order_a = 3
					else:
						print 'Require at least 4 point to fit 3 order.'
			try:
				fit_order_a
			except NameError:
				print 'Fitting b order not assigned.'
			else:
				print 'Segment b Polynomial order = ', fit_order_a
				segment_a = central_wave_sega,shift_sega, fit_order_a, scale_sega
				best_fit_wave_sega,best_fit_shift_sega = Utilities.Get_Best_Fit_curve(spec2,
															segment_a,break_wavelength,'a')

		elif event.key == 'w':
			if fit_order_a >= 0 and fit_order_b >= 0:
				print 'order in a: ', fit_order_a
				print 'order in b: ', fit_order_b
				segment_a = central_wave_sega,shift_sega, fit_order_a, scale_sega
				segment_b = central_wave_segb,shift_segb, fit_order_b, scale_segb

				segment_a_fit, segment_b_fit = Utilities.Perform_Fits(spec2,
												segment_a,segment_b, break_wavelength)
				# For Displaying fits on plot. 
				best_fit_wave_sega,best_fit_shift_sega = segment_a_fit
				best_fit_wave_segb,best_fit_shift_segb = segment_b_fit
			else:
				print 'Fitting orders are not assigned yet. '

		elif event.key == 'D':
			# Delete the last registered offset in segment b
			if not shift_segb:
				print 'Segment b is already empty'
			else:
				central_wave_segb.pop(-1); shift_segb.pop(-1)
				scale_segb.pop(-1)

		elif event.key == 'd':
			# Delete the last registered offset in segment a
			if not shift_sega:
				print 'Segment a is already empty'
			else:
				central_wave_sega.pop(-1); shift_sega.pop(-1)
				scale_sega.pop(-1)

		elif event.key == '?':
			# Show the keys for functions within RelativeCalibration.py
			Utilities.RelativeCalibrationKeyMap()


		#######################################################################
		#														 			  #
		#	Spectrum Control 					 					 		  #
		#																	  #
		#######################################################################
		elif event.key =='right':
			dx += delta_wave_per_pix
			shift_record.append(delta_wave_per_pix)
			my_chi_waveshift[0] = np.sum(shift_record)/delta_wave_per_pix
			my_chi_waveshift[1] = np.sum(shift_record)/delta_wave_per_pix

		elif event.key =='left':
			dx += -delta_wave_per_pix
			shift_record.append(-delta_wave_per_pix)
			my_chi_waveshift[0] = np.sum(shift_record)/delta_wave_per_pix
			my_chi_waveshift[1] = np.sum(shift_record)/delta_wave_per_pix
		
		elif event.key == 'up':
			scale_factor *= 1.1; scale_record.append(1.1)
			
		elif event.key == 'down':
			scale_factor *= 0.9; scale_record.append(0.9)
			
		else:
			dx += 0 
			shift_record.append(0.); scale_record.append(1.)
			my_chi_waveshift[0] = 0

		line.set_xdata(wave2+dx)
		line.set_ydata(flux2*scale_factor)
		
		dfits_segb.set_xdata(central_wave_segb)
		dfits_segb.set_ydata(shift_segb)

		best_fit_line_b.set_xdata(best_fit_wave_segb)
		best_fit_line_b.set_ydata(best_fit_shift_segb)

		chi2_line.set_xdata(dpix)
		chi2_line.set_ydata(chi2)
		
		my_shift_vline.set_xdata(my_chi_waveshift)

		dfits_sega.set_xdata(central_wave_sega)
		dfits_sega.set_ydata(shift_sega)
		
		best_fit_line_a.set_xdata(best_fit_wave_sega)
		best_fit_line_a.set_ydata(best_fit_shift_sega)

		if len(chi2) > 0:
			ax2.set_ylim([min(chi2)- 0.1, max(chi2)+0.1])	

		if len(shift_sega)>0 and len(shift_segb)==0:
			ax1.set_ylim([min(shift_sega)-0.1, max(shift_sega)+0.1])
		elif len(shift_segb)>0 and len(shift_sega)==0:
			ax1.set_ylim([min(shift_segb)-0.1, max(shift_segb)+0.1])
		
		elif len(shift_segb)>0 and len(shift_sega)>0:
			min_b = min(shift_segb); max_b = max(shift_segb)
			min_a = min(shift_sega); max_a = max(shift_sega)
			min_delta_lambda = min(min_a,min_b); max_delta_lambda = max(max_a,max_b)
			
			ax1.set_ylim([min_delta_lambda-0.1, max_delta_lambda+0.1])

		pl.draw()

	civ = fig1.canvas.mpl_connect('key_press_event', shift_spec)
	pl.show()
	
	return

########################################################################################

def print_intro_hints():
	print '--------------------------------------'
	print '----   Focus on Plotting Window	-----'
	print '---- 	 Press ? for keys 	-----'
	print '--------------------------------------'

########################################################################################

def ChooseGrating():		
	"""Choose one greating to align and coadd."""
	grating_option = float(raw_input('Enter keys: G130M = 1, G160M = 2: '))
	
	if grating_option == 1.:
		print '---------------'
		print 'Grating = G130M'
		print '---------------'
		file_list =  qso_path + '/' + 'G130M_filelist'
	
	elif grating_option == 2.:
		print '---------------'
		print 'Grating = G160M'
		print '---------------'
		file_list =  qso_path + '/' + 'G160M_filelist'

	return file_list,grating_option

########################################################################################

def read_filelist():
	"""
	Read in file_list that contains a list 
	of files. 

	OUTPUT: an string array of the names of the files. 
	(not the spectral data)
	"""
	print 'Files in %s: ' % file_list
	spec_array = readfile(qso_path, file_list)
	for i in xrange(len(spec_array)):
		print i+1, ' = ', spec_array[i]
	print '\n'
	return spec_array

########################################################################################

def choose_ref_and_align_spec():
	ref_ind = 0; print 'reference file =  %d' % (ref_ind + 1)
	i = int(raw_input('starting align file =  ')) - 1	
	return ref_ind,i

########################################################################################

if __name__ == '__main__':

	# Full path to the QSO directory
	# It must have subdirectory 'raw' containing x1dsum.fits files.
	qso_path = sys.argv[1]

	if len(sys.argv) != 2:
		print 'python RelativeCalibration.py path_to_QSO_dir'
		sys.exit()
	
	print_intro_hints()
	file_list,grating_option = ChooseGrating()
	spec_array = read_filelist()
	ref_ind, i = choose_ref_and_align_spec()

	# Write reference file into new name like other aligned/rectified files
	new_name = spec_array[ref_ind][:-16] + 'rect_' + spec_array[ref_ind][-16:]
	shutil.copyfile(spec_array[ref_ind],new_name)
	
	# Loop through each file from here
	first_time = True
	while i < len(spec_array):
		if first_time:
			# Print out files names being plotted.
			print '\nSpectra Being plotted:'
			print spec_array[ref_ind]
			print spec_array[i]


		# Define Dummy global variables 
		dx = 0; dwave = 4; scale_factor = 1.0
		dflux_window_up = 0.0; dflux_window_down = 0.0
		dpix = []; chi2 = []
			
		# Obtain plotting parameters...
		temp_wave=np.loadtxt(spec_array[i], usecols=[0])
		line_region=min(temp_wave); break_wavelength=Utilities.break_wave_calc(temp_wave)
		
		# Enter alighment procedure... 
		plotting(spec_array[ref_ind], spec_array[i], line_region, replot=False)
		var=raw_input("(n)ext, or (e)xit program:\n")
		
		if var=='n' or var=='N':
			'Go to next file'
			i = i +1
			first_time = True

		elif var=='e' or var=='E':
			# Exit code without co-adding 
			print 'code exited'
			sys.exit()
		else:
			i=1
			print 'replotting...'
			pass

########################################################################################