########################################################################################
#
# c_specplot.py  	(c) Cameron J. Liang
#						University of Chicago
#     				    jwliang@oddjob.uchicago.edu
#     				    Cameron.Liang@gmail.com
#
########################################################################################

"""
This module applies absolute calibration using Milky Way absorption lines
as zero-point reference. It uses a polynomial fit to correct the wavelength 
array from a wavelength-dependent offset (e.g., dlambda = f(lambda)). 
"""

########################################################################################

import numpy as np
import matplotlib.pyplot as pl

import Utilities

# My default matplotlibrc setting uses tex. tex makes plotting slow
pl.rc('font', family='Bitstream Vera Sans')
pl.rc('text', usetex=False)

#For interactive key_press_event
pl.switch_backend('tkAgg')

########################################################################################

def SelectDataRange(spec,transition_name,redshift,dwave = 10):
	"""
	Select data in some wavepength range dwave centered 
	on the transition wavelength with some redshift. 

	Parameters
	---------------------------------------------------------------------------
	spec: numpy array 
		Multi-dimensional spectral data cube containing wavelength, flux, error, 
		asymmetric flux error (plus and minus)
	transition_name: str
		A name of an ionic transition, according to atom.dat. e.g., 'SiIIa'
	redshift: float
		redshift of an absorber
	dwave: float
		A range of wavelength selected centered on the observed wavelength
		By default, it assumes +/- 10 angstrom from both sides. 
	Returns
	---------------------------------------------------------------------------
	wave: array
		Wavelength array of a selected section 
	flux: array
	sig: array
		error of the flux
	dfp: array
		asymmetric error of flux (dflux_up)
	dfm: array
		asymmetric error of flux (dflux_down)
	"""
	wave,flux,sig,dfp,dfm = spec
	obs_wave = transition_dict[transition_name].wave * (1+redshift)
	inds = np.where( (wave>obs_wave-dwave) & (wave<obs_wave+dwave) )[0]
	return wave[inds],flux[inds],sig[inds],dfp[inds],dfm[inds]

########################################################################################
	
def plot_spec(spec,transition_name, 
				dwave = 10., dflux_window_up = 0.0, dflux_window_down=0.0):
	"""
	Plots a section of a spectrum centered on the transition with 
	centered on the line region 

	Parameters
	---------------------------------------------------------------------------
	spec: numpy array 
		Multi-dimensional spectral data cube containing wavelength, flux, error, 
		asymmetric flux error (plus and minus)	
	transition_name: str
		A name of an ionic transition, according to atom.dat. e.g., 'SiIIa'
	dwave: float
		plotting window from the centroid of the absorption line
		By default it assumes +/- 10 Angstrom
	dflux_window_up: float
		plotting window in the flux direction; mainly used for initialzation in the 
		argument
	dflux_window_down: float
		plotting window in the flux direction; mainly used for initialzation in the 
		argument

	Returns
	---------------------------------------------------------------------------
	centroid_wave: float
		Fitted gaussian centroid of observed wavelength for the input transition; it's 
		used for calculating offset against the rest wavelength
	"""
	wave,flux,error,dfp,dfm = spec
	line_region = np.median(wave)

	fig1 = pl.figure(figsize=(16,8))
	ax1 = fig1.add_subplot(111)
	ax1.clear()
	ax1.set_xlabel(r'$\lambda$ ($\AA$)')
	ax1.set_ylabel(r'$\rm Flux$')

	# Rest central wavelength vertical line
	rest_central_wave = transition_dict[transition_name].wave
	pl.axvline(rest_central_wave,lw=2,ls='--',color = 'b')

	# Original data. 
	pl.step(wave,flux,color = 'k')
	pl.step(wave,error,color = 'r')
	pl.xlim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])

	# Two points setting boundary of data for calculation 
	x_record = []; y_record = []
	points_to_fit, = pl.plot(x_record,y_record,'bo',ms = 8)

	# Gaussian Lines
	gauss_wave = []; gauss_flux = []
	gauss_fit, = pl.plot(gauss_wave,gauss_flux,'b',lw = 1.5)

	# define dummy variables
	dwave = 10;
	big_shift = 0.5;   	   small_shift = 0.1
	big_zoom = 0.5;    	   zoom = 0.1
	dflux_window_up = 0.0; dflux_window_down = 0.0
	flux_zoom = 2e-15; 	   big_flux_zoom = 5e-15

	########################################################################################
	
	def shift_spec(event):
		"""
		Interactive click/plotting Event 

		Parameters
		---------------------------------------------------------------------------
		event: obj
			Mouse clicks or button press when focus on the plotting window

		Returns
		---------------------------------------------------------------------------
		centroid_wave: float
			Fitted gaussian centroid of observed wavelength for the input transition; it's 
			used for calculating offset against the rest wavelength		
		"""
		global transition_name; global centroid_wave
		global line_region; global dwave
		global dflux_window_up; global dflux_window_down
		
		ix, iy = event.xdata, event.ydata
		#######################################################################
		#														 			  #
		#	WINDOW Control 					 					 			  #
		#																	  #
		#######################################################################
		if event.key == '}':
			# Move plotting and spec to right by amount 'big_shift'
			line_region += big_shift
			pl.xlim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == '{':
			# Move plotting and spec to left by amount 'big_shift'
			line_region -= big_shift
			pl.xlim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == ']':
			# Move plotting and spec to right by amount 'small_shift'
			line_region += small_shift
			pl.xlim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == '[':
			# Move plotting and spec to left by amount 'small_shift'
			line_region -= small_shift
			pl.xlim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == '-':
			# Zoom in horizontally in wavelength by amount 'zoom'
			dwave += zoom
			pl.xlim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == '=':
			# Zoom out horizontally in wavelength by amount 'zoom'
			dwave -= zoom
			pl.xlim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == '_':
			# Zoom in horizontally in wavelength by amount 'big_zoom'
			dwave += big_zoom
			pl.xlim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key == '+':
			# Zoom out horizontally in wavelength by amount 'big_zoom'
			dwave -= big_zoom
			pl.xlim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
		
		elif event.key =='b':
			# Zoom (b)ottom in flux by amount 'flux_zoom'
			dflux_window_up += flux_zoom
			pl.ylim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])
		
		elif event.key =='B':
			# UN-Zoom (b)ottom in flux by amount 'flux_zoom'
			dflux_window_up -= flux_zoom
			pl.ylim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])
		
		elif event.key =='t':
			# Zoom in the (t)op in flux by amount 'flux_zoom'
			dflux_window_down -= flux_zoom
			pl.ylim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])			
		
		elif event.key =='T':
			# Zoom out the (t)op in flux by amount 'flux_zoom'
			dflux_window_down += flux_zoom
			pl.ylim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])
		
		elif event.key =='m':
			# Zoom in the bottom in flux by amount 'big_flux_zoom'
			dflux_window_up += big_flux_zoom
			pl.ylim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])
		
		elif event.key =='M':
			# Zoom out the bottom in flux by amount 'big_flux_zoom'
			dflux_window_up -= big_flux_zoom
			pl.ylim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])
		
		elif event.key =='u':
			# Zoom in the top flux by amount 'big_flux_zoom'
			dflux_window_down -= big_flux_zoom
			pl.ylim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])			
		
		elif event.key =='U':
			# Zoom out the top in flux by amount 'big_flux_zoom'
			dflux_window_down += big_flux_zoom
			pl.ylim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])
		
		elif event.key =='r':
			# (r)eplot the spectral window from the original starting point
			dwave = 10; dflux_window_up = 0.0; dflux_window_down = 0.0
			line_region = np.median(wave)
			pl.xlim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[0])
			pl.ylim(Utilities.zoom_region(line_region,dwave,dflux_window_up,dflux_window_down)[1])
		elif event.key == '?':
			print '\n'
			Utilities.printLine()
			print '?	Show keys map (What is shown here)'

			Utilities.printLine()
			print 'WINDOW CONTROL KEYS:'
			print '}		shift to right with 0.5 Angstrom'
			print '{		shift to left with 0.5 Angstrom'
			print ']		shift to right with 0.1 Angstrom'
			print '[		shift to left with 0.1 Angstrom'
			print 'shift +/-	Zoom in/out by 0.5 Angstrom'
			print '+/-		Zoom in/out by 0.1 Angstrom'
			print 'T/t		Zoom top by 1e-15'
			print 'B/b		Zoom bottom by 1e-15'
			print 'U/u		Zoom top by 5e-15'
			print 'M/m		Zoom bottom by 5e-15'
			print 'r		replot'

			Utilities.printLine()
			print 'FITTING SPEC KEYS:'
			print 'a	Add points'
			print 'shift+g	Fit Gaussian'
			Utilities.printLine()

		#######################################################################
		#														 			  #
		#	Fitting Control 					 					 		  #
		#																	  #
		#######################################################################
		elif event.key == 'a':
			# Add 2 Points setting boundary for fitting data

			if len(x_record) < 2:
				x_record.append(ix); y_record.append(iy)
				
			else:
				del x_record[:]; del y_record[:]
				x_record.append(ix); y_record.append(iy)
				
		elif event.key == 'G':
			# Fitt a Gaussian based on the data selected by the two points 
			# in spectrum

			if not x_record: 
				print 'No data selected to fit.'
				pass
			else:
				p1,p2 = np.transpose(np.array([x_record,y_record]))

				x1,y1 = p1; x2,y2 = p2
				temp_spec = Utilities.Select_Data(spec,[x1,x2])
				estimated_cont_level = np.mean((y1,y2))
				gauss_params = Utilities.Fit_Gaussian(temp_spec,estimated_cont_level)

				if gauss_params:
					amplitude,centroid_wave,sigma_width = gauss_params

					# Apparent column density 
					logN = Utilities.ComputeAppColumn(temp_spec,transition_name)

					# Print out results of gaussian fit 
					Utilities.Print_LineInfo(gauss_params,logN,transition_name)

					# Make the plot to show goodness of fit
					gauss_flux = Utilities.Gaussian_function(temp_spec[0],amplitude,centroid_wave,sigma_width)
					gauss_wave = temp_spec[0];
					gauss_fit.set_xdata(gauss_wave)
					gauss_fit.set_ydata(gauss_flux + estimated_cont_level)
				else:
					pass

		points_to_fit.set_xdata(x_record)
		points_to_fit.set_ydata(y_record)

		pl.draw() # needed for instant response. 

	civ = fig1.canvas.mpl_connect('key_press_event', shift_spec)
	pl.show()

	# Exit function properly
	try:
		# Test if centroid_wave exits
		centroid_wave
	except NameError:
		print 'No Gaussian fitted to any lines; returning...'
		return np.nan
	else:
		if abs(centroid_wave - transition_dict[transition_name].wave) > 2: 
			print 'Note: %f has a large deviation > 2 angstrom' % centroid_wave
			print 'Results not recorded.'
			return np.nan
		else:
			return centroid_wave

########################################################################################

if __name__ == '__main__':

	# Define some constants. 
	dwave = 10; dflux_window_up = 0; dflux_window_down = 0
	transition_dict = Utilities.ReadTransitionData()

	# Read Input Spectrum
	spec_path = raw_input('Full path to QSO spectrum directory:\n')
	fname = raw_input('File name of spectrum:\n')
	input_fname = spec_path + '/' + fname
	spec = np.loadtxt(input_fname,unpack=True)
	
	# Define arrays for collection of fits 
	transitions = np.array([]); transition_rest_wave = []; 
	obs_centroid = []; dlambda = [];  
	rectifing = True
	while rectifing:
		# Ask users to enter option
		# Letters in parathesese are the keys to enter in terminal.
		var = raw_input('(f)it new line, choose offset to (d)elete, (p)rint current offsets,'
						'(w)rite new spec, (e)xit program:\n')
		
		if var == 'f' or var == 'F':
			Utilities.PrintIonsNames()
			print 'Wavelength range of data = [%f\t%f]' % (np.min(spec[0]),np.max(spec[0]))
			print 'i.e., do not use lines outside of the range.\n'
			transition_name,redshift,line_region = Utilities.UserInputs()
			spec_chunk = SelectDataRange(spec,transition_name,redshift)
			centroid_wave =  plot_spec(spec_chunk,transition_name)
			
			# Record the fits if centroid_wave was successfully obtained.
			if centroid_wave:
				if np.isnan(centroid_wave):
					continue
				elif centroid_wave in obs_centroid:
					print 'Observed centroid identical to previous ones... Not recorded.'
					print 'Probably because centroid has not been updated from previous fit.'
					continue
				else:
					rest_wave = transition_dict[transition_name].wave
					transitions = np.append(transitions,transition_name)
					transition_rest_wave.append(rest_wave)
					obs_centroid.append(centroid_wave)
					dlambda.append(rest_wave-centroid_wave)
		
		# Print current fitted offsets for viewing purposes. 
		elif var == 'p' or var == 'P':
			if len(dlambda) == 0:
				print 'No lines and shifts recorded\n'
				continue
			Utilities.print_dwave_fit(transitions,transition_rest_wave,obs_centroid,dlambda)
			Utilities.plot_dwave_fit(transition_rest_wave,dlambda,np.median(spec[0]))
		
		# Delete any of the current offsets if needed. 
		elif var == 'd' or var == 'D':
			Utilities.print_dwave_fit(transitions,transition_rest_wave,obs_centroid,dlambda)
			
			transition_to_delete = raw_input('Transition name to delete: ')
			if transition_to_delete in transitions:
				index = np.where(transitions == transition_to_delete)[0]
				transitions = np.delete(transitions,index); transition_rest_wave.pop(index)
				obs_centroid.pop(index); dlambda.pop(index)
				print 'Deleted fit to %s\n' % transition_to_delete
				Utilities.print_dwave_fit(transitions,transition_rest_wave,obs_centroid,dlambda)
			
			elif transition_to_delete not in transitions:
				print '\'%s\' not in the list; back to main command...\n' % transition_to_delete
				continue
			
			else:
				print 'Nothing entered; back to main command...\n'
				continue

		# Write a rectified spectrum based of the offset and a chosen fit. 
		elif var == 'w' or var == 'W':
			if len(dlambda) == 0:
				print 'No lines and shifts recorded for fitting\n'
				continue
			pivot_wave = np.median(spec[0]);
			fit_order = Utilities.determine_fit_order(dlambda)
			fit_params = np.polyfit(transition_rest_wave-pivot_wave,dlambda,fit_order)
			Utilities.write_shiftspec(spec,fit_params,spec_path)
		
		# Exit program properly. 
		elif var == 'e' or var == 'E':
			print 'Exited program\n'
			exit()
		
		else:
			print '\'%s\' is not an option; enter key again.\n' % var

########################################################################################
