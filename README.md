# pycos
###################################################################################
#
#   		(c) Cameron J. Liang
#		    University of Chicago
#     		    jwliang@oddjob.uchicago.edu
#     	            Cameron.Liang@gmail.com
#
###################################################################################


 PyCOS: Python Cosmic Origin Spectrograph Reduction Pipeline
 Copyright (C) 2015 Cameron J. Liang
 
Permission is hereby granted, free of charge, to any person obtaining a copy of this 
software and associated documentation files (the "Software"), to deal in the Software 
without restriction, including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons 
to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or 
substantial portions of the Software.


THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.

###################################################################################


PyCOS assumes standard Python libraries: numpy, scipy, matplotlib and pyfits. 
I recommend installing the Enthought Canopy which comes with all the scientific 
libraries you need. You can find it here: https://www.enthought.com/products/canopy/


###################################################################################

How TO Use the Reduction Pipeline PyCOS

1. Run Process_x1dsum.py as follow:

   python Process_x1dsum.py full_path_to_qso_directory

   — This produces x1dsum ascii files with poisson flux error from x1dsum.fits files

2. Run RelativeCalibration.py:

   python RelativeCalibration.py full_path_to_qso_directory
   
   — Align different exposures to a reference, and produce the relative wavelength 
   calibration files with prefix “rect_”

3. Run Coaddition.py for each grating.

   python Coaddition.py

   — No additional command lines arguments at the beginning. Supply path and other   	     
   information as prompted.

4. Run AbsoluteCalibration.py for absolute calibration. 

   python AbsoluteCalibration.py

   — This produce final rectified file. 


5. Run Coaddition.py to combine both gratings

###################################################################################


See ./docs/pycos.pdf for more documentation. If you use PyCOS, 
please cite Liang & Chen 2014.

###################################################################################

