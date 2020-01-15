# ---------------------------------------------------------------------------------
#   minispice-> signalTools.py
#   Copyright (C) 2019 Michael Winters
#   github: https://github.com/mesoic
#   email:  mesoic@protonmail.com
# ---------------------------------------------------------------------------------
#   
#   Permission is hereby granted, free of charge, to any person obtaining a copy
#   of this software and associated documentation files (the "Software"), to deal
#   in the Software without restriction, including without limitation the rights
#   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#   copies of the Software, and to permit persons to whom the Software is
#   furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
#   copies or substantial portions of the Software.
#   
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.
#

#!/usr/bin/env python 
from scipy import signal
import numpy as np

# Create signal tools namespace
class signalTools:

	def __init__(self):

		pass

	# Method to generate a test pulse for transient simulations
	def pulse(self, config):

		# Construct time array 
		time  = np.linspace(0, config["period"], config["npoints"], endpoint=False)

		# Cache the sampling interval
		delta = time[1] - time[0]

		# Use the signal library to create a pulse signal
		waveform = (config["amplitude"] / 2) * ( signal.square(2 * np.pi * (1.0/config["period"]) * time - np.pi / 2.0, config["duty"]) + 1.0)

		# Return waveform, time, and sampling interval (delta)
		_signal = {
			"waveform"	: waveform,
			"time"		: time,
			"delta" 	: delta
		} 

		return _signal

	# Method to generate a sin wave (one cycle)
	def sinwave(self, config):

		# Construct time array 
		time  = np.linspace(0, config["period"], config["npoints"], endpoint=False)

		# Cache the sampling interval
		delta = time[1] - time[0]

		# Build sin wave over period
		waveform = [ config["amplitude"] * np.sin(  2 * np.pi * (1.0/config["period"]) * _ ) for _ in time ]

		# Return waveform, time, and sampling interval (delta)
		_signal = {
			"waveform"	: waveform,
			"time"		: time,
			"delta" 	: delta
		} 

		return _signal


# A class with two methods which to calclates the Discrete Fourier Transform
# and Inverse Discrete Fourier Transform matricies. 
class DiscreteFourierTransform:

	def __init__(self):

		pass

	# Class to hold the transform object
	class Transform: 

		# User defines frequency and number of harmonics
		def __init__(self, freq = 1e9, n = 20):

			# Calculate n harmonics for initialized frequency
			self.harmonics = {
				"n"	 	: int(n),			  
				"freq"  : float(freq),
				"omega" : [ _ * ( 2 * np.pi * freq ) for _ in range( int(n) ) ],
				"tau"   : [ _ / (  freq * float(n) ) for _ in range( int(n) ) ],	
			} 

			# Build discrete fourier transform matrices
			self.build()

		# Return an empty matrix of zeros for matrices
		def zeros(self):

			return np.zeros(
				shape = ( self.harmonics["n"], self.harmonics["n"] ), 
				dtype='complex'
			)

		# Build the transform matrices   
		def build(self):	

			# Cache number of harmonics
			n = float( self.harmonics["n"] )

			# Twiddle factor
			W = np.exp( np.complex( 0,  (-2.0 * np.pi / n ) ) )

			# Calculate the discrete fourier transform matrix
			self.dft, self.idft = self.zeros(), self.zeros()
			
			for i in range( int(n) ):
			
				for j in range( int(n) ):
					
					# Forward transform element	
					self.dft[i][j]  = ( 1.0 / n ) * np.power( W, i*j ) 

					# Inverse transform element				
					self.idft[j][i] = np.power( W,  i*j ).conj()

		
		# Calculate Discrete Fourier Transform
		def DFT(self, vt): 
			
			return np.dot(self.dft, vt)

		# Calculate Inverse Discrete Fourier Transform 
		def IDFT(self, vf): 
			
			return np.dot(self.idft, vf)
			
		# Return DFT matrix
		def get_dft(self):

			return self.dft

		# Return IDFT matrix
		def get_idft(self):

			return self.idft

		# Get number of harmonics
		def get_n(self): 

			return self.harmonics["n"]

		# Get harmonic by index (angular frequency)
		def get_omega(self, _): 
		   
			return self.harmonics["omega"][_]

		# Get harmonic by index (period)
		def get_tau(self, _): 
		   
			return self.harmonics["tau"][_]	

	# A class which forms the dual vector of a vector f using the DFT/IDFT method. 
	# User specifys whether the vector passed is in the time or frequency domain. 
	class Dual:

		def __init__(self, signal, transform, domain):

			# Initialize DFT object for frequency
			self.DFT = transform

			# If initializing a vector time domain
			if domain == 'time':
				
				self.time = signal
				self.freq = self.DFT.DFT(signal)

				# Cache period
				self.tau = self.DFT.harmonics["tau"]

			# If initializing a vector frequency domain
			if domain == 'freq':
				
				self.time = self.DFT.IDFT(signal)
				self.freq = signal

				# Cache period
				self.tau = self.DFT.harmonics["tau"]

		# Method to synthesize a signal for a given set of Fourier coefficients
		def sampling(self, npoints = 128):
			
			# Signl vector 
			signal = np.zeros(npoints, dtype = complex) 

			# Time domain
			period = np.linspace(0, 1.0 / self.DFT.harmonics['freq'], npoints)

			# Build waveform
			for i, c in enumerate( self.freq ):

				# Calculate contribution from harmonic
				harmonic = [ c * np.exp( np.complex(0, self.DFT.get_omega(i) * t) ) for t in period ]

				# Signal is a linear comibnation of harmonics
				signal = np.add( signal, harmonic )

			return period, signal