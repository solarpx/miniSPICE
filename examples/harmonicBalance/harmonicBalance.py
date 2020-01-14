# ---------------------------------------------------------------------------------
# 	minispice -> examples/harmonicBalance.py
#	Copyright (C) 2020 Michael Winters
#	github: https://github.com/mesoic
#	email:  mesoic@protonmail.com
# ---------------------------------------------------------------------------------
#	
#	Permission is hereby granted, free of charge, to any person obtaining a copy
#	of this software and associated documentation files (the "Software"), to deal
#	in the Software without restriction, including without limitation the rights
#	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#	copies of the Software, and to permit persons to whom the Software is
#	furnished to do so, subject to the following conditions:
#	
#	The above copyright notice and this permission notice shall be included in all
#	copies or substantial portions of the Software.
#	
#	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#	SOFTWARE.
#

#!/usr/bin/env python 
from matplotlib import pyplot as plt
import numpy.linalg as la
import numpy as np

# Import diode models
from minispice.nonlinear import componentModels
from minispice.Converter import *

# Import DFT class
from minispice.signalTools import DiscreteFourierTransform 

# Class to simulate diode resistor circuit using the harmonic balance method
class harmonicBalance:

	# Initialization vector, frequency, convergence
	def __init__(self, amplitude = 1.0, freq = 1e9, order = 20):

		# Initialize DFT object for frequency
		self.DFT = DiscreteFourierTransform( freq, order )

		# Set harmonic order
		self.order = order
		self.freq  = freq

		# Initialize component models
		self.diode  = componentModels.diode()

		# Generate source and initial conditions
		self.generate_source(amplitude)

	def generate_source(self, amplitude):
		
		# Set up a sin wave voltage source in freq domian [0.0, 1.0, 0.0, 0.0 ... ]
		_source = [ complex(0, -1.0*amplitude) if _ == 1 else complex(0, 0) for _ in range(self.order) ]

		self.source = self.DFT.dual( _source, self.freq, domain="freq")

	def initialize(self):
		
		# Create test vector to initialize the harmonic balance
		_V = [ complex(0.0) for _ in range(self.order) ]
		
		return self.DFT.dual(_V , self.freq, domain="freq")

	def solve(self, epsilon = 0.01, conv = 1e-4, max_iterations = 1e3):	

		# Everything up to here is the linear network

		# Linear admittance matrices
		Y11 = self.DFT.zeros()
		Y12 = self.DFT.zeros()

		# Diagonal matrix of harmonics
		OMEGA = self.DFT.zeros()

		# Define source admittance
		G = 1./50.

		for i in range( self.order ):
		
			for j in range( self.order ):
		
				Y11[i][j] = complex( 1.0*G) if i == j else complex(0.0)	   

				Y12[i][j] = complex(-1.0*G) if i == j else complex(0.0)	   

				OMEGA[i][j] = complex(0.0,i * 2 * np.pi * self.freq) if i == j else complex(0.0) 


		# Everything after here iterates the nonlinearity
		self.V = self.initialize()

		# Convergence comparison
		self.prev = self.source
		self.step = 0

		# Iteration loop
		while True:
			
			# Calculate the nonlinear signal (time domain)			
			_Fv = [ self.diode.f(v) for v in self.V.time ]

			Fv = self.DFT.dual(_Fv, self.freq, "time")

			# Calculate the error function (frequency domain)
			_ef  = np.dot(Y11, self.V.freq) + np.dot(Y12, self.source.freq) + Fv.freq
			
			ef = self.DFT.dual(_ef, self.freq, "freq")

			# Calculate delta
			# delta = np.sum( ( self.prev.time.real - self.V.time.real )**2 ) / self.order
			delta = sum( ef.time.real**2 ) / self.order
			
			# Print convergence condition
			if ( self.step % 20 ) == 0:
	
				print("Conv: %s"%delta)

			# Stop on max_iterations
			if ( self.step ) > max_iterations:

				print( "Maximum number of iterations %s", self.step)

				return self.source, self.V				

			# If convergence criteria is met stop iteration
			if sum( ef.time.real**2 ) / self.order < conv:

				return self.source, self.V

			# Otherwise calculate the Jacobian and iterate 
			else:   

				# Conductance and capacitance terms
				dG = self.DFT.zeros()
				dC = self.DFT.zeros()

				for i in range( self.order ):

					for j in range( self.order ):

						## Resistive term in Jacobian
						dG[i][j] = complex( self.diode.df( self.V.time[i] ) ) if i == j else complex(0.0)

						## Capacitive term in Jacobian
						dC[i][j] = complex( self.diode.c( self.V.time[i] ) ) if i == j else complex(0.0)


				# Jacobian: Resistive subterm
				JR = np.dot( self.DFT.dft, np.dot(dG, self.DFT.idft) ) # [F][G][F^-1]

				# Jacobian: Capacitive subterm
				JC = np.dot( self.DFT.dft, np.dot(dC, self.DFT.idft) ) # [F][C][F^-1]

				# Jacobian: Total
				J = Y11 + JR + np.dot(OMEGA, JC)

				# Invert this matrix and multiply it by the error vector
				# to obtain change in voltage (frequency domain).
				_dV = np.dot( la.inv(J), ef.freq )

				# Cache previous value of V
				self.prev  = self.V

				## Add to original vector V and calculate its dual.
				self.V = self.DFT.dual(self.V.freq - epsilon * _dV, self.freq, "freq")
			
				# Increment iteration step counter
				self.step += 1

if __name__ == "__main__":
	
	# Initialize harmonic balance object
	HB = harmonicBalance( amplitude = 1.0, freq = 1e8, order = 41 )

	# Run the simulation
	Vsource, Vdiode = HB.solve( epsilon = 0.01, conv = 1e-5, max_iterations = 2048 )

	# Plot Results
	fig = plt.figure(1)
	ax0 = plt.subplot(111)
	ax0.set_title("Diode Resistor Circuit: Waveform")
	ax0.set_xlabel("Time Step (n)")
	ax0.set_ylabel("Voltage (V)")
	h0, = ax0.plot( Vsource.time.real )
	h1, = ax0.plot( Vdiode.time.real )
	ax0.legend([h0, h1],["Input waveform", "Diode response"])

	fig = plt.figure(2)
	ax0 = plt.subplot(111)
	ax0.set_title("Diode Resistor Circuit: Nonlinear Elements")
	ax0.set_xlabel("Time Step $(n)$")
	ax0.set_ylabel("Diode Conductance Waveform $G$")
	h0, = ax0.semilogy( [ HB.diode.df(_) for _ in Vdiode.time.real], color = "tab:blue" )
	ax1 = ax0.twinx()
	ax1.set_ylabel("Diode Capacitance Waveform $C$")
	h1, = ax1.plot( [ HB.diode.c(_) for _ in Vdiode.time.real], color = "tab:orange" )

	ax0.legend([h0, h1],["Conductance", "Capacitance"])

	plt.show()

