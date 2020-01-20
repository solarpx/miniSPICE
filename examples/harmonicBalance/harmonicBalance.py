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
import minispice.discreteFourierTransform as DFT

# This class implements the harmonic balance method for solving the large 
# signal response for circuits containing nonlinear elements connected to 
# an arbitrary source impedance.
class harmonicBalance:

	# Initialization vector, frequency, convergence
	def __init__( self, config, nonlinear ):

		# Store configuration variables
		self.ampl = config["amplitude"] 
		self.freq = config["frequency"]
		self.order = config["order"]

		# Solver configuration variables
		self.epsilon = config["epsilon"]
		self.maxiter = config["maxiter"]
		self.converge = config["converge"]

		# Calculate transform matrices
		self.Transform = DFT.Transform(self.freq, self.order)

		# Initialize nonlinear component model(s)
		self.nonlinear = nonlinear

		# Generate source and initial conditions
		self.source = self.source()

	# Source dual (freq domain): [0.0, 1.0i, 0.0, 0.0 ... ]
	def source(self):
		
		_source = [ complex(0, -1.0 * self.ampl) if _ == 1 else complex(0, 0) for _ in range(self.order) ]

		return DFT.Dual( _source, self.Transform, domain="freq")

	# Signal dual (freq domain): [0.0, 0.0, 0.0, 0.0 ... ]
	def signal(self):
		
		_signal = [ complex(0.0) for _ in range(self.order) ]
		
		return DFT.Dual( _signal, self.Transform, domain="freq")

	# Harmonic balance solver
	def solve( self, source_impedance ):	

		# Linear admittance matrices
		Y11 = self.Transform.zeros()
		Y12 = self.Transform.zeros()

		# Diagonal matrix of harmonics
		iOMEGA = self.Transform.zeros()

		# Define source admittance
		G = 1./source_impedance

		for i in range( self.order ):
		
			for j in range( self.order ):
		
				Y11[i][j] = complex( 1.0 * G) if i == j else complex(0.0)	   

				Y12[i][j] = complex(-1.0 * G) if i == j else complex(0.0)	   

				iOMEGA[i][j] = complex( 0.0, self.Transform.get_omega(i) ) if i == j else complex(0.0) 

		# Everything after here iterates the nonlinearity
		self.signal = self.signal()


		# Convergence comparison
		self.step = 0

		# Iteration loop
		while True:
			
			# 1) Calculate the nonlinear signal (time domain)			
			_Fv = [ self.nonlinear.f(v) for v in self.signal.time ]

			Fv = DFT.Dual(_Fv, self.Transform, "time")

			# 2) Check if harmonics are balanced (frequency domain)
			_ef  = np.dot(Y12, self.source.freq) + np.dot(Y11, self.signal.freq) + Fv.freq
			
			ef = DFT.Dual(_ef, self.Transform, "freq")

			# Calculate convergence criteria 
			delta = ( sum( np.abs(ef.time.real) ) / self.order )

			# Print convergence condition
			if ( self.step % 20 ) == 0:
	
				print("Conv: %s"%delta)

			# Stop on max_iterations
			if ( self.step ) > self.maxiter:

				print( "Maximum number of iterations %s", self.step)

				return self.source, self.signal				

			# If convergence criteria is met stop iteration
			if ( delta ) < self.converge:

				return self.source, self.signal

			# Otherwise calculate the Jacobian and iterate 
			else:   

				# Conductance and capacitance terms
				dG = self.Transform.zeros()
				dC = self.Transform.zeros()

				for i in range( self.order ):

					for j in range( self.order ):

						# Resistive term in Jacobian (time domain)
						dG[i][j] = complex( self.nonlinear.df( self.signal.time[i] ) ) if i == j else complex(0.0)

						# Capacitive term in Jacobian (time domain)
						dC[i][j] = complex( self.nonlinear.c( self.signal.time[i] ) ) if i == j else complex(0.0)


				# Jacobian: Resistive subterm
				JR = np.dot( self.Transform.get_dft(), np.dot(dG, self.Transform.get_idft()) ) # [F][G][F^-1]

				# Jacobian: Capacitive subterm
				JC = np.dot( self.Transform.get_dft(), np.dot(dC, self.Transform.get_idft()) ) # [F][C][F^-1]

				# Jacobian: Total
				J = Y11 + JR + np.dot(iOMEGA, JC)

				# Invert this matrix and multiply it by the error vector
				# to obtain change in voltage (frequency domain).
				_dV = np.dot( la.inv(J), ef.freq )

				## Add to original vector V and calculate its dual.
				self.signal = DFT.Dual(self.signal.freq - self.epsilon * _dV, self.Transform, "freq")
			
				# Increment iteration step counter
				self.step += 1

# Main program
if __name__ == "__main__":
	
	# Initialize harmonic balance object
	config = {
		"amplitude" : 1.2,
		"frequency" : 5.0e9,
		"order"		: 64,
		"epsilon"	: 0.001,
		"converge"	: 1e-9,
		"maxiter"	: 8192
	}

	# Run the simulation
	HB = harmonicBalance(config, componentModels.diode() )

	Vsource, Vsignal = HB.solve(source_impedance = 10.0)

	# Figure 1: Waveform
	fig = plt.figure(1)
	ax0 = plt.subplot(111)
	ax0.set_title("%s - Resistor Circuit (Waveform)"%HB.nonlinear.name )
	ax0.set_xlabel("Time (s)")
	ax0.set_ylabel("Voltage (V)")

	h0, = ax0.plot( np.real(Vsource.tau), np.real(Vsource.time) )
	h1, = ax0.plot( np.real(Vsignal.tau), np.real(Vsignal.time) )

	ax0.legend([h0, h1],["Input waveform", "Diode response"])

	# Figure 2: Harmonic composition
	fig = plt.figure(2)
	ax0 = plt.subplot(111)
	ax0.set_title("%s - Resistor Circuit (Harmonics)"%HB.nonlinear.name )
	ax0.set_xlabel("Harmonic (n)")
	ax0.set_ylabel("Coefficient Amplitude |c|")

	h0, = ax0.semilogy( np.abs(Vsignal.freq), 'o' )

	# Figure 3: Nonlinear conductance and capacitance
	fig = plt.figure(3)
	ax0 = plt.subplot(111)
	ax0.set_title("%s - Resistor Circuit (Nonlinear Elements)"%HB.nonlinear.name )
	ax0.set_xlabel("Time $(s)$")
	ax0.set_ylabel("Conductance $(G)$")
	ax1 = ax0.twinx()
	ax1.set_ylabel("Capacitance $(C)$")
	
	Gnl = [ HB.nonlinear.df(_) for _ in np.real(Vsignal.time) ]
	Cnl = [ HB.nonlinear.c(_)  for _ in np.real(Vsignal.time) ]

	h0, = ax0.semilogy(np.real(Vsource.tau), Gnl , color = "tab:blue" )
	h1, = ax1.plot(np.real(Vsource.tau), Cnl , color = "tab:orange" )

	ax0.legend([h0, h1],["Conductance", "Capacitance"])

	plt.show()
