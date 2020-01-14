# ---------------------------------------------------------------------------------
# 	minispice -> activeAnalysis.py
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
from .Converter import *
import numpy as np

# This is a container class which contained important equations for calculating
# stability and gain data from transistor S-parameters. This is used primarily
# for designing microwave amplifiers.
class amplAnalysis:

	# Stability is calculated 
	def __init__(self, sparams):

		# Extract s-parameters
		self.sparams = sparams

		# Alias individual s-parameters
		self.s11 = self.sparams[0][0]
		self.s12 = self.sparams[0][1]
		self.s21 = self.sparams[1][0]
		self.s22 = self.sparams[1][1]

		# Alias D and K by calling them
		self.D = self.D()
		self.K = self.K()

		# Caluclate stability circles
		self.isc = self.inputStabilityCircle()
		self.osc = self.outputStabilityCircle()

	# Analysis of S-parameters	
	def analysis(self):	

		# Sparameters as phasors
		print("\nS-Parameters")
		print("\t|s11| = %f : <s11 = %f"%phasor(self.s11, "deg") )
		print("\t|s12| = %f : <s12 = %f"%phasor(self.s12, "deg") )
		print("\t|s21| = %f : <s21 = %f"%phasor(self.s21, "deg") )
		print("\t|s22| = %f : <s22 = %f"%phasor(self.s22, "deg") )
	
		# Stability parameters
		print("\nStability Parameters")
		p = phasor(self.D, "deg")
		print("\t|D| = %f : <D = %f : K = %f"%( p[0],  p[1], self.K ) )

		# Calculate maximum gain in case of stable device
		if ( self.K > 1.0 ):
			
			Gain = self.maxTransducerGain()
			print("\nMaximum Transducer Gain (Device Stable)")
			print("\tGmx = %f = %f dB : gmx = %f"%(Gain["Gmx"], Gain["Gmx_dB"], Gain["gmx"]) )

		# Calculate maximum gain in case of unstable device
		else:	

			Gain = self.maxAvailableGain()
			print("\nMaximum Available Gain (Device Unstable)")
			print("\tGmx = %f = %f dB : gmx = %f"%(Gain["Gmx"], Gain["Gmx_dB"], Gain["gmx"]) )

		# Stability circle data
		print("\nInput Stability Circle")
		p = phasor(self.isc['c'], "deg")
		print("\t|c| = %f : <c = %f : r = %f"%( p[0],  p[1], self.isc['r'] ) ) 

		print("\nOutput Stability Circle")
		p = phasor(self.osc['c'], "deg")
		print("\t|c| = %f : <c = %f : r = %f"%( p[0],  p[1], self.osc['r'] ) ) 
		
		# stdout buffer
		print("\n")

	## Stability Considerations K and Delta
	def D(self):
		return ( self.s11 * self.s22 ) - ( self.s21 * self.s12 )

	def K(self):
		Ka = ( 1 - (abs(self.s11)**2 ) - ( abs(self.s22)**2) + ( abs(self.D)**2) )
		Kb = ( 2*abs(self.s12 * self.s21) )
		return Ka/Kb

	# Calculates points on circle for given center and radius 
	def Circle(self, _c, _r, _npts=512):
		return [ _r*np.exp( complex(0,angle) ) + _c for angle in np.linspace(0, 2*math.pi, _npts) ]

	# Input stability circle 
	def inputStabilityCircle(self):

		# Center
		cs = (self.s11 - ( self.D * self.s22.conj())).conj() / ( abs(self.s11)**2 - abs(self.D)**2 )

		# Radius
		rs = abs( self.s12 * self.s21 ) / ( abs(self.s11)**2 - abs(self.D)**2 )

		# Return data
		return {"c" : cs, "r" : rs, "data" : self.Circle(cs, rs)}

	# Output stability circle
	def outputStabilityCircle(self):
		
		# Center
		cl = (self.s22 - (self.D * self.s11.conj())).conj() / ( abs(self.s22)**2 - abs(self.D)**2 )

		# Radius
		rl = abs( self.s12 * self.s21 ) / ( abs(self.s22)**2 - abs(self.D)**2 )

		# Return data
		return {"c" : cl, "r" : rl, "data" : self.Circle(cl, rl)}

	# Constant gain circles for specified gain in dB
	def constantGainCircle(self, Gp_dB):

		# Calculate normalized scalar gain
		Gp = fromDb(Gp_dB)
		gp = Gp/( abs(self.s21)**2 )

		# Center
		cpa = gp * ( self.s22 - self.D*self.s11.conj() ).conj()
		cpb = 1 + gp * ( abs(self.s22)**2 - abs(self.D)**2 )
		cg  = cpa/cpb
	    
		# Radius
		rpa = math.sqrt( 1 - (2 * self.K * gp * abs(self.s12 * self.s21) ) + (gp * abs(self.s12 * self.s21))**2 )
		rpb = abs(1 + gp * (abs(self.s22)**2 - abs(self.D)**2) )
		rg  = rpa/rpb

		# Return data
		return {"c" : cg, "r" : rg, "data" : self.Circle(cg, rg)}

	# Maximum transducer gain
	def maxTransducerGain(self):

		# Maximum transducer gain
		Gmx = ( abs(self.s21)/abs(self.s12) ) * ( self.K - math.sqrt(self.K**2 -1) )
		Gmx_dB = todB(Gmx)

		# Normalized transducer gain
		gmx = Gmx / abs(self.s21 * self.s21.conj())
		
		return {"Gmx" : Gmx, "Gmx_dB" : Gmx_dB, "gmx" : gmx}

	# Maximum available gain
	def maxAvailableGain(self):

		# Maximum available gain
		Gmx = ( abs(self.s21)/abs(self.s12) )
		Gmx_dB = todB(Gmx)			

		# Normalized maximum available gain
		gmx = Gmx / abs(self.s21 * self.s21.conj())

		return {"Gmx" : Gmx, "Gmx_dB" : Gmx_dB, "gmx" : gmx}

	# Method to calculate source conjugate 
	def conjugateCircleData(self, data):
		return [ ( self.s11 + (self.s12 * self.s21 * _)/(1- (self.s22 * _)) ).conj() for _ in data ]
