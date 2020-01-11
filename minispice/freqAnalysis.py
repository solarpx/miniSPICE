# ---------------------------------------------------------------------------------
# 	minispice-> freqAnalysis.py
#	Copyright (C) 2019 Michael Winters
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

# Classes for array manipulation
import numpy as np
import collections
import math
import copy
import re

# For plotting
import matplotlib.pyplot as plt

# Imprt node matrix
from .nodeMatrix import nodeMatrix
from .Converter import *

# Class to construct y-matrix for a list of frequencies
class freqAnalysis: 
	
	# Method to initialize directly
	def __init__(self, data, freq, components):

		# Data is dict of admittance matrices
		self.data = data

		# List of frequencies to simulate
		self.freq = freq
		
		# Dictionary to hold components
		self.components = components
	
	# Overload constructor via @classmethod	
	@classmethod
	def fromFile(cls, path, freq):
	
		# Read all components into dict
		components = {}

		# While parsing, we will determine the size of admittace matrix
		size = 0

		with open(path, 'r') as f:
		
			for _line in f:

				# Split component data
				_comp = _line.split()

				# Extract nodes from list
				_nodes = [ int(n) for n in _comp[1:len(_comp)-1] ]

				# Augment size if necessary
				if max(_nodes) > size:
					size = max(_nodes)

				# Parse data into dictionary 
				try:
					components[ _comp[0] ] = {"nodes" : _nodes, "value" : float(_comp[-1]) }

				except:
					components[ _comp[0] ] = {"nodes" : _nodes, "value" : str(_comp[-1]) }

		# Create a dictionary for admittance matrices
		data = collections.OrderedDict()
				
		# Loop through frequencies and initialize components
		for f in freq:

			# Initialize node matrix for frequency f
			ymatrix = nodeMatrix(size, f)
		
			# Loop through components and add them in
			for _comp, _conf in components.items():

				# Passive components
				if re.match(r'R|C|L', _comp) is not None:
	
					ymatrix.addPassive(
						_comp, 
						_conf["nodes"][0], 
						_conf["nodes"][1], 
						_conf["value"]
					)

				# Transistor with model file
				elif re.match(r'Q', _comp) is not None: 
					
					ymatrix.addTransistor(
						_comp, 
						_conf["nodes"][0], 
						_conf["nodes"][1], 
						_conf["nodes"][2], 
						_conf["value"]
					)

				# Case of a VCCS (transconductance)			
				elif re.match(r'G', _comp) is not None:
					
					ymatrix.addVCCS(
						_comp, 
						_conf["nodes"][0], 
						_conf["nodes"][1], 
						_conf["nodes"][2], 
						_conf["nodes"][3], 
						_conf["value"]
					)

				# Pass 
				else:
					pass

			# Assign data 
			data[f] = ymatrix
			
		return cls(data, freq, components)

	# Method to return S-parameters for two nodes over all frequencies
	def Sparameters(self, n1, n2):

		sdata = collections.OrderedDict()

		for f, ymatrix in self.data.items(): 
		
			sdata[f] = ytos( ymatrix.toTwoport(n1,n2) )	

		return sdata	

	# Return a single matrix from simulation
	def getMatrix(self, freq ):
		return self.data[ freq ] if freq in self.data.keys() else None

	# Method to return entire data structure	
	def getData(self):
		return self.data

	# Calculate voltage gain between two nodes
	def calcVoltageGain(self, n1, n2):	
		return [ np.abs( self.data[f].voltageGain(n1,n2) ) for f, ymatrix in self.data.items() ]

	# Calculate voltage phase shift for each frequency
	def calcVoltagePhase(self, n1, n2):	
		return [ np.angle( self.data[f].voltageGain(n1,n2), deg = True ) for f, ymatrix in self.data.items() ]

	# Calculate gain of effective twoport network connected to Rs and Rl
	def calcNetworkGain(self, n1, n2, Rs, Rl):
		return [ np.abs( self.data[f].networkGain(n1,n2, Rs, Rl) ) for f, ymatrix in self.data.items() ]


# Helper class to plot frequency analysis data
class plotAnalysis: 

	def __init__(self, arg = "lin"):
		
		# All plots will have the following
		self.fig = plt.figure()

		# Create axes
		self.ax0 = self.fig.add_subplot(111)
		self.ax0.set_xlabel("Frequency (Hz)")
		self.ax0.grid(True, which="both",ls="-", color='0.65')

		# Plotting style
		self.arg = arg

	# Plotting method
	def plotGain(self, freq, gain, arg=None):
		
		# Linear frequency 
		if self.arg == "lin":		
			self.ax0.plot(freq, gain)
			self.ax0.set_ylabel("Gain")

		elif self.arg == "lin(dB)":
			self.ax0.plot(freq, [todB(_) for _ in gain])
			self.ax0.set_ylabel("Gain (dB)")
		
		# Log frequency 
		elif self.arg == "log":
			self.ax0.semilogx(freq, gain)
			self.ax0.set_ylabel("Gain")
			
		elif self.arg == "log(dB)":
			self.ax0.semilogx(freq, [todB(_) for _ in gain])
			self.ax0.set_ylabel("Gain (dB)")
	
	# Plotting method
	def plotPhase(self, freq, phase, arg=None):

		# Linear frequency 
		if self.arg == "lin":		
			self.ax0.plot(freq, gain)
			self.ax0.set_ylabel("Phase (deg)")

		# Log frequency 
		elif self.arg == "log":
			self.ax0.semilogx(freq, gain)
			self.ax0.set_ylabel("Phase (deg)")

	# Wrap mpl show plots
	def show(self):
		plt.show()

	


	#### IMPEDANCE CALCULATIONS DO NOT BELONG HERE ####
	## MOVE TO nodeMatrix.py ##	
	# # Method to plot input impedance for a given load impedance (rl = 50 Ohm)
	# def plotInputImpedance(self,n1,n2,rl,arg):
		
	# 	g  = lambda r : complex(1/r)
	# 	zMag, zPhase = [], []
		 
	# 	# Calculate input impedance
	# 	for f, ymatrix in self.data.items():
			
	# 		twoport = ymatrix.toTwoport(n1,n2)
	# 		delta = np.linalg.det(twoport)
	# 		tmp = ( twoport[1,1] + g(rl) )/( delta + (twoport[0,0]*g(rl)) )

	# 		# Append data to array	
	# 		zMag.append( np.abs(tmp) )
	# 		zPhase.append( np.angle(tmp, deg=True) )
		
	# 	# Long ugly plotting methods			
	# 	if arg=="lin":
	# 		titleStr="mag(Zin)"
	# 		plt.figure(1)
	# 		plt.title(titleStr)
	# 		plt.plot(self.freq, zMag)
	# 		plt.xlabel("Frequency")
	# 		plt.grid(True, which="both",ls="-", color='0.65')
	# 		plt.ylabel("|Zin|")
			
	# 		titleStr="phase(Zin)"
	# 		plt.figure(2)
	# 		plt.title(titleStr)
	# 		plt.plot(self.freq, zPhase)
	# 		plt.xlabel("Frequency")
	# 		plt.grid(True, which="both",ls="-", color='0.65')
	# 		plt.ylabel("<Zin (deg)")
	# 		plt.show()


	# 	if arg=="log":
	# 		titleStr="mag(Zin)"
	# 		plt.figure(1)
	# 		plt.title(titleStr)
	# 		plt.semilogx(self.freq, zMag)
	# 		plt.xlabel("Frequency")
	# 		plt.grid(True, which="both",ls="-", color='0.65')
	# 		plt.ylabel("|Zin|")
			
	# 		titleStr="phase(Zin)"
	# 		plt.figure(2)
	# 		plt.title(titleStr)
	# 		plt.semilogx(self.freq,zPhase)
	# 		plt.xlabel("Frequency")
	# 		plt.grid(True, which="both",ls="-", color='0.65')
	# 		plt.ylabel("<Zin (deg)")
	# 		plt.show()	


	# # Method to plot output impedance for a given source impedance (rs = 50 Ohm)
	# def plotOutputImpedance(self, n1, n2, rs, arg):

	# 	g  = lambda r : complex(1/r)
	# 	zMag, zPhase = [], []
		 
	# 	# Calculate output impedance
	# 	for f, ymatrix in self.data.items():
			
	# 		twoport = ymatrix.toTwoport(n1, n2)
	# 		delta = np.linalg.det(twoport)
	# 		tmp = ( twoport[0,0] + g(rs) )/( delta + (twoport[1,1]*g(rs)) ) 

	# 		# Append data to array	
	# 		zMag.append( np.abs( tmp ) )
	# 		zPhase.append( np.angle(tmp, deg=True) )
		   
	# 	# Long ugly plotting methods
	# 	if arg=="lin":
	# 		titleStr="mag(Zout)"
	# 		plt.figure(3)
	# 		plt.title(titleStr)
	# 		plt.plot(self.freq, zMag)
	# 		plt.xlabel("Frequency")
	# 		plt.grid(True, which="both",ls="-", color='0.65')
	# 		plt.ylabel("|Zout|")
			
	# 		titleStr="phase(Zout)"
	# 		plt.figure(4)
	# 		plt.title(titleStr)
	# 		plt.plot(self.freq, zPhase)
	# 		plt.xlabel("Frequency")
	# 		plt.grid(True, which="both",ls="-", color='0.65')
	# 		plt.ylabel("<Zout (deg)")

	# 	if arg=="log":
	# 		titleStr="mag(Zout)"
	# 		plt.figure(3)
	# 		plt.title(titleStr)
	# 		plt.semilogx(self.freq, zMag)
	# 		plt.xlabel("Frequency")
	# 		plt.grid(True, which="both",ls="-", color='0.65')
	# 		plt.ylabel("|Zout|")
			
	# 		titleStr="phase(Zout)"
	# 		plt.figure(4)
	# 		plt.title(titleStr)
	# 		plt.semilogx(self.freq, zPhase)
	# 		plt.xlabel("Frequency")
	# 		plt.grid(True, which="both",ls="-", color='0.65')
	# 		plt.ylabel("<Zout (deg)")
	# 		plt.show()
