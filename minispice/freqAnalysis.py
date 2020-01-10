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
	def __init__(self, ygroup, freqList):
		self.ygroup = ygroup
		self.freqList = freqList

	# Overload constructor via @classmethod	
	@classmethod
	def fromFile(cls, path, freqList):
	
		# Read all components into data		
		with open(path, 'r') as f:
			data = [line.split() for line in f]

		# Get size of ymatrix from file
		size = 0
		for j,lst in enumerate(data):
			
			# Case of Transitor in data file
			if len(lst) > 4:
				if size < int(lst[1]):
					size = int(lst[1])
				if size < int(lst[2]):
					size = int(lst[2])
				if size < int(lst[3]):
					size = int(lst[3])

			# Case of Passives in data file
			else:
				if size < int(lst[1]):
					size = int(lst[1])
				if size < int(lst[2]):
					size = int(lst[2])

		ygroup = []
		for i, freq in enumerate(freqList):

			tmp = nodeMatrix(size,freq)
			for j, lst in enumerate(data): 
				
				# Exclude header lines
				if re.match(r'#!',str(lst[0])) is not None: 
					pass 

				# Otherwise loop through devices
				else:

					# Case of passives
					if re.match(r'R|C|L',str(lst[0])) is not None:
						tmp.addPassive(str(lst[0]), int(lst[1]), int(lst[2]), float(lst[3]))
		
					# Case of transistor
					elif re.match(r'Q', str(lst[0])) is not None: 
						tmp.addTransistor(str(lst[0]), int(lst[1]), int(lst[2]), int(lst[3]), str(lst[4]))

					# Case of a VCCS (transconductance)			
					elif re.match(r'G',str(lst[0])) is not None:
						tmp.addVCCS(str(lst[0]), int(lst[1]), int(lst[2]), int(lst[3]), int(lst[4]), float(lst[5]) )

					else:
						pass

			ygroup.append(tmp)
			
		return cls(ygroup, freqList)

	def plotGain(self,n1,n2, arg=None):
		
		# Obtain nodal gain for each object in y-group
		gain, phase = [],[]
		for i,ymatrix in enumerate(self.ygroup):
			gain.append(np.abs(self.ygroup[i].voltageGain(n1,n2)))
			phase.append(np.angle(self.ygroup[i].voltageGain(n1,n2),deg=True))

		# Long ugly plotting methods
		if arg=="lin":
			titleStr="vm("+str(n2)+")/"+"vm("+str(n1)+")"
			plt.figure(1)
			plt.title(titleStr)
			plt.plot(self.freqList,gain)
			plt.xlabel("Frequency")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.ylabel("Voltage Gain")
			
			titleStr="vp("+str(n2)+")"
			plt.figure(2)
			plt.plot(self.freqList, phase)
			plt.title(titleStr)
			plt.xlabel("Frequency")
			plt.ylabel("Voltage Phase Difference")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.show()

		if arg=="log":
			titleStr="vm("+str(n2)+")/"+"vm("+str(n1)+")"
			plt.figure(1)
			plt.title(titleStr)
			plt.semilogx(self.freqList, gain)
			plt.xlabel("Frequency")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.ylabel("Voltage Gain")
			
			titleStr="vp("+str(n2)+")"
			plt.figure(2)
			plt.semilogx(self.freqList, phase)
			plt.title(titleStr)
			plt.xlabel("Frequency")
			plt.ylabel("Voltage Phase Difference")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.show()

		if arg=="dB":
			titleStr="vm("+str(n2)+")/"+"vm("+str(n1)+")"
			plt.figure(1)
			plt.title(titleStr)
			plt.semilogx(self.freqList, [todB(_) for _ in gain])
			plt.xlabel("Frequency")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.ylabel("Voltage Gain")
			
			titleStr="vp("+str(n2)+")"
			plt.figure(2)
			plt.semilogx(self.freqList, phase)
			plt.title(titleStr)
			plt.xlabel("Frequency")
			plt.ylabel("Voltage Phase Difference")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.show()


	def inputImpedance(self,n1,n2,rl,arg):
		zMag,zPhase = [],[]
		g  = lambda r : complex(1/r) 
		for i,ymatrix in enumerate(self.ygroup):
			twoport=ymatrix.toTwoport(n1,n2)
			tmp = (np.linalg.det(twoport)+twoport[0,0])/(twoport[1,1]+g(rl))
			zMag.append(np.abs(1/tmp))
			zPhase.append(np.angle((1/tmp),deg=True))
		
		# Long ugly plotting methods			
		if arg=="lin":
			titleStr="mag(Zin)"
			plt.figure(1)
			plt.title(titleStr)
			plt.plot(self.freqList, zMag)
			plt.xlabel("Frequency")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.ylabel("Input Impedance Magnitude")
			
			titleStr="phase(Zin)"
			plt.figure(2)
			plt.title(titleStr)
			plt.plot(self.freqList, zPhase)
			plt.xlabel("Frequency")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.ylabel("Input Impedance Phase")
			plt.show()

		if arg=="log":
			titleStr="mag(Zin)"
			plt.figure(1)
			plt.title(titleStr)
			plt.semilogx(self.freqList, zMag)
			plt.xlabel("Frequency")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.ylabel("Input Impedance Magnitude")
			
			titleStr="phase(Zin)"
			plt.figure(2)
			plt.title(titleStr)
			plt.semilogx(self.freqList,zPhase)
			plt.xlabel("Frequency")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.ylabel("Input Impedance Phase")
			plt.show()	

	def outputImpedance(self,n1,n2,rs,arg):
		zMag,zPhase = [],[]
		g  = lambda r : complex(1/r) 
		for i,ymatrix in enumerate(self.ygroup):
			twoport=ymatrix.toTwoport(n1,n2)
			tmp = (np.linalg.det(twoport)+twoport[1,1])/(twoport[0,0]+g(rs))
			zMag.append(np.abs(1/tmp))
			zPhase.append(np.angle((1/tmp),deg=True))
		   
		# Long ugly plotting methods
		if arg=="lin":
			titleStr="mag(Zout)"
			plt.figure(3)
			plt.title(titleStr)
			plt.plot(self.freqList,zMag)
			plt.xlabel("Frequency")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.ylabel("Output Impedance Magnitude")
			
			titleStr="phase(Zout)"
			plt.figure(4)
			plt.title(titleStr)
			plt.plot(self.freqList,zPhase)
			plt.xlabel("Frequency")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.ylabel("Output Impedance Phase")

		if arg=="log":
			titleStr="mag(Zout)"
			plt.figure(3)
			plt.title(titleStr)
			plt.semilogx(self.freqList,zMag)
			plt.xlabel("Frequency")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.ylabel("Output Impedance Magnitude")
			
			titleStr="phase(Zout)"
			plt.figure(4)
			plt.title(titleStr)
			plt.semilogx(self.freqList,zPhase)
			plt.xlabel("Frequency")
			plt.grid(True, which="both",ls="-", color='0.65')
			plt.ylabel("Output Impedance Phase")
			plt.show()
