# ---------------------------------------------------------------------------------
# 	minispice-> nodeMatrix.py
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
import os

# Node admittace matrix class
class nodeMatrix: 

	def __init__(self, size, freq): 
		self.ymatrix = np.zeros(shape=(size,size),dtype=complex)
		self.freq = float(freq)
		self.size = size

	def showMatrix(self):
		print(self.ymatrix)

	def addPassive(self, name, n1, n2, value):
		# Convert to float automatically
		value = float(value)

		# Resistances (g = 1/r)
		g = lambda r: complex(1/r) 
		if re.match(r'R\d*',name) is not None:
			if n1 == 0:  
				n2-=1
				self.ymatrix[n2,n2] += g(value)
			elif n2 == 0:  
				n1-=1
				self.ymatrix[n1,n1] += g(value)
			else:
				n1-=1
				n2-=1
				self.ymatrix[n1,n1] += g(value)
				self.ymatrix[n1,n2] -= g(value)
				self.ymatrix[n2,n1] -= g(value)
				self.ymatrix[n2,n2] += g(value)

		# Capacitances (bc = iwC)
		bc = lambda c : complex(0,2*math.pi*self.freq*c) 
		if re.match(r'C\d*',name) is not None:
			value = complex(value,0)
			if n1 == 0:  
				n2-=1
				self.ymatrix[n2,n2] += bc(value)
			elif n2 == 0:  
				n1-=1
				self.ymatrix[n1,n1] += bc(value)
			else:
				n1-=1
				n2-=1
				self.ymatrix[n1,n1] += bc(value)
				self.ymatrix[n1,n2] -= bc(value)
				self.ymatrix[n2,n1] -= bc(value)
				self.ymatrix[n2,n2] += bc(value)

		# Inductances (bl = iwL)
		bl = lambda l : complex(0,(-1/(2*math.pi*self.freq*l))) 
		if re.match(r'L\d*',name) is not None:
			value = complex(value,0)
			if n1 == 0:  
				n2-=1
				self.ymatrix[n2,n2] += bl(value)
			elif n2 == 0:  
				n1-=1
				self.ymatrix[n1,n1] += bl(value)
			else:
				n1-=1
				n2-=1
				self.ymatrix[n1,n1] += bl(value)
				self.ymatrix[n1,n2] -= bl(value)
				self.ymatrix[n2,n1] -= bl(value)
				self.ymatrix[n2,n2] += bl(value)

	
	# Transconductances
	def addVCCS(self,name,n1,n2,n3,n4,value):
	
		#3,0,2,0
		value = complex(value,0)

		if re.match(r'G\d*',name) is not None:
	
			# Case of a twoport
			if n2 == 0 and n4 == 0:

				n1-=1
				n3-=1
				self.ymatrix[n1,n3] += value
	
			# General case
			else:
				n1-=1
				n2-=1
				n3-=1
				n4-=1
				self.ymatrix[n2,n3]-=value
				self.ymatrix[n1,n3]+=value
				self.ymatrix[n1,n4]-=value
				self.ymatrix[n2,n4]+=value
	

	# Method for adding a transistor
	def addTransistor(self,name, nb, nc, ne, model): 
		
		nb-=1
		nc-=1
		ne-=1

		# Simple Model with only b and rbe
		g = lambda r : complex(1/r) 
		if model == 'simple':
			# Extract transistor parameters out of params dict			
			params = self.getModel('simple')
			b=float(params['b'])
			rbe=float(params['rbe'])
			# Base
			self.ymatrix[nb,nb]+=g(rbe)
			self.ymatrix[nb,nc]+=0
			self.ymatrix[nb,ne]-=g(rbe)
			# Collector
			self.ymatrix[nc,nb]+=b*g(rbe)
			self.ymatrix[nc,nc]+=0
			self.ymatrix[nc,ne]-=b*g(rbe)
			# Emitter
			self.ymatrix[ne,nb]-=(b+1)*g(rbe)
			self.ymatrix[ne,nc]+=0
			self.ymatrix[ne,ne]+=(b+1)*g(rbe)

		# Intrinsic transistor pi model
		g  = lambda r : complex(1/r) 
		bc = lambda c : complex(0,2*math.pi*self.freq*c) 
		if model == 'hybridpi':
			# Extract transistor parameters out of params dict			
			params = self.getModel('hybridpi')
			
			gm = float(params['gm'])
			r0 = float(params['rce'])
			rpi= float(params['rbe'])
			cpi= float(params['cbei'])
			cmu= float(params['cbc'])

			# Base
			self.ymatrix[nb,nb]+=(g(rpi)+bc(cpi)+bc(cmu))
			self.ymatrix[nb,nc]+=0
			self.ymatrix[nb,ne]-=(g(rpi)+bc(cpi)+bc(cmu)) 
			# Collector
			self.ymatrix[nc,nb]+=(gm-bc(cmu))
			self.ymatrix[nc,nc]+=g(r0)
			self.ymatrix[nc,ne]+=(bc(cmu)-g(r0)-gm)
			# Emitter
			self.ymatrix[ne,nb]-=((g(rpi)+bc(cpi)+gm))
			self.ymatrix[ne,nc]-=g(r0)
			self.ymatrix[ne,ne]+=(g(rpi)+bc(cpi)+g(r0)+gm)


		#Transistor with base spreading resistance
		g  = lambda r : complex(1/r) 
		bc = lambda c : complex(0,2*math.pi*self.freq*c)
		if model == 'hybridpix':
			# Extract transistor parameters out of params dict			
			params = self.getModel('hybridpix')
			
			gm = float(params['gm'])
			rce = float(params['rce'])
			rbe= float(params['rbe'])
			cbe= float(params['cbe'])
			cbc= float(params['cbc'])
			rbb = float(params['rbb'])

			# Construct the CE admittance parameters
			y11 = (g(rbe)+bc(cbe)+bc(cbc))
			y12 = 0
			y21 = (gm-bc(cbc))
			y22 = g(rce)
			rbbDce = (rbb*((y11*y22)-(y21*y12)))
			s = (1/(1+y11*rbb))

			# Base
			self.ymatrix[nb,nb]+=((y11)*s)
			self.ymatrix[nb,nc]+=((y12)*s)
			self.ymatrix[nb,ne]-=((y11+y12)*s)
			# Collector
			self.ymatrix[nc,nb]+=((y21)*s)
			self.ymatrix[nc,nc]+=((y22+rbbDce)*s)
			self.ymatrix[nc,ne]-=((y21+y22+rbbDce)*s)
			# Emitter
			self.ymatrix[ne,nb]-=((y11+y21)*s)
			self.ymatrix[ne,nc]-=((y12+y22+rbbDce)*s)
			self.ymatrix[ne,ne]+=((y11+y22+y12+y21+rbbDce)*s)


	# Method to extract params from a *.model file 
	def getModel(self,name):
		params = {}
		path = os.getcwd()+ os.path.sep + name + '.model'
		
		with open(path, 'r') as f:
			data = [line.split() for line in f]
		
		for i,lst in enumerate(data):
			params[str(lst[0])] = float(lst[1])
		
		return params

	# Method to calculate cofactors Dij
	def cofactorN(self,i,j): 
		try:
			if i<1 or j<1:
				print("Invalid Cofactor Index")
				raise ValueError

		except ValueError:
			return None

		toPower=i+j
		i-=1
		j-=1

		# Store copy of Y matrix and remove rows i and j
		# Need deep copy so operations on A dont change ymatrix
		A = copy.deepcopy(self.ymatrix)
		A[i,:] = np.NaN
		A[:,j] = np.NaN
		
		Am = np.ma.masked_invalid(A)
		Am = np.ma.compressed(Am)
		Am = np.array(self.chunks(Am,self.size-1))
		return np.linalg.det(Am)*((-1)**toPower)
	
	# Method to calculate cofactor Dii,jj
	def cofactorD(self,i,j): 
		
		#Check bounds of cofactor index
		try:
			if i<1 or j<1:
				print("Invalid Cofactor Index")
				raise ValueError
		except ValueError:
			return None

		#In case of twoport return 1
		try: 
			if self.size == 2:
				raise ValueError
		except ValueError:
			return 1

		i-=1
		j-=1
		# Store copy of Y matrix and delete rows using NaN mask
		A = copy.deepcopy(self.ymatrix)
		A[i,:] = np.NaN
		A[j,:] = np.NaN
		A[:,i] = np.NaN
		A[:,j] = np.NaN

		# Compress out the remaning data
		Am = np.ma.masked_invalid(A)
		Am = np.ma.compressed(Am)
		Am = np.array(self.chunks(Am,self.size-2))
		
		# Am is now the submatrix we desire generate determinant
		return np.linalg.det(Am)

	# Method which calculates cofactors and returns the corresponding twoport parameters
	def toTwoport(self, n1, n2):

		# Initialize 2x2 matrix of zeros
		twoport = np.zeros(shape=(2,2),dtype=complex)

		# Calculate Cofactors and generate Y_twoport (p.122)
		twoport[1,1] += self.cofactorN(n1,n1)
		twoport[1,0] += self.cofactorN(n1,n2)*(-1)
		twoport[0,1] += self.cofactorN(n2,n1)*(-1)
		twoport[0,0] += self.cofactorN(n2,n2)
		Delta = self.cofactorD(n1,n2)
		twoport=twoport/Delta
		return twoport

	# Calculate node gain
	def voltageGain(self, n1, n2):

		#The voltage gain is given by the ratio of cofactors
		return self.cofactorN(n1,n2) / self.cofactorN(n1,n1)

	# Calculate gain in a network
	def networkGain(self, n1, n2, Zs = 50., Zl = 50.):

		# Compress admittance matrix to twoport
		tp = self.toTwoport(n1,n2)

		# Source and load admittances
		Ys = complex(1./Zs)
		Yl = complex(1./Zl)

		# calculate twoport gain between nodes
		num = 4.0 * Ys.real * Yl.real * np.abs(tp[1,0])**2 
		den = np.abs( ( tp[0,0] + Ys ) * ( tp[1,1] + Yl ) - tp[0,1] * tp[1,0] )**2

		return num / den

	# Calculate twoport input impedace given a certain load 
	def inputImpedance(self, n1, n2, Zl = 50.):

		# Compress admittance matrix to twoport
		tp = self.toTwoport(n1,n2)
		
		# Calculate twoport and determinant
		det = np.linalg.det(tp)

		# Load admittance
		Yl = complex(1./Zl)

		# Calculate input impedance
		num = tp[1,1] + Yl
		den = det + tp[0,0] * Yl

		return num/den

	# Calculate twoport output impedace given a certain source
	def outputImpedance(self, n1, n2, Zs = 50.):

		# Compress admittance matrix to twoport
		tp = self.toTwoport(n1,n2)
		
		# Calculate twoport and determinant
		det = np.linalg.det(tp)

		# Load admittance
		Ys = complex(1./Zs)

		# Calculate input impedance
		num = tp[0,0] + Ys
		den = det + tp[1,1] * Ys

		return num/den

	# Calculate transfer function
	def transferFunction(self, n1, n2):

		# Compress admittance matrix to twoport
		tp = self.toTwoport(n1,n2)

		return tp[1,0]/tp[1,1]

	# Breaks a list (l) into chunks of length (n)	
	def chunks(self, l, n):

		return [list(l[i:i+n]) for i in range(0, len(l), n)]
