


from .Converter import *
import numpy as np


class activeAnalysis:

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
			
			# Maximum transducer gain
			Gmx = ( abs(self.s21)/abs(self.s12) ) * ( self.K - math.sqrt(self.K**2 -1) )
			Gmxdb = todB(Gmx)

			# Normalized transducer gain
			gmx = Gmx / abs(self.s21 * self.s21.conj())
		
			print("\nMaximum Gain (Device Stable)")
			print("\tGmx = %f = %f dB : gmx = %f"%(Gmx, Gmxdb, gmx))

		# Calculate maximum gain in case of unstable device
		else:	

			# Maximum available gain
			Gmx = ( abs(self.s21)/abs(self.s12) )
			Gmxdb = todB(Gmx)			

			# Normalized maximum available gain
			gmx = Gmx / abs(self.s21 * self.s21.conj())
			
			print("\nMaximum Gain (Device Unstable)")
			print("\tGmx = %f = %f dB : gmx = %f"%(Gmx, Gmxdb, gmx))


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

	# Gain circles
	def GainCircle(self, gp):

		# Center
		cpa = gp*(s22-D*s11.conj()).conj()
		cpb = 1+gp*(abs(s22)**2 - abs(D)**2)
		cg  = _cpa/_cpb
	    
		# Radius
		rpa = math.sqrt(1 - (2*K*gp*abs(s12*s21)) + (gp*abs(s12*s21))**2)
		rpb = abs(1 + gp*(abs(s22)**2 - abs(D)**2))
		rg  = rpa/rpb

		# Return data
		return {"c" : cg, "r" : rg, "data" : self.CSircle(cl, rl)}