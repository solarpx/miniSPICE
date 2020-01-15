# ---------------------------------------------------------------------------------
#   minispice -> nonlinear/componentModels.py
#   Copyright (C) 2020 Michael Winters
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
import numpy as np

# Class to represent a diode
class diode:

	def __init__(self): 

		self.name = "Diode"

		self.G   = 0.33333		# Diode conductance
		self.n   = 1.3			# Ideality factor
		self.Vt  = 0.026		# Thermal voltage (kbT@300K")
		self.I0  = 1e-11		# Reverse bias current
		self.Cj	 = 1.2e-12		# Zero bias junction capacitance
		self.phi = 1.3 			# Diode built in potential
		self.g   = 0.5 			# Capicitance factor

	# Diode current equation
	def f(self, _v):
	
		return ( self.I0 * np.exp( _v / ( self.n * self.Vt) ) - self.I0 )

	# Diode current derivative
	def df(self, _v): 

		return ( self.I0 / ( self.n * self.Vt) ) * np.exp( _v / ( self.n * self.Vt) ) 

	# Diode capacitance equation
	def c(self, _v):

		return self.Cj * np.power( ( 1 - ( _v / self.phi) ), -1.0 * self.g )	

	# Diode capacitance derivative
	def dc(self, _v): 

		return ( self.Cj * self.g / self.phi ) * np.power( ( 1 - ( _v / self.phi) ), -1.0 * (self.g  + 1.0 ) )
		
# Class to represent a cross diode
class cross_diode:

	def __init__(self): 

		self.name = "Cross Diode"

		self.diode = diode()

	# Diode current equation
	def f(self, _v):
	
		return self.diode.f(_v) + self.diode.f(-_v)

	# Diode current derivative
	def df(self, _v): 

		return self.diode.df(_v) + self.diode.df(-_v)

	# Diode capacitance equation
	def c(self, _v):

		return self.diode.c(_v) + self.diode.c(-_v)

	# Diode capacitance derivative
	def dc(self, _v): 

		return self.diode.dc(_v) + self.diode.dc(-_v)


# Class to represent a cross diode
class vdp_conductance:

	def __init__(self): 

		self.name = "VDP Conductance"

	# Diode current equation
	def f(self, _v):

		return np.power(_v, 3.0) / 3.0 - _v

	def df(self, _v):

		return np.power(_v, 2.0) - 1

	def c(self, _v):

		return 0.0

	def dc(self, _v):

		return 0.0

# Nonlinear transistor model
class HFET:

    def __init__(self):

        self.a = -0.080
        self.b = 3.3
        self.c = 0.41
        self.d = 10.0
        self.g = 0.056
        self.phi = 2.4 * ( np.pi / 180.)

    ## f1(ugd+ , ugs+)*f2(vgs-vgd) - f1(ugs- , ugd-)*f2(vgd-vgs)
    def f(self, _Vgd, _Vgs):
        
        # Source drain current difference
        delta = _Vgs - _Vgd

        # (ugd+, ugs+)
        _pgd, _pgs = self._plus(_Vgd, _Vgs) 

        # (ugd-, ugs)    
        _mgd, _mgs = self._minus(_Vgd, _Vgs) 
    
        # Calculate subcurrents
        _fp = self.f1(_pgd, _pgs) * self.f2(  1.0 * delta )
        _fm = self.f1(_mgs, _mgd) * self.f2( -1.0 * delta )

        # Return current
        return  self.g * ( _fp - _fm )

    # Calculate (ugd+, ugs+) terms
    def _plus(self, _Vgd, _Vgs): 

        _pgs =  np.cos( self.phi ) * _Vgs - np.sin( self.phi ) * _Vgd 
        _pgd =  np.cos( self.phi ) * _Vgd + np.sin( self.phi ) * _Vgs 
        
        return _pgd, _pgs

    # Calculate (ugd-, ugs-) terms
    def _minus(self, _Vgd, _Vgs): 

        _mgs =  np.sin( self.phi ) * _Vgd + np.cos( self.phi ) * _Vgs 
        _mgd =  np.cos( self.phi ) * _Vgd - np.sin( self.phi ) * _Vgs 

        return _mgd, _mgs

    def f1(self, _u1, _u2 ):

        const = np.exp( -self.b * ( _u2 + self.c) )

        return (1.0 + self.a * _u1)*(1.0 - np.tanh(const))

    def f2(self, _v ):

        const = np.exp(-self.d * _v)

        return (1.0 - np.tanh(const))