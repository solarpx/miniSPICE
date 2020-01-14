# ---------------------------------------------------------------------------------
#   minispice -> nonlinear/companionModels.py
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

# General nonlinear resistor companion model. 
# This is equivalent to Newtown's method
class nonlinearR:

	# This class is initialized with models. 
	# Models must define f(v) and its derivative.
	def __init__(self, _model):

		self.model = _model

	# Companion model: linear conductance.
	def gm(self, _vm):	

		return self.model.df( _vm ) 

	# Companion model: constant current source
	def im(self, _vm):

		return ( self.model.df( _vm ) * _vm  - self.model.f( _vm ) )


# General nonlinear capacitance companion model. 
# This is equivalent to Newtown's method
class nonlinearC:

	# This class is initialized with models. 
	# Models must define c(v) and its derivative.
	def __init__(self, _model):

		self.model = _model

	# Companion model: linear conductance.
	def gm(self, _vm, _vn, T):	

		return ( self.model.c(_vm) + self.model.dc(_vm) * (_vm - _vn) ) / T 

	# Companion model: constant current source.
	def im(self, _vm, _vn, T):

		return ( self.model.c(_vm) * _vn  + self.model.dc(_vm) * (_vm - _vn) * _vm ) / T
		