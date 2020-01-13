# ---------------------------------------------------------------------------------
#   minispice-> diodeTransient.py
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

from minispice.nonlinear import companionModels
from minispice.nonlinear import componentModels
from minispice.Converter import *

# Signal tools for pulse
from minispice.signalTools import signalTools
from matplotlib import pyplot as plt

# Example class for transient simulation of diode. A diode is modeled as 
# a nonlinear conductance and nonlinear capacitance in parallel. 
class diode_transient: 

	def __init__(self, signal):

		# The signal we will simulate
		self.signal = signal

		# Initialize component models
		self.diode  = componentModels.diode()

		# Initialize companion models
		self.diodeR = companionModels.nonlinearR( self.diode )
		self.diodeC = companionModels.nonlinearC( self.diode )
	
	def solve(self, source_impedacnce, conv = 1e-6):

		diode_waveform = []

		# Initilize transient solver with first signal value
		vm = self.signal['waveform'][0]

		# Perform iteration over companion model
		for source_voltage in self.signal['waveform']:

			# vn is determined by final step of previous iteration (Eq 8.20)
			# if n = 0 then vn is determined by the start value.
			vn = diode_waveform[-1] if len(diode_waveform) > 1 else vm

			# Perform Newton iteration at each timestep to solve diode voltage
			while True:
			
				num = ( self.diodeR.im(vm) + self.diodeC.im(vm, vn, self.signal['delta']) + nortonI( source_voltage, source_impedacnce ) ) 
				den = ( self.diodeR.gm(vm) + self.diodeC.gm(vm, vn, self.signal['delta']) + nortonG( source_impedacnce ) )

				_vm = num / den

				if np.abs( _vm - vm ) <= conv:  

					diode_waveform.append(_vm)

					break 

				else: 

					vm = _vm

		# Return solution waveform
		return diode_waveform

# Main program
if __name__ == "__main__":


	# Create pulse parameters
	config = {
		"amplitude" : 1.0,
		"period"	: 1e-8,
		"npoints"	: 1e3, 
		"duty"		: 0.25
	}


	# Generate square pulse using signalTools
	signal = signalTools().pulse(config)

	# Create analysis object
	analysis = diode_transient( signal )
	
	# Create some test resistances(1 Ohm to 1024 Ohm)
	source_impedacnce = [ 2.0**_power for _power in range(11) ]

	# Create axes object	
	fig = plt.figure()
	ax0 = fig.add_subplot(111)
	ax0.set_xlabel("Time (s)")
	ax0.set_ylabel("Voltage (V)")
	ax0.set_title("Diode Resistor Circuit : Transient Response")

	# Loop through all source impeances
	for _ in source_impedacnce:

		diode_waveform = analysis.solve(source_impedacnce = _, conv = 1e-9 )

		h0, = ax0.plot(signal['time'], diode_waveform, color="tab:orange")

	# Plot the input waveform
	h1, = ax0.plot(signal['time'], signal['waveform'], color="tab:blue")
	
	# Create a simple legend 
	ax0.legend([h1,h0],["Input waveform", "Diode response\n $1\Omega \leq Z_s \leq 1024\Omega$"])

	plt.show()
