# ---------------------------------------------------------------------------------
#   minispice-> diodeBiasPoint.py
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

# Example class to calculate the bias point of diode for various
# source impedances
class diode_bias_point:

	def __init__(self, signal):

		# The signal we will simulate
		self.signal = signal

		# Initialize companion models
		self.diode  = componentModels.diode()

		# Initialize companion models
		self.diodeR = companionModels.nonlinearR( self.diode )

	def solve(self, source_impedance, conv = 0.0):
		
		diode_waveform = []

		# Initilize solver
		vm = self.signal['waveform'][0]

		# Perform iteration over companion model
		for source_voltage in self.signal['waveform']:

			# Perform iteration over companion model for nonlinear resistor
			while True:

				num = self.diodeR.im(vm) + nortonI(source_voltage, source_impedance) 
				den = self.diodeR.gm(vm) + nortonG(source_impedance)

				_vm = num / den


				if np.abs( _vm - vm ) <= conv:  

					diode_waveform.append(_vm)

					break 

				else: 

					vm = _vm

		# Return solution waveform
		return diode_waveform

if __name__ == "__main__":


	# Create pulse parameters
	config = {
		"amplitude" : 1.2,
		"period"	: 1e-8,
		"npoints"	: 1e3, 
	}

	# Generate square pulse using signalTools
	signal = signalTools().sinwave(config)

	# Create analysis object
	analysis = diode_bias_point( signal )

	# Create figure for waveform response
	fig0 = plt.figure()
	ax0  = fig0.add_subplot(111)
	ax0.set_xlabel("Time (s)")
	ax0.set_ylabel("Voltage (V)")
	ax0.set_title("Diode Resistor Circuit : Waveform")

	# Create figure for diode bias point
	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111)
	ax1.set_xlabel("Source Voltage (V)")
	ax1.set_ylabel("Diode Voltage (V)")
	ax1.set_title("Diode Resistor Circuit : Bias Point")

	# Create some test resistances(1 Ohm to 1024 Ohm)
	source_impedacnce = [ 2.0**_power for _power in range(11) ]

	# Loop through all source impeances
	for _ in source_impedacnce:
		
		diode_waveform = analysis.solve(source_impedance = _, conv = 1e-15 )
		
		h0, = ax0.plot(signal['time'], diode_waveform, color="tab:orange")
		h1, = ax1.plot(signal['waveform'], diode_waveform, color="tab:blue")

	# Plot the input waveform
	h2, = ax0.plot(signal['time'], signal['waveform'], color="tab:blue")

	# Create a simple legend for both plots
	ax0.legend([h2,h0] , ["Input waveform", "Diode response\n $1\Omega \leq Z_s \leq 1024\Omega$"])
	ax1.legend([h1] , ["Diode bias point \n $1\Omega \leq Z_s \leq 1024\Omega$"])

	# Show the plots
	plt.show()
