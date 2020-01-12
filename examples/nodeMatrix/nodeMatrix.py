# ---------------------------------------------------------------------------------
# 	minispice-> examples/nodeMatrix.py
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

#!/usr/bin/env python 
import numpy as np

# Import modules
from minispice.nodeMatrix import nodeMatrix
from minispice.plotAnalysis import plotAnalysis

# These examples provide an illustration of how to interact
# directly with the nodeMatrix class.

# Example 1: Resistive pi-Network
# Direct initialization of node admittance matrix
if False:
	
	size = 2   # Size of matrix
	freq = 1.0 # Frequency to simulate
	y = nodeMatrix(size = 2, freq = 1)

	# Add resistors
	y.addPassive("R", 0, 1, 1.0) # 1 Ohm between node 1 and ground
	y.addPassive("R", 1, 2, 2.0) # 2 Ohm between node 1 and node 2
	y.addPassive("R", 2, 0, 1.0) # 1 Ohm between node 2 and ground

	# Show matrix
	y.showMatrix()
	
	# Calculate voltage gain from node 1 to node 2
	print(y.voltageGain(1,2))


# Example 2: Frequency analysis with a reactive component
if True:

	# Admittance matrix size
	size = 2

	# Initialize a list of frequencies (0.1Hz to2kHz)
	freq = np.linspace(.1, 2, 10000)

	# Create arraus to store voltage gain and phase
	mag, phase = [],[]

	for i, f in enumerate(freq):
		y = nodeMatrix(size, f)

		# Resistors
		y.addPassive("R", 0, 1, 1.0) # 1 Ohm between node 1 and ground
		y.addPassive("R", 2, 0, 1.0) # 1 Ohm between node 1 and node 2
		y.addPassive("R", 1, 2, 1.0) # 1 Ohm between node 2 and ground

		# Capatitor
		y.addPassive("C", 1, 2, 2.0) # 2 Farad between node 1 and node 2

		# Calculate voltage gain
		_gain = y.voltageGain(1,2)

		# Extraxt magnitude and phase and add to array
		mag.append( np.abs(_gain) )
		phase.append( np.angle(_gain, deg=True) )


	# Create plot analysis object
	plt = plotAnalysis()
	
	# Plot frequency analysis data
	key = "MAG"
	plt.add_figure(key)
	plt.set_xlabel(key, "Frequency (Hz)")
	plt.set_ylabel(key, "Voltage Gain")
	plt.set_title(key, "High Pass Filter" )
	plt.plot(key, freq, mag, "lin")

	key = "PHASE"
	plt.add_figure(key)
	plt.set_xlabel(key, "Frequency (Hz)")
	plt.set_ylabel(key, "Voltage Phase (deg)")
	plt.set_title(key, "High Pass Filter" )
	plt.plot(key, freq, phase, "lin")

	plt.show()
