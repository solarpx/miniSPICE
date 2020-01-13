# ---------------------------------------------------------------------------------
#   minispice-> signalTools.py
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

#!/usr/bin/env python 
from scipy import signal
import numpy as np

# Create signal tools namespace
class signalTools:

	def __init__(self):

		pass

	# Method to generate a test pulse for transient simulations
	def pulse(self, config):

		# Construct time array 
		time  = np.linspace(0, config["period"], config["npoints"], endpoint=False)

		# Cache the sampling interval
		delta = time[1] - time[0]

		# Use the signal library to create a pulse signal
		waveform = (config["amplitude"] / 2) * ( signal.square(2 * np.pi * (1.0/config["period"]) * time - np.pi / 2.0, config["duty"]) + 1.0)

		# Return waveform, time, and sampling interval (delta)
		_signal = {
			"waveform"	: waveform,
			"time"		: time,
			"delta" 	: delta
		} 

		return _signal

	# Method to generate a sin wave (one cycle)
	def sinwave(self, config):

		# Construct time array 
		time  = np.linspace(0, config["period"], config["npoints"], endpoint=False)

		# Cache the sampling interval
		delta = time[1] - time[0]

		# Build sin wave over period
		waveform = [ config["amplitude"] * np.sin(  2 * np.pi * (1.0/config["period"]) * _ ) for _ in time ]

		# Return waveform, time, and sampling interval (delta)
		_signal = {
			"waveform"	: waveform,
			"time"		: time,
			"delta" 	: delta
		} 

		return _signal
