# ---------------------------------------------------------------------------------
# 	minispice -> examples/notchFilter.py
#	Copyright (C) 2020 Michael Winters
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

# Import module
from minispice.freqAnalysis import freqAnalysis
from minispice.plotAnalysis import plotAnalysis

# Create a list of frequencies to solve
freq = np.linspace(1,4e7,1000)

# Import spice file and run anlaysis
analysis = freqAnalysis.fromFile('./notchFilter.cir', freq)
voltage_gain = analysis.calcVoltageGain(1, 2)

# Create plot analysis object
plt = plotAnalysis()

# Plot frequency analysis data (magnigude)
key = "MAG" 
plt.add_figure(key)
plt.set_xlabel(key, "Frequency (Hz)")
plt.set_ylabel(key, "Voltage Gain (magnitude)")
plt.set_title(key, "Notch Filter Frequency Analysis")
plt.plot(key, freq, analysis.abs(voltage_gain), "lin")

# Plot frequency analysis data (phase)
key = "PHASE" 
plt.add_figure(key)
plt.set_xlabel(key, "Frequency (Hz)")
plt.set_ylabel(key, "Voltage Gain (phase)")
plt.set_title(key, "Notch Filter Frequency Analysis")
plt.plot(key, freq, analysis.angle(voltage_gain), "lin")

# Show data
plt.show()
