# ---------------------------------------------------------------------------------
# 	minispice-> examples/cascodeAmp.py
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

# Import module
from minispice.freqAnalysis import freqAnalysis
from minispice.plotAnalysis import plotAnalysis

# Create a list of frequencies to solve
freq = np.logspace(1,12,100)

# Import spice file and run anlaysis
analysis = freqAnalysis.fromFile('./cascodeAmp.cir', freq)
network_gain = analysis.calcNetworkGain(1, 7, Zs = 4e3, Zl = 4.0e3)

# Calculate input impedance and output impedance
input_impedance = analysis.calcInputImpedance(1, 7, Zl = 4.0e3)
output_impedance = analysis.calcOutputImpedance(1, 7, Zs = 4.0e3)

# Create plot analysis object
plt = plotAnalysis()

# Plot voltage gain (dB)
key = "NETWORK"
plt.add_figure(key)
plt.set_xlabel(key, "Frequency (Hz)")
plt.set_ylabel(key, "Network Gain (dB)")
plt.set_title(key, "Cascode Amplifier Frequency Analysis")
plt.plot(key, freq, network_gain, "log(dB)")

# Plot input and output impedance 
key = "Z_MAG"
plt.add_figure(key)
plt.set_xlabel(key, "Frequency (Hz)")
plt.set_ylabel(key, "|Impedance| (Ohm)")
plt.set_title(key, "Zin(blue) : Zout(orange)")
plt.plot(key, freq, analysis.abs(input_impedance), "loglog")
plt.plot(key, freq, analysis.abs(output_impedance), "loglog")

# Plot input and output impedance 
key = "Z_PHASE"
plt.add_figure(key)
plt.set_xlabel(key, "Frequency (Hz)")
plt.set_ylabel(key, "<Impedance (deg)")
plt.set_title(key, "Zin(blue) : Zout(orange)")
plt.plot(key, freq, analysis.angle(input_impedance), "log")
plt.plot(key, freq, analysis.angle(output_impedance), "log")

# Show plots
plt.show()
