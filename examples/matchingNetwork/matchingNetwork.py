# ---------------------------------------------------------------------------------
# 	minispice-> examples/matchingNetwork.py
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
from minispice.amplAnalysis import amplAnalysis
from minispice.Converter import *

# We would like to match our intrinsic transistor at 10GHz for maximum gain. 
# The system impedance is 50Ohm. First we need to extract the yparamters
analysis = freqAnalysis.fromFile("./transistor_model.cir", [1e10])

# Extract the data and compress to twoport
ymatrix = analysis.getMatrix( freq = 1e10 )
sparams = ytos( ymatrix.toTwoport(1,6) )

# Pass parameters and call the analysis
active = amplAnalysis(sparams)
active.analysis()

# Extract stability circle data for plotting
_input  = active.inputStabilityCircle()["data"]
_output = active.outputStabilityCircle()["data"]

# Prepare a plotting object
plt = plotAnalysis()

# Prepare a smith chart
key = "SMITH_CHART"
plt.add_smith(key)
plt.set_title(key, "Active Microwave Circuit Design: Impedance Matching")
plt.plot_smith(key, _input, linestyle=":")
plt.plot_smith(key, _output, linestyle=":")

# Constant gain circles (3.0dB to 8.0dB) (grey)
for dB in [3.0, 4.0, 5.0, 6.0, 7.0, 8.0]:
    plt.plot_smith(key, active.constantGainCircle(dB)["data"] , color="grey", linewidth=0.8)

# 9.0dB gain circle and conjugate (orange)
_gain = active.constantGainCircle(dB)["data"]
_conj = active.conjugateCircleData(_gain)

plt.plot_smith(key, _gain, color="tab:orange" )
plt.plot_smith(key, _conj, color="tab:blue" )

# Now we would like to perform the conjugate match. Pick a Gl on the  
# _gain circle and find the associated value on conjugate circle. 
#   Make sure |GammaL| < 1.0 and |GammaS| < 1.0
GammaL = _gain[350]
GammaS = _conj[350]

# Plot points on the smith chart to visualize. 
plt.plot_smith(key, [GammaL], marker="o", color="tab:orange")
plt.plot_smith(key, [GammaS], marker="o", color="tab:blue")

# Print match data to screen
print("\nConjugate Match (Reflection Coefficients)")
print("\t|GammaL| = %f : <GammaL = %f"%phasor(GammaL, "deg") )
print("\t|GammaS| = %f : <GammaS = %f"%phasor(GammaS, "deg") )

print("\nConjugate Match (Admittance)")
print("\tYs = %s"%(gammatoy(GammaS,1)) )
print("\tYl = %s"%(gammatoy(GammaL,1)) )

# stdout buffer
print("\n")

# Show smith plot
plt.show()
