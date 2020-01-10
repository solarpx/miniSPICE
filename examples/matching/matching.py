#!/usr/bin/env python 
import numpy as np

# Import module
from minispice.activeAnalysis import activeAnalysis
from minispice.freqAnalysis import freqAnalysis
from minispice.plotSmith import plotSmith
from minispice.Converter import *

# We would like to match our intrinsic transistor at 10GHz for maximum gain. 
# The system impedance is 50Ohm. First we need to extract the yparamters
data = freqAnalysis.fromFile("./intrinsic.cir", [1e10])

# Extract the data and compress to twoport
sparams = ytos( data.ygroup[0].toTwoport(1,6) )

# Pass parameters and call the analysis
active  = activeAnalysis(sparams)
active.analysis()

# Plot stability and gain circles
pltS = plotSmith()

# Extract stability circle data for plotting
_input  = active.inputStabilityCircle()["data"]
_output = active.outputStabilityCircle()["data"]

pltS.plot(_input, linestyle=":")
pltS.plot(_output, linestyle=":")

# Constant gain circles (3.0dB to 8.0dB) (grey)
for dB in [3.0, 4.0, 5.0, 6.0, 7.0, 8.0]:
    pltS.plot( active.constantGainCircle(dB)["data"] , color="grey", linewidth=0.8)

# 9.0dB gain circle and conjugate (orange)
_gain = active.constantGainCircle(dB)["data"]
_conj = active.conjugateCircleData(_gain)

pltS.plot( _gain, color="tab:orange" )
pltS.plot( _conj, color="tab:blue" )


# Now we would like to perform the conjugate match. Pick a Gl on the  
# _gain circle and find the associated value on conjugate circle. 
#   Make sure |GammaL| < 1.0 and |GammaS| < 1.0
GammaL = _gain[350]
GammaS = _conj[350]

# Plot points on the smith chart to visualize. 
pltS.plot([GammaL], marker="o", color="tab:orange")
pltS.plot([GammaS], marker="o", color="tab:blue")

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
pltS.show()
