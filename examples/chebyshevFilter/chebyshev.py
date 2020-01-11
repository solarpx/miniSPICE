# ---------------------------------------------------------------------------------
# 	minispice-> examples/plotSmith.py
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
from minispice.freqAnalysis import plotAnalysis
from minispice.plotSmith import plotSmith

# Create a list of frequencies to solve
freq = np.linspace(100, 4.2e3, 1000)

# This example compares two chebyshev filters (order 2 and order 4)
analysis2 = freqAnalysis.fromFile('./chebyshev2.cir', freq)
analysis4 = freqAnalysis.fromFile('./chebyshev4.cir', freq)

# The filters are designed for source and load impdances: 
#	Rs = 1e3 Ohm and Rl = 5e3 Ohm
network_gain2 = analysis2.calcNetworkGain(1, 3, 1e3, 5e3)
network_gain4 = analysis4.calcNetworkGain(1, 5, 1e3, 5e3)

# Plot frequency analysis data
plt = plotAnalysis("lin(dB)")
plt.plotGain(freq, network_gain2)
plt.plotGain(freq, network_gain4)
plt.show()

# The SPICE files do not include the source and load impedances.
# It is thus possible to look examine the behavioud of a bare LC 
# network in the Smith Chart.

# First we need to get the Sparameters of both networks
Sparams2 = analysis2.Sparameters(1, 3)
Sparams4 = analysis4.Sparameters(1, 5)


# We can then look at the parmeters in the Smith chart for both filters
#	s[0,0] = s11
#	s[0,1] = s12
#	s[1,0] = s21
#	s[2,2] = s22
#
pltSmith = plotSmith()
pltSmith.plot( [ s[1,1] for f, s in Sparams2.items() ], color="tab:blue" )
pltSmith.plot( [ s[1,1] for f, s in Sparams4.items() ], color="tab:orange" )
pltSmith.show()
