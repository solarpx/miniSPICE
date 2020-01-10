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

# Create a list of frequencies to solve
freqList = np.logspace(1,12,100)

# Import spice file and run anlaysis
data = freqAnalysis.fromFile('./cascodeAmp.cir', freqList)

# Plot voltage gain (dB)
if True:
	data.plotGain(1,7,"dB")

# Calculate input impedance Rl = 50
if True: 
	data.plotInputImpedance(1, 7, 50., "log")	

# Calculate output impedance Rs = 50
if True:
	data.plotOutputImpedance(1, 7, 50., "log")