#!/usr/bin/env python 
import numpy as np

# Import module
from minispice.freqAnalysis import freqAnalysis

# Create a list of frequencies to solve
freqList = np.linspace(1,4e7,1000)

# Import spice file and run anlaysis
data = freqAnalysis.fromFile('./notchFilter.cir', freqList)
handle = data.plotGain(1,2,"lin")