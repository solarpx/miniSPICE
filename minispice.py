########################
# A miniSPICE emulator #
########################

# File IO handling
from __future__ import with_statement
from io import *
import os

# Classes for array manipulation
import re
import numpy as np
import scipy as scipy 
import math as math
import copy

# For plotting
import matplotlib.pyplot as plt

# Base Class for nodal admittance matrix
class nodeMatrix(object): 

    def __init__(self,size,freq): 
        self.ymatrix = np.zeros(shape=(size,size),dtype=complex)
        self.freq = float(freq)
        self.size = size

    def addPassive(self,name,n1,n2,value):
        # Convert to float automatically
        value = float(value)

        # Resistances (g = 1/r)
        g = lambda r: complex(1/r) 
        if re.match(r'R\d*',name) is not None:
            if n1 == 0:  
                n2-=1
                self.ymatrix[n2,n2] += g(value)
            elif n2 == 0:  
                n1-=1
                self.ymatrix[n1,n1] += g(value)
            else:
                n1-=1
                n2-=1
                self.ymatrix[n1,n1] += g(value)
                self.ymatrix[n1,n2] -= g(value)
                self.ymatrix[n2,n1] -= g(value)
                self.ymatrix[n2,n2] += g(value)

        # Capacitances (bc = iwC)
        bc = lambda c : complex(0,2*math.pi*self.freq*c) 
        if re.match(r'C\d*',name) is not None:
            value = complex(value,0)
            if n1 == 0:  
                n2-=1
                self.ymatrix[n2,n2] += bc(value)
            elif n2 == 0:  
                n1-=1
                self.ymatrix[n1,n1] += bc(value)
            else:
                n1-=1
                n2-=1
                self.ymatrix[n1,n1] += bc(value)
                self.ymatrix[n1,n2] -= bc(value)
                self.ymatrix[n2,n1] -= bc(value)
                self.ymatrix[n2,n2] += bc(value)

        # Inductances (bl = iwL)
        bl = lambda l : complex(0,(-1/(2*math.pi*self.freq*l))) 
        if re.match(r'L\d*',name) is not None:
            value = complex(value,0)
            if n1 == 0:  
                n2-=1
                self.ymatrix[n2,n2] += bl(value)
            elif n2 == 0:  
                n1-=1
                self.ymatrix[n1,n1] += bl(value)
            else:
                n1-=1
                n2-=1
                self.ymatrix[n1,n1] += bl(value)
                self.ymatrix[n1,n2] -= bl(value)
                self.ymatrix[n2,n1] -= bl(value)
                self.ymatrix[n2,n2] += bl(value)

    
    # Method for adding a transistor
    def addTransistor(self,name, nb, nc, ne, model): 
        
        nb-=1
        nc-=1
        ne-=1

        # Simple Model with only b and rbe
        g = lambda r : complex(1/r) 
        if model == 'simple':
            # Extract transistor parameters out of params dict            
            params = self.getModel('simple')
            b=float(params['b'])
            rbe=float(params['rbe'])
            # Base
            self.ymatrix[nb,nb]+=g(rbe)
            self.ymatrix[nb,nc]+=0
            self.ymatrix[nb,ne]-=g(rbe)
            # Collector
            self.ymatrix[nc,nb]+=b*g(rbe)
            self.ymatrix[nc,nc]+=0
            self.ymatrix[nc,ne]-=b*g(rbe)
            # Emitter
            self.ymatrix[ne,nb]-=(b+1)*g(rbe)
            self.ymatrix[ne,nc]+=0
            self.ymatrix[ne,ne]+=(b+1)*g(rbe)

        # Intrinsic transistor pi model
        g  = lambda r : complex(1/r) 
        bc = lambda c : complex(0,2*math.pi*self.freq*c) 
        if model == 'hybridpi':
            # Extract transistor parameters out of params dict            
            params = self.getModel('hybridpi')
            
            gm = float(params['gm'])
            r0 = float(params['rce'])
            rpi= float(params['rbe'])
            cpi= float(params['cbei'])
            cmu= float(params['cbc'])

            # Base
            self.ymatrix[nb,nb]+=(g(rpi)+bc(cpi)+bc(cmu))
            self.ymatrix[nb,nc]+=0
            self.ymatrix[nb,ne]-=(g(rpi)+bc(cpi)+bc(cmu)) 
            # Collector
            self.ymatrix[nc,nb]+=(gm-bc(cmu))
            self.ymatrix[nc,nc]+=g(r0)
            self.ymatrix[nc,ne]+=(bc(cmu)-g(r0)-gm)
            # Emitter
            self.ymatrix[ne,nb]-=((g(rpi)+bc(cpi)+gm))
            self.ymatrix[ne,nc]-=g(r0)
            self.ymatrix[ne,ne]+=(g(rpi)+bc(cpi)+g(r0)+gm)


        #Transistor with base spreading resistance
        g  = lambda r : complex(1/r) 
        bc = lambda c : complex(0,2*math.pi*self.freq*c) 
        if model == 'hybridpix':
            # Extract transistor parameters out of params dict            
            params = self.getModel('hybridpix')
            
            gm = float(params['gm'])
            rce = float(params['rce'])
            rbe= float(params['rbe'])
            cbe= float(params['cbe'])
            cbc= float(params['cbc'])
            rbb = float(params['rbb'])

            # Construct the CE admittance parameters
            y11 = (g(rbe)+bc(cbe)+bc(cbc))
            y12 = 0
            y21 = (gm-bc(cbc))
            y22 = g(rce)
            rbbDce = (rbb*((y11*y22)-(y21*y12)))
            s = (1/(1+y11*rbb))

            # Base
            self.ymatrix[nb,nb]+=((y11)*s)
            self.ymatrix[nb,nc]+=((y12)*s)
            self.ymatrix[nb,ne]-=((y11+y12)*s)
            # Collector
            self.ymatrix[nc,nb]+=((y21)*s)
            self.ymatrix[nc,nc]+=((y22+rbbDce)*s)
            self.ymatrix[nc,ne]-=((y21+y22+rbbDce)*s)
            # Emitter
            self.ymatrix[ne,nb]-=((y11+y21)*s)
            self.ymatrix[ne,nc]-=((y12+y22+rbbDce)*s)
            self.ymatrix[ne,ne]+=((y11+y22+y12+y21+rbbDce)*s)


    # Method to extract params from a *.model file 
    def getModel(self,name):
        params = {}
        path=os.getcwd()+'/model/'+name+'.model'
        with open(path, 'r') as f:
            data = [line.split() for line in f]
        for i,lst in enumerate(data):
            params[str(lst[0])] = float(lst[1])
        return params

    # Method to calculate cofactors Dij
    def cofactorN(self,i,j): 
        try:
            if i<1 or j<1:
                print "Invalid Cofactor Index"
                raise ValueError

        except ValueError:
            return None

        toPower=i+j
        i-=1
        j-=1

        # Store copy of Y matrix and remove rows i and j
        # Need deep copy so operations on A dont change ymatrix
        A = copy.deepcopy(self.ymatrix)
        A[i,:] = np.NaN
        A[:,j] = np.NaN
        
        Am = np.ma.masked_invalid(A)
        Am = np.ma.compressed(Am)
        Am = np.array(chunks(Am,self.size-1))
        return np.linalg.det(Am)*((-1)**toPower)
    
    # Method to calculate cofactor Dii,jj
    def cofactorD(self,i,j): 
        
        #Check bounds of cofactor index
        try:
            if i<1 or j<1:
                print "Invalid Cofactor Index"
                raise ValueError
        except ValueError:
            return None

        #In case of twoport return 1
        try: 
            if self.size == 2:
                raise ValueError
        except ValueError:
            return 1

        i-=1
        j-=1
        # Store copy of Y matrix and delete rows using NaN mask
        A = copy.deepcopy(self.ymatrix)
        A[i,:] = np.NaN
        A[j,:] = np.NaN
        A[:,i] = np.NaN
        A[:,j] = np.NaN

        # Compress out the remaning data
        Am = np.ma.masked_invalid(A)
        Am = np.ma.compressed(Am)
        Am = np.array(chunks(Am,self.size-2))
        
        # Am is now the submatrix we desire generate determinant
        return np.linalg.det(Am)

    # Method which calculates cofactors and returns the corresponding twoport parameters
    def toTwoport(self, n1, n2):

        # Initialize 2x2 matrix of zeros
        twoport = np.zeros(shape=(2,2),dtype=complex)

        # Calculate Cofactors and generate Y_twoport (p.122)
        twoport[1,1] += self.cofactorN(n1,n1)
        twoport[1,0] += self.cofactorN(n1,n2)*(-1)
        twoport[0,1] += self.cofactorN(n2,n1)*(-1)
        twoport[0,0] += self.cofactorN(n2,n2)
        Delta = self.cofactorD(n1,n2)
        twoport=twoport/Delta
        return twoport

    def voltageGain(self, n1, n2):
        #The voltage gain is given by the ratio of cofactors
        return self.cofactorN(n1,n2)/self.cofactorN(n1,n1)


# Class to construct y-matrix for a list of frequencies
class frequencyAnalysis(object): 
    
    def __init__(self, ygroup, freqList):
        self.ygroup = ygroup
        self.freqList = freqList

    @classmethod
    def fromFile(cls, path, freqList):
        # Read all components into data        
        with open(path, 'r') as f:
            data = [line.split() for line in f]

        # Get size of ymatrix from file
        size = 0
        for j,lst in enumerate(data):
            
            # Case of Transitor in data file
            if len(lst) > 4:
                if size < int(lst[1]):
                    size = int(lst[1])
                if size < int(lst[2]):
                    size = int(lst[2])
                if size < int(lst[3]):
                    size = int(lst[3])

            # Case of Passives in data file
            else:
                if size < int(lst[1]):
                    size = int(lst[1])
                if size < int(lst[2]):
                    size = int(lst[2])

        ygroup = []
        for i, freq in enumerate(freqList):
            tmp = nodeMatrix(size,freq)
            for j, lst in enumerate(data): 
                # Case of passives
                if re.match(r'R|C|L',str(lst[0])) is not None:
                    tmp.addPassive(str(lst[0]), int(lst[1]), int(lst[2]), float(lst[3]))
                # Case of transistor
                if re.match(r'Q', str(lst[0])) is not None: 
                    tmp.addTransistor(str(lst[0]), int(lst[1]), int(lst[2]), float(lst[3]), str(lst[4]))

            ygroup.append(tmp)
        return cls(ygroup, freqList)

    def plotGain(self,n1,n2, arg=None):
        
        # Obtain nodal gain for each object in y-group
        gain, phase = [],[]
        for i,ymatrix in enumerate(self.ygroup):
            gain.append(np.abs(self.ygroup[i].voltageGain(n1,n2)))
            phase.append(np.angle(self.ygroup[i].voltageGain(n1,n2),deg=True))

        # Long ugly plotting methods
        if arg=="log":
            titleStr="vm("+str(n2)+")/"+"vm("+str(n1)+")"
            plt.figure(1)
            plt.title(titleStr)
            plt.semilogx(self.freqList,gain)
            plt.xlabel("Frequency")
            plt.grid(True, which="majorminor",ls="-", color='0.65')
            plt.ylabel("Voltage Gain")
            
            titleStr="vp("+str(n2)+")"
            plt.figure(2)
            plt.semilogx(self.freqList,phase)
            plt.title(titleStr)
            plt.xlabel("Frequency")
            plt.ylabel("Voltage Phase Difference")
            plt.grid(True, which="majorminor",ls="-", color='0.65')
            plt.show()

        if arg=="lin":
            titleStr="vm("+str(n2)+")/"+"vm("+str(n1)+")"
            plt.figure(1)
            plt.title(titleStr)
            plt.plot(self.freqList,gain)
            plt.xlabel("Frequency")
            plt.grid(True, which="majorminor",ls="-", color='0.65')
            plt.ylabel("Voltage Gain")
            
            titleStr="vp("+str(n2)+")"
            plt.figure(2)
            plt.plot(self.freqList,phase)
            plt.title(titleStr)
            plt.xlabel("Frequency")
            plt.ylabel("Voltage Phase Difference")
            plt.grid(True, which="majorminor",ls="-", color='0.65')
            plt.show()


    def inputImpedance(self,n1,n2,rl,arg):
        zMag,zPhase = [],[]
        g  = lambda r : complex(1/r) 
        for i,ymatrix in enumerate(self.ygroup):
            twoport=ymatrix.toTwoport(n1,n2)
            tmp = (np.linalg.det(twoport)+twoport[0,0])/(twoport[1,1]+g(rl))
            zMag.append(np.abs(1/tmp))
            zPhase.append(np.angle((1/tmp),deg=True))
        
        # Long ugly plotting methods
        if arg=="log":
            titleStr="mag(Zin)"
            plt.figure(1)
            plt.title(titleStr)
            plt.semilogx(self.freqList,zMag)
            plt.xlabel("Frequency")
            plt.grid(True, which="majorminor",ls="-", color='0.65')
            plt.ylabel("Input Impedance Magnitude")
            
            titleStr="phase(Zin)"
            plt.figure(2)
            plt.title(titleStr)
            plt.semilogx(self.freqList,zPhase)
            plt.xlabel("Frequency")
            plt.grid(True, which="majorminor",ls="-", color='0.65')
            plt.ylabel("Input Impedance Phase")
            plt.show()
            
        if arg=="lin":
            titleStr="mag(Zin)"
            plt.figure(1)
            plt.title(titleStr)
            plt.plot(self.freqList,zMag)
            plt.xlabel("Frequency")
            plt.grid(True, which="majorminor",ls="-", color='0.65')
            plt.ylabel("Input Impedance Magnitude")
            
            titleStr="phase(Zin)"
            plt.figure(2)
            plt.title(titleStr)
            plt.plot(self.freqList,zPhase)
            plt.xlabel("Frequency")
            plt.grid(True, which="majorminor",ls="-", color='0.65')
            plt.ylabel("Input Impedance Phase")
            plt.show()

    def outputImpedance(self,n1,n2,rs,arg):
        zMag,zPhase = [],[]
        g  = lambda r : complex(1/r) 
        for i,ymatrix in enumerate(self.ygroup):
            twoport=ymatrix.toTwoport(n1,n2)
            tmp = (np.linalg.det(twoport)+twoport[1,1])/(twoport[0,0]+g(rs))
            zMag.append(np.abs(1/tmp))
            zPhase.append(np.angle((1/tmp),deg=True))
           

        # Long ugly plotting methods
        if arg=="log":
            titleStr="mag(Zout)"
            plt.figure(3)
            plt.title(titleStr)
            plt.semilogx(self.freqList,zMag)
            plt.xlabel("Frequency")
            plt.grid(True, which="majorminor",ls="-", color='0.65')
            plt.ylabel("Output Impedance Magnitude")
            
            titleStr="phase(Zout)"
            plt.figure(4)
            plt.title(titleStr)
            plt.semilogx(self.freqList,zPhase)
            plt.xlabel("Frequency")
            plt.grid(True, which="majorminor",ls="-", color='0.65')
            plt.ylabel("Output Impedance Phase")
            plt.show()

        if arg=="lin":
            titleStr="mag(Zout)"
            plt.figure(3)
            plt.title(titleStr)
            plt.plot(self.freqList,zMag)
            plt.xlabel("Frequency")
            plt.grid(True, which="majorminor",ls="-", color='0.65')
            plt.ylabel("Output Impedance Magnitude")
            
            titleStr="phase(Zout)"
            plt.figure(4)
            plt.title(titleStr)
            plt.plot(self.freqList,zPhase)
            plt.xlabel("Frequency")
            plt.grid(True, which="majorminor",ls="-", color='0.65')
            plt.ylabel("Output Impedance Phase")


def chunks(l, n):
    return [list(l[i:i+n]) for i in range(0, len(l), n)]

#############################################
# Various examples of usage of minispice.py #
#############################################
#
# Note ... these go in a separate file
# Then minispice is invoked via the following
#
# import minispice 
#
#
# End Comment

#############################################
# Example pi-Network: Direct Initialization #
#############################################
#freq = 1
#size = 2     
#y = nodeMatrix(size, freq)

#name = "R"
#y.addPassive(name,0,1,1)
#y.addPassive(name,1,2,1)
#y.addPassive(name,2,0,1)
#print x.ymatrix
#print x.toTwoport(1,2)
#print x.voltageGain(1,2)


####################################################
# Example pi-Network with C: Direct Initialization #
####################################################
#size = 2     
#freqList = np.linspace(.1, 2, 10000)
#gain,phase = [],[]
#for i,freq in enumerate(freqList):
#    x = nodeMatrix(size, freq)
#    name = "R"
#    x.addPassive(name,0,1,1)
#    x.addPassive(name,2,0,1)
#    x.addPassive(name,1,2,1)
#    name = "C"
#    x.addPassive(name,1,1,2)
#    gain.append(np.abs(x.voltageGain(1,2)))
#    phase.append(np.angle(x.voltageGain(1,2),deg=True))

#plt.figure(1)
#plt.plot(freqList,gain)
#plt.figure(2)
#plt.plot(freqList,phase)

#############################################
# Example Transistor: Direct Initialization #
#############################################
#size = 5
#freq = 1
#params = {'b': 5, 'rbe': 1}
#y = nodeMatrix(size,freq)
#y.addTransistor("Q",1,3,2,"simple")

####################################################
# For Frequency Analysis: Initialization from File #
####################################################
#freqList = np.linspace(.1, 2, 200)
#path='/home/thoth/Desktop/python/minispice/circuits/piWithC.cir'
#data = frequencyAnalysis.fromFile(path, freqList)
#data.plotGain(1,2)


####################################################################
# For Frequency Analysis: Initialization from File With Transistor #
####################################################################
#freqList = np.linspace(.1, 2, 200)
#path='/home/thoth/Desktop/python/minispice/circuits/commonEmitter.cir'
#x = frequencyAnalysis.fromFile(path, freqList)
#print x.ygroup[0].ymatrix
#x.plotGain(1,3)

####################################
# For Frequency Analysis: Cascode  #
####################################
#freqList = np.logspace(0,10,1000)
#path='/home/thoth/Desktop/python/minispice/circuits/cascode_nosl.cir'
#data = frequencyAnalysis.fromFile(path, freqList)
#data.inputImpedance(1,7,4000.,"log")
#data.outputImpedance(1,7,4000.,"log")
#data.plotGain(1,7,"log")

####################################
# For Frequency Analysis: BandStop #
####################################
#freqList = np.linspace(1e6,40e6,1000)
#path='/home/thoth/Desktop/python/minispice/circuits/bandStop.cir'
#data = frequencyAnalysis.fromFile(path, freqList)
#matrix=data.ygroup[0]
#data.plotGain(1,2)
#print matrix.ymatrix
#print matrix.cofactorD(1,2)
#print matrix.cofactorN(1,1)
