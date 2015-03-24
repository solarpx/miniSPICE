#!/usr/bin/env python 
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.mpl as mpl
import numpy.linalg as la
import matplotlib.cm as cm
import numpy as np
import math
import os

## Custom modules
import microwaveConverter as mwc
import plotSmith as ps
import minispice as ms
import ndfit

def _complex(a, phi):
    return a*np.exp(complex(0,phi*math.pi/180))
    
def getCircle(r,c):
    angles = np.linspace(0,2*math.pi,500)
    circle = []
    for a in angles:
        circle.append(r*np.exp(complex(0,a)) + c)
    return circle

def getGainCircle(gp,D,K):

    ## Gain circles
    _cpa = gp*(s22-D*s11.conj()).conj()
    _cpb = 1+gp*(abs(s22)**2 - abs(D)**2)
    _cg  = _cpa/_cpb
    
    _rpa = math.sqrt(1 - (2*K*gp*abs(s12*s21)) + (gp*abs(s12*s21))**2)
    _rpb = abs(1 + gp*(abs(s22)**2 - abs(D)**2))
    _rg  = _rpa/_rpb

    return getCircle(_rg,_cg), _rg,_cg

## We would like to match our transistor at
## 10GHz for maximum gain. The system impedance
## is 50Ohm. First we need to get the yparamters

freqList = [1e10]
path=os.getcwd()+"/intrinsic.cir"
data = ms.frequencyAnalysis.fromFile(path, freqList)
yint =data.ygroup[0].toTwoport(1,6)

## Convert to S parameters
sint = mwc.ytos(yint,50.0)
s11,s12,s21,s22 = sint[0,0],sint[0,1],sint[1,0],sint[1,1]
print sint

#s11 = _complex(0.5,-180)
#s12 = _complex(0.08,30)
#s21 = _complex(2.5,70)
#s22 = _complex(0.8,-100)
               
## Stability Considerations
D = (s11*s22-s21*s12)
K = (1-(abs(s11)**2)-(abs(s22)**2)+(abs(D)**2))/(2*abs(s12*s21))
print "|s11| = %s : <s11 = %s"%(str(abs(s11)),str(np.angle(s11)*180/math.pi))
print "|s12| = %s : <s12 = %s"%(str(abs(s12)),str(np.angle(s12)*180/math.pi))
print "|s21| = %s : <s21 = %s"%(str(abs(s21)),str(np.angle(s21)*180/math.pi))
print "|s22| = %s : <s22 = %s"%(str(abs(s22)),str(np.angle(s22)*180/math.pi))
print "|Delta| = %s : <D = %s : K = %s"%(str(abs(D)),str(np.angle(D)*180/math.pi),str(K))

#################################
#      Stability Circles        #
#################################
rs = abs(s12*s21)/(abs(s11)**2-abs(D)**2)
cs = (s11-(D*s22.conj())).conj()/(abs(s11)**2-abs(D)**2)

rl = abs(s12*s21)/(abs(s22)**2-abs(D)**2)
cl = (s22-(D*s11.conj())).conj()/(abs(s22)**2-abs(D)**2)

inputC  = getCircle(rs,cs)
outputC = getCircle(rl,cl)
print "rs = %s : |cs| = %s : <cs =%s"%(str(abs(rs)),str(abs(cs)),str(np.angle(cs)*180/math.pi))
print "rl = %s : |cl| = %s : <cl =%s"%(str(abs(rl)),str(abs(cl)),str(np.angle(cl)*180/math.pi))

##############################
#    Gain considerations     #
##############################

# Stable case
if K>1:
    print "---------Device Stable---------"
    Gmx = (abs(s21)/abs(s12))*(K - math.sqrt(K**2 -1))
    gmx = (1./abs(s12*s21))*(K - math.sqrt(K**2 -1))
    Gmxdb  = 10*np.log10(Gmx)
    print "Gmx %f = %f dB : gmx = %f"%(Gmx, Gmxdb, gmx)

# Unstable case
else:
    print "---------Device Unstable---------"
    Gmx = (abs(s21)/abs(s12))
    gmx = (1./abs(s12*s21))
    Gmxdb  = 10*np.log10(Gmx)
    print "Gmx %f = %f dB : gmx = %f"%(Gmx, Gmxdb, gmx)

# Define the desired gain
GdB  = 9.0 # <-------dB 
G   = np.power(10.0, GdB/10.)
g    = G/(abs(s21)**2)
print "Gop %f = %f dB : gop = %f"%(G, GdB, g)

## Calculate the 10dB gain circle
gainC,rg,cg = getGainCircle(g,D,K)
print "rg = %s : |cg| = %s : <cg =%s"%(str(abs(rg)),str(abs(cg)),str(np.angle(cg)*180/math.pi))

#################################
#      Plot Output Circles      #
#################################
plt = ps.plotSmith()
plt.plotSmith(inputC,'data','--','b')
plt.plotSmith(outputC,'data','--','r')
plt.plotSmith(gainC,'data','-','r')


#################################
#     More Gain Circles         #
#################################
_GdB = [3.,4.,5.,6.,7.,8.] # <-------dB 
_G   = [10**(i/10) for i in _GdB] 
_g   = [i/(abs(s21)**2) for i in _G]

for gp in _g: 
    _gainC,_rg,_cg = getGainCircle(gp,D,K)
    plt.plotSmith(_gainC,'data','--','0.65','1')

###############################
#   Perform Conjugate Match   #
###############################

# select a value of Gamma_l on or inside 
# the gain circle. Be sure that it is in
# the stable region of the smith chart.
Gl = gainC
Gs = [(s11 + (s12*s21*i)/(1-(s22*i))).conj() for i in Gl]

## Example Check
#__Gl = _complex(0.1,97)
#__Gs = (s11 + (s12*s21*__Gl)/(complex(1)-s22*i))
#print "|Gs| = %s  :  <Gs  = %s"%(str(abs(__Gs)),str(np.angle(__Gs)*180/math.pi))

print "------------------------------------------------"

# Now simply pick a Gamma_l and find the associated
# Yl Gamma_s and Ys. This is our match.

_Gl = gainC[300]
_Gs = (s11 + (s12*s21*_Gl)/(complex(1.0)-s22*_Gl)).conj()
print "|Gs| = %s  :  <Gs  = %s"%(str(abs(_Gs)),str(np.angle(_Gs)*180/math.pi))
print "|Gl| = %s  :  <Gl  = %s"%(str(abs(_Gl)),str(np.angle(_Gl)*180/math.pi))

print "------------------------------------------------"
print "match Zs = %s"%(str(mwc.gammatoz(_Gs,1)))
print "match Zl = %s"%(str(mwc.gammatoz(_Gl,1)))
print "------------------------------------------------"
print "Check to find a Gamma_s is in the stable region!"

plt.plotSmith(Gs,'data','-','b')
plt.plotSmith([_Gs],'data','o','b')
plt.plotSmith([_Gl],'data','o','r')
plt.show()
