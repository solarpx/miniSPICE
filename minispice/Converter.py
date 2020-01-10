#!/usr/bin/env python 
import numpy as np 
import math

# Data check method
def dataCheck(matrix):
    if not isinstance(matrix, np.ndarray):
        return 0 
    elif np.shape(matrix) is (2,2):
        return 0
    else: 
        return 1

# Converters for y-parameters
def ytoz(y):
    if dataCheck(y):
        delta = np.linalg.det(y)
        z00 = y[1][1]/delta
        z01 = -y[0][1]/delta
        z10 = -y[1][0]/delta
        z11 = y[0][0]/delta
        return np.array([[z00,z01],[z10,z11]],dtype='complex')
    else: 
        return None

def ytot(y):
    if dataCheck(y):
        delta = np.linalg.det(y)
        t00 = -y[1][1]/y[1][0]
        t01 = -1/y[1][0]
        t10 = -delta/y[1][0]
        t11 = -y[0][0]/y[1][0]
        return np.array([[t00,t01],[t10,t11]],dtype='complex')
    else: 
        return None

def ytoh(y):
    if dataCheck(y):
        delta = np.linalg.det(y)
        h00 = 1/y[0][0]
        h01 = -y[0][1]/y[0][0]
        h10 = y[1][0]/y[0][0]
        h11 = delta/y[0][0]
        return np.array([[h00,h01],[h10,h11]],dtype='complex')
    else: 
        return None

# Converters for z-parameters
def ztoy(z):
    if dataCheck(z):
        delta = np.linalg.det(z)
        y00 = z[1][1]/delta
        y01 = -z[0][1]/delta
        y10 = -z[1][0]/delta
        y11 = z[0][0]/delta
        return np.array([[y00,y01],[y10,y11]],dtype='complex')
    else: 
        return None

def ztot(z):
    if dataCheck(z):
        delta = np.linalg.det(z)
        t00 = z[0][0]/z[1][0]
        t01 = delta/z[1][0]
        t10 = 1/z[1][0]
        t11 = z[1][1]/z[1][0]
        return np.array([[t00,t01],[t10,t11]],dtype='complex')
    else: 
        return None

def ztoh(z):
    if dataCheck(z):
        delta = np.linalg.det(z)
        h00 = delta/z[1][1]
        h01 = z[0][1]/z[1][1]
        h10 = -z[1][0]/z[1][1]
        h11 = 1/z[1][1]
        return np.array([[h00,h01],[h10,h11]],dtype='complex')
    else: 
        return None

# Converters for T parameters
def ttoz(t):
    if dataCheck(t):
        delta = np.linalg.det(t)
        z00 = t[0][0]/t[1][0]
        z01 = delta/t[1][0]
        z10 = 1/t[1][0]
        z11 = t[1][1]/t[1][0]
        return np.array([[z00,z01],[z10,z11]],dtype='complex')
    else: 
        return None

def ttoy(t):
    if dataCheck(t):
        delta = np.linalg.det(t)
        y00 = t[1][1]/t[0][1]
        y01 = -delta/t[0][1]
        y10 = -1/t[0][1]
        y11 = t[0][0]/t[0][1]
        return np.array([[y00,y01],[y10,y11]],dtype='complex')
    else: 
        return None

def ttoh(t):
    if dataCheck(t):
        delta = np.linalg.det(t)
        h00 = t[0][1]/t[1][1]
        h01 = -delta/t[1][1]
        h10 = -1/t[1][1]
        h11 = t[1][0]/t[1][1]
        return np.array([[h00,h01],[h10,h11]],dtype='complex')
    else: 
        return None

# Converters for H parameters
def htoz(h):
    if dataCheck(h):
        delta = np.linalg.det(h)
        z00 = delta/h[1][1]
        z01 = h[0][1]/h[1][1]
        z10 = -h[1][0]/h[1][1]
        z11 = 1/h[1][1]
        return np.array([[z00,z01],[z10,z11]],dtype='complex')
    else: 
        return None

def htoy(h):
    if dataCheck(h):
        delta = np.linalg.det(h)
        y00 = 1/h[0][0]
        y01 = -h[0][1]/h[0][0]
        y10 = h[1][0]/h[0][0]
        y11 = delta/h[0][0]
        return np.array([[y00,y01],[y10,y11]],dtype='complex')
    else: 
        return None

def htot(h):
    if dataCheck(h):
        delta = np.linalg.det(h)
        t00 = -delta/h[1][0]
        t01 = -h[0][0]/h[1][0]
        t10 = -h[1][1]/h[1][0]
        t11 = -1/h[1][0]
        return np.array([[t00,t01],[t10,t11]],dtype='complex')
    else: 
        return None


# S-parameters to z and y no argument assumes 50 Ohms
def stoy(s, z0 = 50.):
    if dataCheck(s): 

        # Line impedances
        y0= 1/float(z0)
        delta = (1+s[0][0])*(1+s[1][1])-s[0][1]*s[1][0]
        
        # Calculate y
        y00 = (((1-s[0][0])*(1+s[1][1]) + (s[0][1]*s[1][0]))*y0)/delta
        y01 = -2*s[0][1]*y0/delta
        y10 = -2*s[1][0]*y0/delta
        y11 = (((1+s[0][0])*(1-s[1][1]) + (s[0][1]*s[1][0]))*y0)/delta
        return np.array([[y00,y01],[y10,y11]],dtype='complex')
    else: 
        return None
        
def ytos(y, z0 = 50.):
    if dataCheck(y):         
    
        # Normalize y to line impedances   
        delta = (1+z0*y[0][0])*(1+z0*y[1][1])-(z0*z0*y[1][0]*y[0][1])

        # Calculate s
        s00 = ((1-z0*y[0][0])*(1+z0*y[1][1]) + (z0*z0*y[0][1]*y[1][0]))/delta
        s01 = -2*y[0][1]*z0/delta
        s10 = -2*y[1][0]*z0/delta
        s11 = ((1+z0*y[0][0])*(1-z0*y[1][1]) + (z0*z0*y[0][1]*y[1][0]))/delta
        return np.array([[s00,s01],[s10,s11]],dtype='complex')
    else: 
        return None
 
def stoz(s, z0 = 50.):
    if dataCheck(s): 

        # Line impedances
        delta = (1-s[0][0])*(1-s[1][1])-s[0][1]*s[1][0]
        
        # Calculate y
        z00 = ((1+s[0][0])*(1-s[1][1])+(s[0][1]*s[1][0]))*z0/delta
        z01 = 2*s[0][1]*z0/delta
        z10 = 2*s[1][0]*z0/delta
        z11 = ((1-s[0][0])*(1+s[1][1])+(s[0][1]*s[1][0]))*z0/delta
        return np.array([[z00,z01],[z10,z11]],dtype='complex')
    else: 
        return None


def ztos(z, z0 = 50.):
     if dataCheck(z):         
    
        # Normalize y to line impedances   
        delta = (z[0][0]+z0)*(z[1][1]+z0)-(z[1][0]*z[0][1])

        # Calculate s
        s00 = ((z[0][0]-z0)*(z[1][1]+z0)-(z[0][1]*z[1][0]))/delta
        s01 = 2*z[0][1]*z0/delta
        s10 = 2*z[1][0]*z0/delta
        s11 = ((z[0][0]+z0)*(z[1][1]-z0)-(z[0][1]*z[1][0]))/delta
        return np.array([[s00,s01],[s10,s11]],dtype='complex')
     else: 
        return None

# S to transfer scattering parameters
def stor(s): 
    if dataCheck(s):
        delta = np.linalg.det(s)
        r00 = -delta/s[1][0]
        r01 = s[0][0]/s[1][0]
        r10 = -s[1][1]/s[1][0]
        r11 = 1/s[1][0] 
        return np.array([[r00,r01],[r10,r11]],dtype='complex')
    else:
        return None

def rtos(r): 
    if dataCheck(r):
        delta = np.linalg.det(r)
        s00 = r[0][1]/r[1][1]
        s01 = delta/r[1][1]
        s10 = 1/r[1][1]
        s11 = -r[1][0]/r[1][1] 
        return np.array([[s00,s01],[s10,s11]],dtype='complex')
    else:
        return None

# Z to Gamma and Gamma to Z
def gammatoz(gamma,z0=50.0):
    return complex(z0*(1.0+complex(gamma))/(1.0-complex(gamma)))

def ztogamma(z,z0=50.0):
    return complex((z-z0)/(z+z0))

# Y to Gamma and Gamma to Y
def gammatoy(gamma,y0=0.02):
    return complex(y0*(1-complex(gamma))/(1+complex(gamma)))

def ytogamma(y,y0=0.02):
    return complex((y-y0)/(y+y0))


# A few smith chart things
def seriesL(x,f,z0=50.):
    w = 2*math.pi*f
    return x*z0/w

def seriesC(x,f,z0=50.):
    w = 2*math.pi*f
    return 1./(w*z0*x)

def shuntL(b,f,z0=50.):
    w = 2*math.pi*f
    return z0/(w*b)

def shuntC(b,f,z0=50.):
    w = 2*math.pi*f
    return b/(w*z0)

# Convert complex number to phasor
def phasor(c, _type="rad"):    

    # Radians
    if _type == "rad":
        return ( abs(c), np.angle(c) )

    # Degrees
    if _type == "deg":
        return ( abs(c), np.angle(c)*180./math.pi )

# scalar to dB
def todB(val):
    return 10*np.log10(val)

# dB to scalar
def fromDb(val):
    return np.power(10.0, val/10.)

# A test program
if __name__=="__main__":
    z = np.array([[1,2],[3,4]],dtype='complex')
    print(htoz(ytoh(ttoy(htot(ytoh(stoy(rtos(stor(ztos(z))))))))))
