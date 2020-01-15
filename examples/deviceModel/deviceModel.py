# ---------------------------------------------------------------------------------
#   minispice -> examples/deviceModel.py
#   Copyright (C) 2020 Michael Winters
#   github: https://github.com/mesoic
#   email:  mesoic@protonmail.com
# ---------------------------------------------------------------------------------
#   
#   Permission is hereby granted, free of charge, to any person obtaining a copy
#   of this software and associated documentation files (the "Software"), to deal
#   in the Software without restriction, including without limitation the rights
#   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#   copies of the Software, and to permit persons to whom the Software is
#   furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
#   copies or substantial portions of the Software.
#   
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.
#

#!/usr/bin/env python 
from matplotlib import pyplot as plt
import numpy.linalg as la
import numpy as np
import os

# Import device model
from minispice.nonlinear import componentModels

class deviceModel: 

    def __init__( self, size = 50 ):

        # Initalize transistor model 
        self.device = componentModels.HFET()

        # Size of model
        self.size = size

        # Build model data
        self.build()

    # Build model over domain
    def build(self):
    
        # Gate and drain voltages
        self.Vgs = np.linspace(-1.0, 1.0, self.size)
        self.Vds = np.linspace( 0.0, 3.0, self.size)

         # Array to hold the grid
        self.VGS, self.VDS = np.meshgrid( self.Vgs, self.Vds )

        # Arrays to hold current
        self.IDS = np.zeros( (self.size, self.size) )
        
        # Loop through voltages
        for i, _vds in enumerate(self.Vds): 
        
            for j, _vgs in enumerate(self.Vgs): 
            
                # Model looks takes gate drain and gate source voltage
                _vgd = (_vgs -_vds)
                
                # Calculate current
                self.IDS[i][j] = self.device.f( _vgd , _vgs ) 

    # Show model data
    def show(self):
        
        # Plot id(vd) curves (const vg)
        fig = plt.figure(1)
        ax0 = fig.add_subplot(111)
        ax0.set_xlabel("Drain Voltage $(V_{ds})$")
        ax0.set_ylabel("Drain Current $(I_{ds})$")
        ax0.set_title("Nonlinear HFET Model $I_{ds}(V_{ds} : V_{gs} = const)$")

        for i in range(self.size):
        
            ax0.plot(self.Vds, self.IDS[:, i], color="tab:blue" )

        # Plot id(vg) curves (const vd)
        fig = plt.figure(2)        
        ax0 = fig.add_subplot(111)
        ax0.set_xlabel("Gate Voltage $(V_{gs})$")
        ax0.set_ylabel("Drain Current $(I_{ds})$")
        ax0.set_title("Nonlinear HFET Model $I_{ds}(V_{gs} : V_{ds} = const)$")

        for i in range( self.size):
        
            ax0.plot( self.Vgs,  self.IDS[i, :], color="tab:blue" )    

        # Plot contour
        fig = plt.figure(3)
        ax0 = fig.add_subplot(111)        

        h = ax0.contourf(self.VDS, self.VGS, self.IDS )
        ax0.set_title("Nonlinear HFET Model (contour)")
        ax0.set_xlabel('Drain Voltage $(V_{ds})$')
        ax0.set_ylabel('Gate Voltage $(V_{gs})$')
        cbar = fig.colorbar(h)
        cbar.ax.set_ylabel('Drain Current $(I_{ds})$')
        fig.tight_layout()

        plt.show()    
        

if __name__ == "__main__":

    # Initialize and show the model
    model = deviceModel()
    model.show()
   