#!/usr/bin/env python 

#######################################################################
#
# $Id: PlotSmith.py,v 1.3 2008/08/18 06:13:14 lorenz Exp $
# 
# Copyright (C) 2004-2013 Petr Lorenz, Michael Winters
# All rights reserved.
# 
# ##### BEGIN LICENSE BLOCK #####
#
# Version: MPL 1.1 / GPL 2.0 / LGPL 2.1
#
# The contents of this file are subject to the Mozilla Public License
# Version 1.1 (the "License"); you may not use this file except in
# compliance with the License. You may obtain a copy of the License at
# http://www.mozilla.org/MPL/
# 
# Software distributed under the License is distributed on an "AS IS"
# basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
# License for the specific language governing rights and limitations
# under the License.
# 
# Copyright (C) 2013. All Rights Reserved.
#
# Contributor(s): Michael Winters (GTK version/fully rebuilt)
#               : Sherif Sayed Ahmed (2008 wxPython version)
#               : Petr Lorenz (initial author wxPython version 2008)
#
# Alternatively, the contents of this file may be used under the terms of
# either the GNU General Public License Version 2 or later (the  "GPL"), or the
# GNU Lesser General Public License Version 2.1 or later (the "LGPL"), in which
# case the provisions of the GPL or the LGPL are applicable instead of those
# above.  If you wish to allow use of your version of this file only under the
# terms of either the GPL or the LGPL and not to allow others to use your
# version of this file under the MPL, indicate your decision by deleting the
# provisions above and replace them with the notice and other provisions
# required by the GPL or LGPL.  If you do not delete the provisions above, a
# recipient may use your version of this file under the terms of any one of the
# MPL License, the GPL or the LGPL.
#
# ##### END LICENSE BLOCK #####
#
#######################################################################
# 
# $Log: PlotSmith.py,v $
#
# Revision 2.1  2013/04/11 00:58:30  winters
#   - Rebuilt in GTKAgg backend
#   - Cleaned up axes with GTK
#   - Added methods for subplotting 
#   - Revised internal structure of the plotter
#
# Revision 1.3  2008/08/18 06:13:14  lorenz
# Changes related to 0.7.5 testing version:
#   - new docking window manager
#   - changed calculation of S-Parameters
#   - triangulation of polygons
#
# Revision 1.2  2008-04-07 05:09:35  lorenz
# Smith diagram integration.
#
# Revision 1.1  2008-03-28 06:18:27  lorenz
# Initial commit.
#
#
#######################################################################
#
# March 2013
#
# 2D plotting with matplotlib
#
# Import GTK
import pygtk 
pygtk.require('2.0')
import gtk

# Import matplotlib gtk backend
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
from matplotlib.figure import Figure
## For Plotting
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
rc('text', usetex=True)
rc('xtick', labelsize=20) 
rc('ytick', labelsize=20)
rc('legend', fontsize=20) 


# to calculate the phases
from numpy import *
import microwaveConverter as mwc

class plotSmith:

  # Destroy callback. Destroys gtk loop when 
  # "x" button is clicked.
  def destroy(self, widget, data=None):
    gtk.main_quit()

  def __init__(self):

    # GTK startup
    self.win = gtk.Window(gtk.WINDOW_TOPLEVEL)
    self.win.connect("destroy", self.destroy)
    self.win.set_default_size(850,600)
    self.win.set_title("Smith Chart")
        
    # create the matplotlib-based Plot and pack the
    # figure canvas into the window
    self.figure = Figure()
    self.figure_canvas = FigureCanvas(self.figure)

    # Create Vbox to pack items
    self.main_vbox = gtk.VBox()
    self.win.add(self.main_vbox)

    # Pack the figure canvas into vbox
    self.main_vbox.pack_start(self.figure_canvas, True, True)
    self.toolbar = NavigationToolbar(self.figure_canvas,self.win)
    self.main_vbox.pack_start(self.toolbar, False, False)
    
    # Set an empty list for plot handles
    self.hlist,self.tlist=[],[]

  def show(self):
    self.win.show_all()
    gtk.main()

  def set_legend(self):
    # Note that it is outside the axes. 
    self.axes.legend(self.hlist,self.tlist, bbox_to_anchor=(1.05, 1), loc=2)
    
  def plotSmith(self,s_data,string="data",_ls=None,c=None,lw=None):
    # Make a larger window 
    self.win.set_default_size(800,800)
    
    if _ls is None:
      _ls = '-'
    if c is None:
      c='k'
    if lw is None:
      lw=2
      
    # Generate an axes object to plot into 
    self.axes = self.figure.add_axes([0.1, 0.1, 0.8, 0.8],polar=True)
 
    # Draw the Smith Chart circles
    self.plotsmithcircles(self.axes)
    self.axes.hold(True)

    # Plot the data and set r = 1
    theta,r = [angle(s) for s in s_data],[abs(s) for s in s_data]
    h1 = self.axes.plot(theta,r,_ls, color=c, linewidth=float(lw))
    
    # Append to list for legend
    self.hlist.append(h1)
    self.tlist.append(string)

    # Show freq limits
    #self.axes.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
    #self.axes.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)
    self.axes.set_rmax(1)   
     
  def plotPolar(self, s_data, string="data"):
    # Generate an axes object to plot into 
    self.axes = self.figure.add_axes([0.0, 0.0, 1.0, 1.0],polar=True)
    self.axes.hold(True)
    
    # Plot the data and set r = 1
    theta,r = [angle(s) for s in s_data],[abs(s) for s in s_data]
    h1 = self.axes.plot(theta,r,linewidth=3)
    
    # Append to list for legence
    self.hlist.append(h1)
    self.tlist.append(string)
   
    # Show freq limits
    self.axes.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
    self.axes.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)

  def plotTwoport(self, s, string="data"):

    # Make a larger window 
    self.win.set_default_size(1300,1000)

    # Set up the axes and subplots
    # s11,s22 in smith chart
    # s12,s21 in polar chart
    self.axes1 = self.figure.add_subplot(2,2,1,projection='polar')
    self.plotsmithcircles(self.axes1,True)
    
    self.axes2 = self.figure.add_subplot(2,2,2,projection='polar')
    self.axes3 = self.figure.add_subplot(2,2,3,projection='polar')
    
    self.axes4 = self.figure.add_subplot(2,2,4,projection='polar')
    self.plotsmithcircles(self.axes4, True)

    # Turn hold on for all axes
    self.axes1.hold(True)
    self.axes2.hold(True)
    self.axes3.hold(True)
    self.axes4.hold(True)

    # Get the s parameters out 
    s11 = [i[0][0] for i in s]
    s12 = [i[0][1] for i in s]
    s21 = [i[1][0] for i in s]
    s22 = [i[1][1] for i in s]   

    # Calculate the angle and phase for sij
    theta,r = [angle(s) for s in s11],[abs(s) for s in s11]
    h1 = self.axes1.plot(theta,r,linewidth=3)
    # Show freq limits
    self.axes1.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
    self.axes1.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)

    theta,r = [angle(s) for s in s12],[abs(s) for s in s12]
    h1 = self.axes2.plot(theta,r,linewidth=3)
    # Show freq limits
    self.axes2.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
    self.axes2.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)


    theta,r = [angle(s) for s in s21],[abs(s) for s in s21]
    h1 = self.axes3.plot(theta,r,linewidth=3)
    # Show freq limits
    self.axes3.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
    self.axes3.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)
 
    theta,r = [angle(s) for s in s22],[abs(s) for s in s22]
    h1 = self.axes4.plot(theta,r,linewidth=3)
    # Show freq limits
    self.axes4.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
    self.axes4.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)
    
    # Set up for working legend
    self.hlist.append(h1)
    self.tlist.append(string)
    self.axes = self.axes2

  def plotInput(self, s, string="data"): 

    # Make a larger window 
    self.win.set_default_size(1300,500)

    # Set up the axes and subplots s11,s12 
    self.axes1 = self.figure.add_subplot(1,2,1,projection='polar')
    self.plotsmithcircles(self.axes1,True)
       
    self.axes2 = self.figure.add_subplot(1,2,2,projection='polar')
    
    # Turn hold on for all axes
    self.axes1.hold(True)
    self.axes2.hold(True)

    # Get the s parameters out 
    s11 = [i[0][0] for i in s]
    s12 = [i[0][1] for i in s]

    # Calculate the angle and phase for sij
    theta,r = [angle(s) for s in s11],[abs(s) for s in s11]
    h1 = self.axes1.plot(theta,r,linewidth=3)
    # Show freq limits
    self.axes1.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
    self.axes1.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)

    theta,r = [angle(s) for s in s12],[abs(s) for s in s12]
    h1 = self.axes2.plot(theta,r,linewidth=3)
    # Show freq limits
    self.axes2.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
    self.axes2.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)
    
    # Set up for working legend
    self.hlist.append(h1)
    self.tlist.append(string)
    self.axes = self.axes2

  def plotOutput(self,s,string="data"):

    # Make a larger window 
    self.win.set_default_size(1300,500)

    # Set up the axes and subplots s11,s12 
    self.axes1 = self.figure.add_subplot(2,1,1,projection='polar')
    #self.plotsmithcircles(self.axes1,True)
    
    self.axes2 = self.figure.add_subplot(2,1,2,projection='polar')
     
    # Turn hold on for all axes
    self.axes1.hold(True)
    self.axes2.hold(True)

    # Get the s parameters out 
    s21 = [i[1][0] for i in s]
    s22 = [i[1][1] for i in s]

    # Calculate the angle and phase for sij
    theta,r = [angle(s) for s in s21],[abs(s) for s in s11]
    h1 = self.axes1.plot(theta,r,linewidth=3)
    # Show freq limits
    self.axes1.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
    self.axes1.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)

    theta,r = [angle(s) for s in s22],[abs(s) for s in s12]
    h1 = self.axes2.plot(theta,r,linewidth=3)
    # Show freq limits
    self.axes2.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
    self.axes2.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)

    # Set up for working legend
    self.hlist.append(h1)
    self.tlist.append(string)
    self.axes = self.axes2

  def plotsmithcircles(self,axes, reduced=False):

    # Disable the polar angle marks
    axes.set_rgrids([1],['',''])
    axes.set_thetagrids([0,180],['',''])
    
    # Define linetypes
    _typemajor,_typeminor,_minor, _major = 'k','k:',0.5,1

    # plot R circles  
    if reduced:
      r_array = (0.0,0.2,0.5,1.0,2.0,5.0)
      x_array = (0.2,0.5,1.0,2.0,5.0)
      val = 1.15

    else:
      r_array = (0.0,0.2,0.4,0.7,1.0,1.4,2.0,3.,5.0)
      x_array = (0.2,0.4,0.6,0.8,1.0,1.4,2.0,3.0,5.0,10.0)
      val = 1.08
    
    z = zeros((1,2000),complex)
    z.imag = arange(-50,50,0.05) 
    x_shift = 1e-3 # can be 1e-3 but not zero
    for r in r_array:
      z.real = r
      gamma_r = transpose((z-1)/(z+1))
      
      if r==1:
        axes.plot(angle(gamma_r), abs(gamma_r), _typemajor,linewidth = float(_major) ,zorder=0)
      else:
        axes.plot(angle(gamma_r), abs(gamma_r), _typeminor,linewidth = float(_minor),zorder=0)
      
      if not r==0:
        axes.text(angle((r+x_shift*1j-1)/(r+x_shift*1j+1)),1e-6+abs((r+x_shift*1j-1)/(r+x_shift*1j+1)),
                       "$%.1f$" %(r), backgroundcolor='white',horizontalalignment='center',
                        verticalalignment='center',zorder=1, fontsize=16)      
            
    # plot inductive circles (top half plane)
    z = zeros((1,1000),complex)
    z.real = arange(0, 100, 0.1)
    for x in x_array:
      z.imag = x
      gamma_x = transpose((z-1)/(z+1))
      
      if x==1:
        axes.plot(angle(gamma_x),abs(gamma_x),_typemajor,linewidth=float(_major), zorder=1)
        axes.text(angle((x*1j-1)/(x*1j+1)),val,r'$%.1fi$'%x,horizontalalignment='center',
                  verticalalignment='center',fontsize=16)
      else:
        axes.plot(angle(gamma_x),abs(gamma_x),_typeminor,linewidth=float(_minor),zorder=1)
        axes.text(angle((x*1j-1)/(x*1j+1)),val,r'$%.1fi$'%x,horizontalalignment='center',
                  verticalalignment='center',fontsize=16)
      
    # plot capacitive circles (bottom half plane)
    for x in multiply(-1,x_array):
      z.imag = x
      gamma_x = transpose((z-1)/(z+1))
      
      if x==-1:
        axes.plot(angle(gamma_x),abs(gamma_x),_typemajor,linewidth=float(_major),zorder=1)
        axes.text(angle((x*1j-1)/(x*1j+1)),val,r'$%.1fi$'%x,horizontalalignment='center',
                  verticalalignment='center',fontsize=16)

      else:
        axes.plot(angle(gamma_x),abs(gamma_x),_typeminor,linewidth=float(_minor),zorder=1)
        axes.text(angle((x*1j-1)/(x*1j+1)),val,r'$%.1fi$'%x,horizontalalignment='center',
                     verticalalignment='center',fontsize=16)
 
    # plot x=0 line
    axes.plot([-pi, 0],[1 , 1], _typemajor,linewidth=0.5,zorder=1)
 
    # plot zero and infinity symbols
    axes.text(0 ,val,'$\infty$',horizontalalignment='center',verticalalignment='center',fontsize=16)
    axes.text(-pi ,val,'$0$',horizontalalignment='center',verticalalignment='center',fontsize=16)

def main():

  freq = logspace(0,9,10000)
  C1,C2,L  = 1e-9,1e-9,5e-6

  y12 = lambda w: complex(0,-w*C1+(1/(w*L)))
  y11 = lambda w: complex(.01,-(1/(w*L))+w*C1)
  

  s11_vec,s12_vec = [],[]
  for f in freq:
    y = array([[y11(f),y12(f)],[y12(f),y11(f)]])
    s11_vec.append(mwc.ytos(y)[0][0])
    s12_vec.append(mwc.ytos(y)[0][1])

  #print [angle(s) for s in s11_vec]
  #print [abs(s) for s in s11_vec]
    
  plt = plotSmith()
  plt.plot(s11_vec)
  plt.plot(s12_vec)

  plt.show()
    
if __name__ == '__main__':
  main()
  

    



