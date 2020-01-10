#!/usr/bin/env python 
import numpy as np
import math
import sys

# Import QT backends
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QComboBox, QCheckBox, QLabel, QMessageBox, QSizePolicy, QApplication
from PyQt5.QtGui import QIcon

# Import matplotlib QT backends
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt

# Import mwconverter
import Converter as mwc

# Object for plotting smith charts
class plotSmith:

	def __init__(self):

		# Create figure
		self.figure  = plt.figure(figsize=(6,6))

		# Add polar axes
		self.axes = self.figure.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)

		# Define position for text labels
		self.params = {
			'style_major'	: 'k',
			'style_minor'	: 'k:',
			'lw_major'		: 0.7,
			'lw_minor'		: 0.5,		
			'text_position'	: 1.08,
			'grid_points'	: 2048
		}
 
		# Draw the Smith Chart circles
		self.plotSmithCircles(self.axes)


	def plotSmithCircles(self, axes, reduced=False):

		# Disable the polar np.angle marks
		axes.set_rlim(0,1)
		axes.set_rgrids([1],['',''])
		axes.set_thetagrids([0,180],['',''])


		# Constant resistance circles
		for r in [0.0,0.2,0.5,1.0,2.0,5.0]:

			# Create an array of zeros		
			z = np.zeros(self.params['grid_points'], complex)

			# Synthesize impedance
			z.real = r
			z.imag = np.linspace(-50.0,50.0, self. params['grid_points']) 

			# Calculate reflection coefficient
			gamma = np.transpose((z-1)/(z+1))
			
			# Plot major gridline circles
			if r == 1.0:
				axes.plot(
					np.angle(gamma), np.abs(gamma), 
					self.params['style_major'],
					lw = self.params['lw_major'], 
					zorder=0
				)

			# Plot minor gridline circles			
			else:
				axes.plot(
					np.angle(gamma), np.abs(gamma), 
					self.params['style_minor'],
					lw = self.params['lw_minor'], 
					zorder = 0
				)

				
		# Constant reactance circles
		for x in [0.2,0.5,1.0,2.0,5.0]:

			# Create an array of zeros		
			z = np.zeros(self.params['grid_points'], complex)

			# Synthesize impedance 
			z.real = np.linspace(0, 100.0, self. params['grid_points']) 
		
			# Need constant capacitive(x) and inductive (-x) circles
			for x in [x, -1.0*x]:

				z.imag = x

				# Calculate reflection coefficient
				gamma = np.transpose((z-1)/(z+1))
				
				# Plot major gridline circles
				if (x == 1.0) or (x == -1.0):
					axes.plot(
						np.angle(gamma), np.abs(gamma), 
						self.params['style_major'],
						lw = self.params['lw_major'], 
						zorder=0
					)

				# Plot minor gridline circles			
				else:
					axes.plot(
						np.angle(gamma), np.abs(gamma), 
						self.params['style_minor'],
						lw = self.params['lw_minor'], 
						zorder = 0
					)

				# Add text labels	
				axes.text(
					np.angle((x*1j-1)/(x*1j+1)),
					self.params['text_position'],
					r'$%.1fi$'%x,
					horizontalalignment='center',
					verticalalignment='center'
				)	


		# Add infinity label (open circuit)
		axes.text(0, self.params['text_position'],
			'$\infty$',
			horizontalalignment='center',
			verticalalignment='center'
		)


		# Add zero label (short circuit)		
		axes.text(-math.pi,self.params['text_position'],
			'$0$',
			horizontalalignment='center',
			verticalalignment='center'
		)
	
	# Wrap plot.plot()
	def plot(self, *args, **kwargs):

		# If we only pass one argument (array of complex numbers) 
		if len(args) == 1.0:
			self.axes.plot(np.angle(args[0]), np.abs(args[0]), **kwargs)

		# Otherwise simply pass along the parameters
		else:
			self.axes.plot(*args, **kwargs)


	# Wrap plot.show()
	def show(self):
		plt.show()	


if __name__ == '__main__':

	# Create a test network
	freq = np.logspace(0,8,1000)

	# Calculate a test network
	C1,C2,L	= 1e-9, 1e-9, 5e-6
	y12 = lambda w: complex(0,-w*C1+(1/(w*L)))
	y11 = lambda w: complex(.01,-(1/(w*L))+w*C1)
	
	
	s11_vec,s12_vec = [],[]
	for f in freq:
		y = np.array([[y11(f),y12(f)],[y12(f),y11(f)]])
		s11_vec.append(mwc.ytos(y)[0][0])
		s12_vec.append(mwc.ytos(y)[0][1])
	
	_plt = plotSmith()
	_plt.plot(s11_vec)
	_plt.plot(np.angle(s12_vec), np.abs(s12_vec))
	_plt.show()



	# def set_legend(self):
	# 	# Note that it is outside the axes. 
	# 	self.axes.legend(self.hlist,self.tlist, bbox_to_anchor=(1.05, 1), loc=2)
		
	# def plotSmith(self,s_data,string="data",_ls=None,c=None,lw=None):
	# 	# Make a larger window 
	# 	self.win.set_default_size(800,800)
		
	# 	if _ls is None:
	# 		_ls = '-'
	# 	if c is None:
	# 		c='k'
	# 	if lw is None:
	# 		lw=2
			
	# 	# Generate an axes object to plot into 
	# 	self.axes = self.figure.add_axes([0.1, 0.1, 0.8, 0.8],polar=True)
 
	# 	# Draw the Smith Chart circles
	# 	self.plotsmithcircles(self.axes)
	# 	self.axes.hold(True)	

	# 	# Plot the data and set r = 1
	# 	theta,r = [np.angle(s) for s in s_data],[np.abs(s) for s in s_data]
	# 	h1 = self.axes.plot(theta,r,_ls, color=c, linewidth=float(lw))
		
	# 	# Append to list for legend
	# 	self.hlist.append(h1)
	# 	self.tlist.append(string)

	# 	# Show freq limits
	# 	#self.axes.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
	# 	#self.axes.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)
	# 	self.axes.set_rmax(1)	 
		 
	# def plotPolar(self, s_data, string="data"):
	# 	# Generate an axes object to plot into 
	# 	self.axes = self.figure.add_axes([0.0, 0.0, 1.0, 1.0],polar=True)
	# 	self.axes.hold(True)
		
	# 	# Plot the data and set r = 1
	# 	theta,r = [np.angle(s) for s in s_data],[np.abs(s) for s in s_data]
	# 	h1 = self.axes.plot(theta,r,linewidth=3)
		
	# 	# Append to list for legence
	# 	self.hlist.append(h1)
	# 	self.tlist.append(string)
	 
	# 	# Show freq limits
	# 	self.axes.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
	# 	self.axes.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)

	# def plotTwoport(self, s, string="data"):

	# 	# Make a larger window 
	# 	self.win.set_default_size(1300,1000)

	# 	# Set up the axes and subplots
	# 	# s11,s22 in smith chart
	# 	# s12,s21 in polar chart
	# 	self.axes1 = self.figure.add_subplot(2,2,1,projection='polar')
	# 	self.plotsmithcircles(self.axes1,True)
		
	# 	self.axes2 = self.figure.add_subplot(2,2,2,projection='polar')
	# 	self.axes3 = self.figure.add_subplot(2,2,3,projection='polar')
		
	# 	self.axes4 = self.figure.add_subplot(2,2,4,projection='polar')
	# 	self.plotsmithcircles(self.axes4, True)

	# 	# Turn hold on for all axes
	# 	self.axes1.hold(True)
	# 	self.axes2.hold(True)
	# 	self.axes3.hold(True)
	# 	self.axes4.hold(True)

	# 	# Get the s parameters out 
	# 	s11 = [i[0][0] for i in s]
	# 	s12 = [i[0][1] for i in s]
	# 	s21 = [i[1][0] for i in s]
	# 	s22 = [i[1][1] for i in s]	 

	# 	# Calculate the np.angle and phase for sij
	# 	theta,r = [np.angle(s) for s in s11],[np.abs(s) for s in s11]
	# 	h1 = self.axes1.plot(theta,r,linewidth=3)
	# 	# Show freq limits
	# 	self.axes1.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
	# 	self.axes1.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)

	# 	theta,r = [np.angle(s) for s in s12],[np.abs(s) for s in s12]
	# 	h1 = self.axes2.plot(theta,r,linewidth=3)
	# 	# Show freq limits
	# 	self.axes2.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
	# 	self.axes2.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)


	# 	theta,r = [np.angle(s) for s in s21],[np.abs(s) for s in s21]
	# 	h1 = self.axes3.plot(theta,r,linewidth=3)
	# 	# Show freq limits
	# 	self.axes3.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
	# 	self.axes3.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)
 
	# 	theta,r = [np.angle(s) for s in s22],[np.abs(s) for s in s22]
	# 	h1 = self.axes4.plot(theta,r,linewidth=3)
	# 	# Show freq limits
	# 	self.axes4.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
	# 	self.axes4.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)
		
	# 	# Set up for working legend
	# 	self.hlist.append(h1)
	# 	self.tlist.append(string)
	# 	self.axes = self.axes2

	# def plotInput(self, s, string="data"): 

	# 	# Make a larger window 
	# 	self.win.set_default_size(1300,500)

	# 	# Set up the axes and subplots s11,s12 
	# 	self.axes1 = self.figure.add_subplot(1,2,1,projection='polar')
	# 	self.plotsmithcircles(self.axes1,True)
			 
	# 	self.axes2 = self.figure.add_subplot(1,2,2,projection='polar')
		
	# 	# Turn hold on for all axes
	# 	self.axes1.hold(True)
	# 	self.axes2.hold(True)

	# 	# Get the s parameters out 
	# 	s11 = [i[0][0] for i in s]
	# 	s12 = [i[0][1] for i in s]

	# 	# Calculate the np.angle and phase for sij
	# 	theta,r = [np.angle(s) for s in s11],[np.abs(s) for s in s11]
	# 	h1 = self.axes1.plot(theta,r,linewidth=3)
	# 	# Show freq limits
	# 	self.axes1.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
	# 	self.axes1.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)

	# 	theta,r = [np.angle(s) for s in s12],[np.abs(s) for s in s12]
	# 	h1 = self.axes2.plot(theta,r,linewidth=3)
	# 	# Show freq limits
	# 	self.axes2.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
	# 	self.axes2.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)
		
	# 	# Set up for working legend
	# 	self.hlist.append(h1)
	# 	self.tlist.append(string)
	# 	self.axes = self.axes2

	# def plotOutput(self,s,string="data"):

	# 	# Make a larger window 
	# 	self.win.set_default_size(1300,500)

	# 	# Set up the axes and subplots s11,s12 
	# 	self.axes1 = self.figure.add_subplot(2,1,1,projection='polar')
	# 	#self.plotsmithcircles(self.axes1,True)
		
	# 	self.axes2 = self.figure.add_subplot(2,1,2,projection='polar')
		 
	# 	# Turn hold on for all axes
	# 	self.axes1.hold(True)
	# 	self.axes2.hold(True)

	# 	# Get the s parameters out 
	# 	s21 = [i[1][0] for i in s]
	# 	s22 = [i[1][1] for i in s]

	# 	# Calculate the np.angle and phase for sij
	# 	theta,r = [np.angle(s) for s in s21],[np.abs(s) for s in s11]
	# 	h1 = self.axes1.plot(theta,r,linewidth=3)
	# 	# Show freq limits
	# 	self.axes1.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
	# 	self.axes1.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)

	# 	theta,r = [np.angle(s) for s in s22],[np.abs(s) for s in s12]
	# 	h1 = self.axes2.plot(theta,r,linewidth=3)
	# 	# Show freq limits
	# 	self.axes2.plot(theta[0],r[0],'oc',linewidth = 5,zorder=10)
	# 	self.axes2.plot(theta[-1],r[-1],'sm',linewidth = 5,zorder=10)

	# 	# Set up for working legend
	# 	self.hlist.append(h1)
	# 	self.tlist.append(string)
	# 	self.axes = self.axes2

	# 

		
#	main()
	

		



