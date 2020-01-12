# ---------------------------------------------------------------------------------
# 	minispice-> plotAnalysis.py
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
import math
import sys

# Import mwconverter
from .Converter import *

# Import matplotlib 
import matplotlib.pyplot as plt

# Helper class to plot frequency analysis data
class plotAnalysis: 

	def __init__(self):

		# Dictionary of figures
		self.figures = {}

	# Method to add a standard figure
	def add_figure(self, key):
		
		# All plots will have the following
		fig = plt.figure()

		# Create some axes
		ax = fig.add_subplot(111)
		ax.grid(True, which="both",ls="-", color='0.65')

		# Add figure to dict
		self.figures[key] = fig

	# Method to add polar plot figures
	def add_polar(self, key):

		fig = plt.figure(figsize=(6,6))

		# Add polar axes
		ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)

		# Add figure to dict
		self.figures[key] = fig


	# Method to add smith chart figures	
	def add_smith(self, key):

		fig = plt.figure(figsize=(6,6))

		# Add polar axes
		ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)

		# Draw the Smith Chart circles
		self.buildSmithChart(ax)

		# Add figure to dict
		self.figures[key] = fig

	# Methods to set axes lables and title
	def set_xlabel(self, key, x_label):
		self.figures[key].axes[0].set_xlabel(x_label)

	def set_ylabel(self, key, y_label):
		self.figures[key].axes[0].set_ylabel(y_label)

	def set_title(self, key, title):
		self.figures[key].axes[0].set_title(title)

	# Plotting method (linear plots)
	def plot(self, key, x_data, y_data, arg, **kwargs):
		
		# Check for valid plotting modes
		try:

			# Linear frequency 
			if arg == "lin":		
				h = self.figures[key].axes[0].plot(x_data, y_data, **kwargs)
				
			elif arg == "lin(dB)":
				h = self.figures[key].axes[0].plot(x_data, [todB(_) for _ in y_data], **kwargs)
				
			# Log frequency 
			elif arg == "log":
				h = self.figures[key].axes[0].semilogx(x_data, y_data, **kwargs)
				
			elif arg == "log(dB)":
				h = self.figures[key].axes[0].semilogx(x_data, [todB(_) for _ in y_data], **kwargs)

			# Double logarithmic
			elif arg == "loglog":
				h = self.figures[key].axes[0].loglog(x_data, y_data, **kwargs)

			else: 
				raise ValueError

			# Return axes handle
			return h
	
		# Invlaid plotting mode
		except ValueError:
			print("Invalid plot configuration (%s)"%arg)
			print("\t Select from [\"lin\", \"lin(dB)\", \"log\", \"log(dB)\", \"loglog\"]")	

	# Plotting method (polar charts)
	def plot_polar(self, key, data, **kwargs):

		# We only pass one argument (array of complex numbers) 
		h = self.figures[key].axes[0].plot(np.angle(data), np.abs(data), **kwargs)

		# Return artist object
		return h	

	# Plotting method (smith charts)
	def plot_smith(self, key, data, **kwargs):

		# We only pass one argument (array of complex numbers) 
		h = self.figures[key].axes[0].plot(np.angle(data), np.abs(data), **kwargs)

		# Return artist object
		return h

	# Build Smith Chart constant resistance and reactance lines
	def buildSmithChart(self, ax):
	
		# Define position for text labels
		params = {
			'style_major'	: 'k',
			'style_minor'	: 'k:',
			'lw_major'		: 0.7,
			'lw_minor'		: 0.5,		
			'text_position'	: 1.08,
			'grid_points'	: 2048
		}
 	
		# Disable the polar np.angle marks
		ax.set_rlim(0,1)
		ax.set_rgrids([1],['',''])
		ax.set_thetagrids([0,180],['',''])

		# Constant resistance circles
		for r in [0.0, 0.2, 0.5, 1.0, 2.0, 5.0]:

			# Create an array of zeros		
			z = np.zeros(params['grid_points'], complex)

			# Synthesize impedance
			z.real = r
			z.imag = np.linspace(-50.0,50.0, params['grid_points']) 

			# Calculate reflection coefficient
			gamma = np.transpose((z-1)/(z+1))
			
			# Plot major gridline circles
			if r == 1.0:
				ax.plot(
					np.angle(gamma), np.abs(gamma), 
					params['style_major'],
					lw = params['lw_major'], 
					zorder=0
				)

			# Plot minor gridline circles			
			else:
				ax.plot(
					np.angle(gamma), np.abs(gamma), 
					params['style_minor'],
					lw = params['lw_minor'], 
					zorder = 0
				)

				
		# Constant reactance circles
		for x in [0.2, 0.5, 1.0, 2.0, 5.0]:

			# Create an array of zeros		
			z = np.zeros(params['grid_points'], complex)

			# Synthesize impedance 
			z.real = np.linspace(0, 100.0, params['grid_points']) 
		
			# Need constant capacitive(x) and inductive (-x) circles
			for x in [x, -1.0*x]:

				z.imag = x

				# Calculate reflection coefficient
				gamma = np.transpose((z-1)/(z+1))
				
				# Plot major gridline circles
				if (x == 1.0) or (x == -1.0):
					ax.plot(
						np.angle(gamma), np.abs(gamma), 
						params['style_major'],
						lw = params['lw_major'], 
						zorder=0
					)

				# Plot minor gridline circles			
				else:
					ax.plot(
						np.angle(gamma), np.abs(gamma), 
						params['style_minor'],
						lw = params['lw_minor'], 
						zorder = 0
					)

				# Add text labels	
				ax.text(
					np.angle((x*1j-1)/(x*1j+1)),
					params['text_position'],
					r'$%.1fi$'%x,
					horizontalalignment='center',
					verticalalignment='center'
				)	

		# Add infinity label (open circuit)
		ax.text(0, params['text_position'],
			'$\infty$',
			horizontalalignment='center',
			verticalalignment='center'
		)

		# Add zero label (short circuit)		
		ax.text(-math.pi,params['text_position'],
			'$0$',
			horizontalalignment='center',
			verticalalignment='center'
		)	

	# Wrap mpl show plots
	def show(self):
		plt.show()
