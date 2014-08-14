import math
import numpy as np
import pylab
import sys
from scipy import interpolate
import random


'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Artificial Lightcurve generation functions - parameters 
control the curve being generated. Can be used to imprint 
Seasonal and observational effects.
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
'''
#--------------------------------------------------------------------
# List of parameters
#--------------------------------------------------------------------

tsteps = 1024
mass = 1.0E8
radius = 100
sigma = 1.0
b = 0.
alpha = 0.01

#Imprint effects

#Starting date of observations
mhjd = 50000

#Maximum shift between the two curves. A shift will be selected at random between range 0 to Max_Shift
#MAXIMUM SHIFT SHOULD BE LESS THAN THE T_STEPS
Max_Shift = 750

#Difference in days between two observations
Days_diff = 1

#Seasonal difference in days between observations
Seasonal_diff = 90

#Probability of a seasonal effect being imprinted
Prob_Seasonal = 0.01

#Probability of a day being non- observable
Prob_Day = 0.10

#Curve Offsets (values added to returned intensities)
Offset1 = 0
Offset2 = 0

#Mean and variance in observational errors of light-curves
Mean_error_1 = 0.00418616252822 
Mean_error_2 = 0.00543952595937
Variance_1 	 = 0.000107790282228 
Variance_2	 = 0.000155465436713

#--------------------------------------------------------------------
#--------------------------------------------------------------------


class l_curve:

	'''
	A Class for generating lightcurves. Intended use a "pair" of light curves. Class structure might be changed soon.
	'''
	def stochasticIntegral(self, tau):
		steps = 5 # Default number of sampling points for integral evaluation
	
		# Normal parameters
		mean = 0.0
		sigma = 1.0
	
		# Our integral sample points
		sample = np.arange(steps)
	
		# Our integration is simply summing up the random dt and the function value

		values = np.exp(-sample/tau) * np.random.normal(mean, sigma, steps)

		return (np.sum(values))

	def initTimescales(self, mass, radius, alpha):
		# Renormalize to base units
		mass /= 1.0E8
		radius /= 100.0

		# The B. Kelly model has three timescales
		tau = np.zeros(3)
	
		# First we compute the light-crossing timescale
		tau[0] = 1.1 * mass * radius
	
		# Second, we compute the orbital timescale
		tau[1] = 104.0 * mass * (radius**1.5)
	
		# Third we calculate the thermal timescale
		diny = 365 # days in a year.
		tau[2] = 4.6 * (0.01/alpha) * mass * (radius**1.5) * diny

		return tau


	def __init__(self,mass, radius, alpha, b, sigma, tsteps):

		# Define parameters:
		ntimes = 3 # Number of timescales
		
		# First get timescales
		tau = self.initTimescales(mass, radius, alpha)
		
		# Initialize arrays
		fluxes  = np.zeros((ntimes,tsteps))
		ssi = np.zeros(ntimes)
		
		# Initial random flux values
		fluxes[:,0] = np.random.normal(b * tau, sigma * (tau/2.0)**0.5)

		# Common calculations    
		bt = b * tau
		exp_tau = np.exp(-5.0/tau)

		for i in range(1, tsteps):
		    # Do our three stochastic integrals
		    for j in range(ntimes):
		        ssi[j] = sigma * self.stochasticIntegral(tau[j])
		    #added abs
		    fluxes[:,i] = exp_tau * fluxes[:, i - 1] + bt * (1.0 - exp_tau) + ssi
	 
		self.Curve = fluxes

	
	def getCurve(self):
		return self.Curve

	def getThermalCurve(self):
		return self.Curve[2]

	def getOrbitalCurve(self):
		return self.Curve[1]

	def getLightCrossCurve(self):
		return self.Curve[0]

	def display(self):
		print self.Curve

	def getError(self,Mean,Variance): 
		error = []		
		for i in range(0, tsteps):
			error.append(random.gauss(Mean,Variance))
		return error	

	def getTimes(self,start): 
		time = []		
		for i in range(0, tsteps):
			time.append(start + i)
		return time

	def addSeasonalEffects(self, Flux1, error1, time1, Flux2, error2, time2, tsteps):
		X = []
		shift = abs(time2[0] - time1[0])
		#Trim The curves 		
		for i in range(0,int(shift)):
			del Flux1[0]
		for i in range(0,int(tsteps)):		
			del Flux2[-1]
		for i in range(0,int(tsteps-shift)):
			del Flux1[-1]
		
		#Append points	
		i=0
		while (i < tsteps):		
		
			p = random.random()
			if(p <= Prob_Seasonal):
				#Remove points
				i = i + Seasonal_diff 
			else:
				#Append points
				X.append([time2[i],Flux1[i],error1[i],Flux2[i], error2[i]])
				i = i + 1
		return X


	def addDailyEffects(self, Light_data):
		i=0		
		while (i < len(Light_data)):		
			p = random.random()
			if(p <= Prob_Day):
				del Light_data[i]
			i = i+1
		return Light_data
	
'''
Plots the data points 
'''
def plot_data(D):
	t  = []
	c1 = []
	c2 = []
	for i in range(0,len(D)):		
		t.append(D[i][0])
		c1.append(D[i][1])
		c2.append(D[i][3])
	pylab.plot(t, c1, 'rx', label='Curve A')
	pylab.plot(t, c2, 'bx', label='Curve B')
	pylab.show()


def generateFile(filename, Data ):

	#Write to file
	output = open(filename, 'w')
	output.write('mhjd\tmag_A\tmagerr_A\tmag_B\tmagerr_B\ttelescope\n')
	output.write('====\t=====\t========\t=====\t========\t=========\n')
	
	for i in range (0,len(Data)):
		output.write(str(D[i][0])+'\t')
		output.write(str(D[i][1])+'\t')
		output.write(str(D[i][2])+'\t')
		output.write(str(D[i][3])+'\t')
		output.write(str(D[i][4])+'\t')
		output.write('None'+'\n')

	output.close()


# Obsolete : TODO: DELETE THIS!

def function1(mass, radius, alpha, b, sigma, tsteps, mhjd):
	flux1 = l_curve(mass, radius, alpha, b, sigma, tsteps)
	X1 = flux1.getCurve()

	flux2 = l_curve(mass, radius, alpha, b, sigma, tsteps)
	X2 = flux2.getCurve()

	
	#Write to file

	output = open('Artificial_LC.rdb', 'w')
	output.write('mhjd\tmag_A\tmagerr_A\tmag_B\tmagerr_B\ttelescope\n')
	output.write('====\t=====\t========\t=====\t========\t=========\n')
	

	lc1 = []
	lc2 = []
	time= []
	for i in range (0,len(X2[0])):
		output.write(str(mhjd)+'\t')
		time.append(mhjd)
		mhjd = mhjd + random.expovariate(1/Mean_diff)
		p = random.random()
		if(p>.95):
			mhjd = mhjd + random.expovariate(1/Seasonal_diff)
		output.write(str(Offset1+X1[0][i])+'\t')
		lc1.append(X1[0][i])
		output.write(str(random.gauss(Mean_error_1,Variance_1)) +'\t')
		output.write(str(Offset2+X2[0][i])+'\t')
		lc2.append(X2[0][i])
		output.write(str(random.gauss(Mean_error_2,Variance_2))+'\t')
		output.write('None'+'\n')


	output.close()

	Smooth1 = interpolate.UnivariateSpline(time, lc1)(time)
	Smooth2 = interpolate.UnivariateSpline(time, lc2)(time)	
	pylab.plot(time, lc1,'bx')
	pylab.plot(time, Smooth1,'r')
	pylab.plot(time, lc2,'rx')
	pylab.plot(time, Smooth2,'b')
	pylab.show()


if '__main__':
	
	#function1(mass, radius, alpha, b, sigma, tsteps, mhjd)

	#Generate a light curve twice of twice the light steps - to allow a margin for shift
	LC = l_curve(mass, radius, alpha, b, sigma, 2*tsteps)

	Flux1 = LC.getThermalCurve()
	#Flux2 = LC.getThermalCurve()
	Flux2 = Flux1
	error1 = LC.getError(Mean_error_1, Variance_1)
	error2 = LC.getError(Mean_error_2, Variance_2)
	time1  = LC.getTimes(mhjd)
	tshift = int(random.uniform(0,Max_Shift))
	print tshift
	time2  = LC.getTimes(mhjd+tshift)
	
	# time2 must be shifted ahead of time1
	S = LC.addSeasonalEffects(Flux1.tolist(), error1, time1, Flux2.tolist(), error2, time2, tsteps)
	D = LC.addDailyEffects(S)
	plot_data(S)
	generateFile('Artificial_LC_'+str(tshift)+'.rdb',S)
	
