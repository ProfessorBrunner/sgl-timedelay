import sys
import pylab
from scipy import interpolate
import random
from dcor import * 
import operator
import numpy as np
'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Has functions to load, segment, plot and estimate the shift
between two light curves. Takes an input Light Curve file 
and returns the time shift.
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
'''
'''
Class for light curve segment
'''
class Segment(object):
	def __init__(self, intensity_function, intensity_data, time_data):
		self.intensity_function = intensity_function
		self.time_start  = min(time_data)
		self.time_end 	 = max(time_data)
		self.time_data   = time_data
		self.intensity	 = intensity_data

	'''
	Shifts a light curve segment
	'''	
	def shift_segment(self, days):
		self.time_start  = self.time_start + days
		self.time_end 	 = self.time_end + days
		
		for x in range(0,len(self.time_data)):
			self.time_data[x] = self.time_data[x] + days

		self.intensity_function = interpolate.UnivariateSpline(self.time_data, self.intensity,k=3, s=0.04)

	'''
	Returns the intensities corresponding to temporal overlap of light curves
	'''
	def intersection(self, segment):
		#if no overlap return empty
		if ( (self.time_end <= segment.getStartTime()) or (self.time_start >= segment.getEndTime()) ) :	
			return [] , []
		else: 
			start, end = -1, -1
			#find the common region			
			if(segment.getEndTime()<= self.time_end) and (self.time_start >= segment.getStartTime()):
				start	= self.time_start
				end		= segment.getEndTime()
			elif(segment.getEndTime()<= self.time_end) and (self.time_start <= segment.getStartTime()):
				start	= segment.getStartTime()
				end		= segment.getEndTime()
			elif(segment.getEndTime()>= self.time_end) and (self.time_start <= segment.getStartTime()):
				start	= segment.getStartTime()
				end		= self.time_end
			elif(segment.getEndTime()>= self.time_end) and (self.time_start >= segment.getStartTime()):
				start	= self.time_start
				end		= self.time_end
			
			common_segment1 = []
			common_segment2 = []	
			for t in range (int(start), int(end)):	
				common_segment1.append(self.intensity_function(t))
				common_segment2.append(segment.getIntensity(t))
			return common_segment1 , common_segment2

	'''
	Returns the intensitiy at a given time
	'''
	def getIntensity(self, time):
		return self.intensity_function(time)
	
	'''
	Returns the start time
	'''
	
	def getStartTime(self):
		return self.time_start

	'''
	Returns the end time
	'''
	
	def getEndTime(self):
		return self.time_end
	
	'''
	Returns the time across which a segment spans
	'''
	
	def getSegmentLength(self):
		return self.time_end - self.time_start

	'''
	Prints a segments parameters
	'''
	
	def display(self):
		print 'Segment start:'
		print self.time_start
		print 'Segment end:'
		print self.time_end

'''
Reads a file and returns a lists of elements having the format
[ Time, [[Lightcurve Ai, Error in Ai] for i curves]]
Header of input file is ignored (first two lines) 
'''
def load_data():
	if (len(sys.argv)!=2):
		print 'Please specify one input file'
		exit(0)
	time = []
	Lcurves = []
	infile = open(sys.argv[1], 'r')
	lines = infile.readlines()
	Data = []	
	for i in range (2,len(lines)):
		elements = lines[i].split()
		time.append(float(elements[0]))
		temp = []
		for x in range(1, int(len(elements)/2)):		
			temp.append([float(elements[2*x-1]),float(elements[2*x])])		
		Lcurves.append(temp)	
	infile.close()
	return time, Lcurves

'''
Plots the data points 
'''
def plot_data_points(time, lc_a, lc_b, avg_a, avg_b):
	pylab.figure()	
	pylab.plot(time, [l + avg_a for l in lc_a], 'rx' , label='Curve A')
	pylab.plot(time, [l + avg_b for l in lc_b], 'bx', label='Curve B')


'''
Plots a Light Curve Segment 
'''
def plot_segment(LC, avg):
	t_range = []
	intensity = []	
	for t in range(int(LC.getStartTime()),int(LC.getEndTime())): 	
		t_range.append(t)
		intensity.append(avg + LC.getIntensity(t))		
	pylab.plot(t_range, intensity,'g')

'''
Show all the plots
'''
def show_plots():
	pylab.legend(loc=4)
	pylab.xlabel('Modified Heliocentric Julian Date (MHJD)')
	pylab.ylabel('Magnitude (relative)')
	pylab.title(sys.argv[1])
	pylab.show()

'''
Find shift between two Light-Curve Segments
'''
def find_shift(LC_1, LC_2, Max_Shift, Shift_Step):
	
	Corr_Data = []
	i = 0

	#Shift Light Curve 1		
	while (i < Max_Shift):
		i = i + Shift_Step		
		print 'Shift '+ str(i) + '...'		
		
		Intersection_vector_1 = []
		Intersection_vector_2 = []
		
		#Shift curve 1 forward		
		for s in LC_1:
			s.shift_segment(Shift_Step)

		#Find intersections
		for s1 in LC_1:		
			for s2 in LC_2:
				c1, c2 = s1.intersection(s2)
				if(c1 != [] and c2 != []):
					for x in range (0,len(c1)):
						Intersection_vector_1.append(c1[x])
						Intersection_vector_2.append(c2[x])	
		print 'Overlap Vector Length: ' + str(len(Intersection_vector_1))
		cVal = dcor(np.asarray(Intersection_vector_1),np.asarray(Intersection_vector_2)).find_correlation()
		Corr_Data.append([i, cVal])
		print 'Correlation: ' + str(cVal)
		
		 


	for s in LC_1:
		s.shift_segment(-1 * i)

		
	#Shift Light Curve 2
	i = 0
	while (i < Max_Shift):
		i = i + Shift_Step 		
		print 'Shift -'+ str(i) + '...'		
	
		#Shift curve 2 forward		
		Intersection_vector_1 = []
		Intersection_vector_2 = []
		
		for s in LC_2:
			s.shift_segment(Shift_Step)

		#Find intersections
		for s1 in LC_1:		
			for s2 in LC_2:
				c1, c2 = s1.intersection(s2)
				if(c1 != [] and c2 != []):
					for x in range (0,len(c1)):
						Intersection_vector_1.append(c1[x])
						Intersection_vector_2.append(c2[x])	
		print 'Overlap Vector Length: ' + str(len(Intersection_vector_1))
		cVal = dcor(np.asarray(Intersection_vector_1),np.asarray(Intersection_vector_2)).find_correlation()
		Corr_Data.append([-1 * i, cVal])
		print 'Correlation: ' + str(cVal)
		
	
	#Used only for plotting correctly
	for s in LC_2:
		s.shift_segment(-1 * i)
	
	
	# Find Maximum correlation Correlation 
	print Corr_Data
	Shift = max(Corr_Data, key=operator.itemgetter(1))
	return Shift

'''
Returns a list of light curve segments fitting the points of 2 light curves
'''
def fit_curve(time, lc_a, lc_b, MAX_GAP):

# Perform segmentation
	LC_1 = []
	LC_2 = []
	T_Segments = []
	i = 0
	while (i < (len(time)-1)):
		temp1 = []
		temp2 = []
		tempt = []		
		while (i < (len(time)-1)) and (time[i+1]-time[i]<MAX_GAP):
			temp1.append(lc_a[i])
			temp2.append(lc_b[i])
			tempt.append(time[i])
			i = i + 1
		if(len(temp1)>5):	# Each segment must have at least 5 observational points	

		# s = smoothing parameter, if s = 0 curve passes through all the points.        //k = 3, s = 0.04
		# interpolate.UnivariateSpline(TimeScale, Intensity, k {Degree of Spline}, s {Smoothning Factor} )
			LC_1.append( Segment(interpolate.UnivariateSpline(tempt, temp1, k=3, s=0.04), temp1, tempt) ) 
			LC_2.append( Segment(interpolate.UnivariateSpline(tempt, temp2, k=3, s=0.04), temp2, tempt) ) 
		
			T_Segments.append(tempt)
		i = i + 1

	return LC_1, LC_2, T_Segments

if '__main__':
	
	# Load Data from file	
	time, Lcurves = load_data()

	# Introduce error by mutiplying given error with a standard normal distribution
	
	for i in range (0, len(Lcurves)):
		for x in range(0, len(Lcurves[0]) ):		
			Lcurves[i][x][0] = Lcurves[i][x][0] + random.gauss(0,1) * Lcurves[i][x][1]
	
	
	# Compute mean of light curves and force it to zero	(non essential)
	avg = np.zeros(len(Lcurves[0]))
	for i in range (0, len(Lcurves)):
		for x in range(0, len(Lcurves[0]) ):		
			avg[x] = avg[x] + Lcurves[i][x][0]

	for x in range(0, len(Lcurves[0])):		
		avg[x] = avg [x] / len(Lcurves)
	
	for i in range (0, len(Lcurves)):
		for x in range(0, len(Lcurves[0]) ):		
			Lcurves[i][x][0] = Lcurves[i][x][0] - avg[x]
	

	# Perform Segmentation of Light Curve
	MAX_GAP = 60		# In Modified Heliocentric Julian Date	
	
	for x in range (0, len(Lcurves[0])-1):
		print 'Finding shift between curve ' + str (x+1) + ' and ' + str (x+2) 
		
		lc_1 = []
		lc_2 = []
	
		for i in range (0, len(Lcurves)):
			lc_1.append(Lcurves[i][x][0])
			lc_2.append(Lcurves[i][x+1][0])
	
		LC_1 = []
		LC_2 = []

		# Perform Segmentation and Spline fitting
		LC_1, LC_2, T_Segments = fit_curve(time, lc_1, lc_2, MAX_GAP)

		Max_Shift  = int((max(time)-min(time))/20)	
		Shift_Step = 1	

		# Find Distance Correlation 
		Shift = find_shift(LC_1, LC_2, Max_Shift, Shift_Step)

		print 'Observed Shift'
		print Shift

		print 'Shift : '+ str(Shift[0]) + ' mhjd days with correlation ' + str(Shift[1])
	

		#Plotting
		plot_data_points(time, lc_1, lc_2, avg[x], avg[x+1])
		for c1 in LC_1:
			plot_segment(c1, avg[x])
		for c2 in LC_2:
			plot_segment(c2, avg[x+1])
	show_plots()

