import numpy as np
import scipy as sp
import scipy.ndimage
from scipy import misc
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
import sys
from dcor import * 
from polar import *

'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Experimental algorithm for finding the center of a foreground
lensing galaxy.
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
'''

if __name__ == '__main__':

	if (len(sys.argv)!=2):
		print 'Please specify one input file'
		exit(0)
	
	z = misc.imread(sys.argv[1],True)
	m , n = z.shape

	Corr = []
	cVal = 0

	Max = 0
	Max_x = 0
	Max_y = 0 
	Val_sum = 0
	for i in range (364,365):
		for j in range (413,n-1):
			print 'Coordinates: ' + str( i ) + ', ' +str (j)
			Val_sum += z[i][j]												#Find Sum to calculate Mean
			try:			
				polar_image = plot_polar_image(z, (i,j))
				row_mean, row_median, row_index = get_statistics(polar_image)
			
				Model_Fit = fit_profile(row_mean,row_index)
				cVal = dcor(np.asarray(Model_Fit),np.asarray(row_mean)).find_correlation()
				Corr.append([cVal,i,j])				
				print cVal
				if(cVal>Max):
					Max = cVal
					Max_x = i
					Max_y = j 
			except(RuntimeError):
				print 'Optimized parameters not found'				
				pass
				

	Mean =  Val_sum / (float)( m * n )
	print Corr
	print 'MAX'
	print m, n	
	print (Max_x, Max_y, cVal)	
	polar_image = plot_polar_image(z, (Max_x,Max_y))
	row_mean, row_median, row_index = get_statistics(polar_image)
	Model_Fit = fit_profile(row_mean,row_index)



	'''
	Plotting
	'''
	Residual = []
	for i in range (0,len(row_mean)):
		Residual.append(abs(Model_Fit[i] - row_mean[i]))

	plt.figure()
	plt.plot(row_index,row_mean,label='Mean Row Intensity')
	plt.plot(row_index,row_median,label='Median Row Intensity')
	plt.plot(row_index,Model_Fit,label='Sersic Profile')
	plt.plot(row_index,Residual,label='Residual Profile')
	plt.legend()
	plt.xlabel('Row Number')
	plt.ylabel('Intensity')
	plt.show()
