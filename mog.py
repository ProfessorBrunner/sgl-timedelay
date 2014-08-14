# -*- coding: utf-8 -*-
import scipy.optimize as opt
import numpy as np
import pylab as plt
from scipy import misc
from scipy import misc
import sys
from polar import *

'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Mixture of gaussians mimicing different types of astronomical
lensing profiles. Scipy.optimize is used for fitting.
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
'''

# Generalized gaussian
def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g.ravel()


# Theta = 0 (Reduce a  parameter)
def normpdf_2D((x,y), amplitude, mu_x, mu_y, sigma_x, sigma_y, offset):
    xo = float(mu_x)
    yo = float(mu_y)    
    u = ( ((x-mu_x) * (x-mu_x)) / (abs(sigma_x) * abs(sigma_x)) ) + ( ((y-mu_y) * (y-mu_y)) / (abs(sigma_y) * abs(sigma_y)) )
    y = offset + amplitude * np.exp(-u/2) / (np.sqrt(2*np.pi) * abs(sigma_x) * abs(sigma_y) )
    return y.ravel()

'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Exponential Approximation by mixture of Gaussians 
	Qexp (ξ) = exp(−α^exp [|ξ| − 1])
Approximation by a mixture of 4, 6 and 8 gaussians with relative weights and
Variance as given in "Replacing standard galaxy profiles with mixtures of 
Gaussians - David W. Hogg & Dustin Lang
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
'''
def exp_mix_2D_4((x,y), amplitude, offset, x0, y0, k):
    w = [0.09733,1.12804,4.99846,5.63632]
    v = [0.12068,0.32730,0.68542,1.28089]
    y = (w[0] * normpdf_2D((x,y), amplitude, x0, y0, k * v[0], k * v[0], offset) + w[1] * normpdf_2D((x,y), amplitude, x0, y0, k * v[1], k * v[1], offset) + w[2] * normpdf_2D((x,y), amplitude, x0, y0, k * v[2], k * v[2], offset) + w[3] * normpdf_2D((x,y), amplitude, x0, y0, k * v[3], k * v[3], offset) )/sum(w)
    return y.ravel()

def exp_mix_2D_6((x,y), amplitude, offset, x0, y0, k):
    w = [0.00735,0.09481,0.63572,2.60077,5.42848,3.16445]
    v = [0.05072,0.13756,0.28781,0.53195,0.91209,1.50157]
    y = (w[0] * normpdf_2D((x,y), amplitude, x0, y0, k * v[0], k * v[0], offset) + w[1] * normpdf_2D((x,y), amplitude, x0, y0, k * v[1], k * v[1], offset) + w[2] * normpdf_2D((x,y), amplitude, x0, y0, k * v[2], k * v[2], offset) + w[3] * normpdf_2D((x,y), amplitude, x0, y0, k * v[3], k * v[3], offset) + w[4] * normpdf_2D((x,y), amplitude, x0, y0, k * v[4], k * v[4], offset) + w[5] * normpdf_2D((x,y), amplitude, x0, y0, k * v[5], k * v[5], offset) )/sum(w)
    return y.ravel()

def exp_mix_2D_8((x,y), amplitude, offset, x0, y0, k):
    w = [0.00077,0.01017,0.07313,0.37184,1.39736,3.56100,4.74338,1.78684]
    v = [0.02394,0.06492,0.13581,0.25095,0.42942,0.69675,1.08885,1.67302]
    y = (w[0] * normpdf_2D((x,y), amplitude, x0, y0, k * v[0], k * v[0], offset) + w[1] * normpdf_2D((x,y), amplitude, x0, y0, k * v[1], k * v[1], offset) + w[2] * normpdf_2D((x,y), amplitude, x0, y0, k * v[2], k * v[2], offset) + w[3] * normpdf_2D((x,y), amplitude, x0, y0, k * v[3], k * v[3], offset) + w[4] * normpdf_2D((x,y), amplitude, x0, y0, k * v[4], k * v[4], offset) + w[5] * normpdf_2D((x,y), amplitude, x0, y0, k * v[5], k * v[5], offset) + w[6] * normpdf_2D((x,y), amplitude, x0, y0, k * v[6], k * v[6], offset) + w[7] * normpdf_2D((x,y), amplitude, x0, y0, k * v[7], k * v[7], offset) )/sum(w)
    return y.ravel()


'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Lux profile approximation by mixture of Gaussians 
	Qexp (ξ) = exp(−α^lux [|ξ| − 1]) for |ξ| < 3
	Qexp (ξ) = exp(−α^lux [|ξ| − 1]) [1 − [|ξ| − 3]^2 ]^2 for 3 < |ξ| < 4
	Qexp (ξ) = 0 for 4 < |ξ|
Approximation by a mixture of 4, 6 and 8 gaussians with relative weights and
Variance as given in "Replacing standard galaxy profiles with mixtures of 
Gaussians - David W. Hogg & Dustin Lang
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
'''
def lux_mix_2D_4((x,y), amplitude, offset, x0, y0, k):
    w = [0.07275,0.86763,4.33214,6.48325]
    v = [0.10938,0.29694,0.62601,1.19571]
    y = (w[0] * normpdf_2D((x,y), amplitude, x0, y0, k * v[0], k * v[0], offset) + w[1] * normpdf_2D((x,y), amplitude, x0, y0, k * v[1], k * v[1], offset) + w[2] * normpdf_2D((x,y), amplitude, x0, y0, k * v[2], k * v[2], offset) + w[3] * normpdf_2D((x,y), amplitude, x0, y0, k * v[3], k * v[3], offset) )/sum(w)
    return y.ravel()

def lux_mix_2D_6((x,y), amplitude, offset, x0, y0, k):
    w = [0.00235,0.03080,0.22336,1.17949,4.33874,5.99821]
    v = [0.03465,0.09405,0.19785,0.37413,0.67894,1.22540]
    y = (w[0] * normpdf_2D((x,y), amplitude, x0, y0, k * v[0], k * v[0], offset) + w[1] * normpdf_2D((x,y), amplitude, x0, y0, k * v[1], k * v[1], offset) + w[2] * normpdf_2D((x,y), amplitude, x0, y0, k * v[2], k * v[2], offset) + w[3] * normpdf_2D((x,y), amplitude, x0, y0, k * v[3], k * v[3], offset) + w[4] * normpdf_2D((x,y), amplitude, x0, y0, k * v[4], k * v[4], offset) + w[5] * normpdf_2D((x,y), amplitude, x0, y0, k * v[5], k * v[5], offset) )/sum(w)
    return y.ravel()

def lux_mix_2D_8((x,y), amplitude, offset, x0, y0, k):
    w = [0.00007,0.00098,0.00736,0.04404,0.24005,1.18175,4.31918,5.97985]
    v = [0.01092,0.02966,0.06241,0.11794,0.21345,0.38155,0.68169,1.22635]
    y = (w[0] * normpdf_2D((x,y), amplitude, x0, y0, k * v[0], k * v[0], offset) + w[1] * normpdf_2D((x,y), amplitude, x0, y0, k * v[1], k * v[1], offset) + w[2] * normpdf_2D((x,y), amplitude, x0, y0, k * v[2], k * v[2], offset) + w[3] * normpdf_2D((x,y), amplitude, x0, y0, k * v[3], k * v[3], offset) + w[4] * normpdf_2D((x,y), amplitude, x0, y0, k * v[4], k * v[4], offset) + w[5] * normpdf_2D((x,y), amplitude, x0, y0, k * v[5], k * v[5], offset) + w[6] * normpdf_2D((x,y), amplitude, x0, y0, k * v[6], k * v[6], offset) + w[7] * normpdf_2D((x,y), amplitude, x0, y0, k * v[7], k * v[7], offset) )/sum(w)
    return y.ravel()


'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Luv profile approximation by mixture of Gaussians 
Approximation by a mixture of 8 and 10 gaussians with relative weights and
Variance as given in "Replacing standard galaxy profiles with mixtures of 
Gaussians - David W. Hogg & Dustin Lang
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
'''
def luv_mix_2D_8((x,y), amplitude, offset, x0, y0, k):
    w = [0.04263,0.24013,0.68591,1.51937,2.83627,4.46467,5.72441,5.60990]
    v = [0.01496,0.03166,0.06471,0.13017,0.26170,0.53592,1.15464,2.89864]
    y = (w[0] * normpdf_2D((x,y), amplitude, x0, y0, k * v[0], k * v[0], offset) + w[1] * normpdf_2D((x,y), amplitude, x0, y0, k * v[1], k * v[1], offset) + w[2] * normpdf_2D((x,y), amplitude, x0, y0, k * v[2], k * v[2], offset) + w[3] * normpdf_2D((x,y), amplitude, x0, y0, k * v[3], k * v[3], offset) + w[4] * normpdf_2D((x,y), amplitude, x0, y0, k * v[4], k * v[4], offset) + w[5] * normpdf_2D((x,y), amplitude, x0, y0, k * v[5], k * v[5], offset) + w[6] * normpdf_2D((x,y), amplitude, x0, y0, k * v[6], k * v[6], offset) + w[7] * normpdf_2D((x,y), amplitude, x0, y0, k * v[7], k * v[7], offset) )/sum(w)
    return y.ravel()

def luv_mix_2D_10((x,y), amplitude, offset, x0, y0, k):
    w = [0.01468,0.09627,0.28454,0.63005,1.19909,2.03195,3.07255,4.10682,4.83948,4.94943]
    v = [0.01190,0.02210,0.03995,0.07117,0.12586,0.22240,0.39593,0.71922,1.37549,3.13117]  
    y = (w[0] * normpdf_2D((x,y), amplitude, x0, y0, k * v[0], k * v[0], offset) + w[1] * normpdf_2D((x,y), amplitude, x0, y0, k * v[1], k * v[1], offset) + w[2] * normpdf_2D((x,y), amplitude, x0, y0, k * v[2], k * v[2], offset) + w[3] * normpdf_2D((x,y), amplitude, x0, y0, k * v[3], k * v[3], offset) + w[4] * normpdf_2D((x,y), amplitude, x0, y0, k * v[4], k * v[4], offset) + w[5] * normpdf_2D((x,y), amplitude, x0, y0, k * v[5], k * v[5], offset) + w[6] * normpdf_2D((x,y), amplitude, x0, y0, k * v[6], k * v[6], offset) + w[7] * normpdf_2D((x,y), amplitude, x0, y0, k * v[7], k * v[7], offset) + w[8] * normpdf_2D((x,y), amplitude, x0, y0, k * v[8], k * v[8], offset) + w[9] * normpdf_2D((x,y), amplitude, x0, y0, k * v[9], k * v[9], offset) )/sum(w)
    return y.ravel()

def exp_mix_2D((x,y), amplitude, mu_x, mu_y, sigma_x, sigma_y, offset, w0 , w1, w2 ):
    y =  w0 * normpdf_2D((x,y), amplitude, mu_x, mu_y, sigma_x, sigma_y, offset) + w1 * normpdf_2D((x,y), amplitude, mu_x, mu_y, sigma_x, sigma_y, offset) + w2 * normpdf_2D((x,y), amplitude, mu_x, mu_y, sigma_x, sigma_y, offset) 
    return y.ravel()


def sersic((x,y),kx,nx,x0,ky,ny,y0,s):
	S = s * math.e ** (-1.0 * kx * (  (x - x0) ** (1.0/nx) ) ) * math.e ** (-1.0 * ky * (  (y - y0) ** (1.0/ny) ) ) 
	return S.ravel()


if '__main__':

	if (len(sys.argv)!=2):
		print 'Please specify one input file'
		exit(0)

	z = misc.imread(sys.argv[1],True)
	m , n = z.shape
	x, y = np.mgrid[:m, :n]

	# Get stats for an initial guess
	Max = 0
	Max_x = 0
	Max_y = 0 
	Val_sum = 0
	for i in range (0,m-1):
		for j in range (0,n-1):
			Val_sum += z[i][j]												#Find Sum to calculate Mean
			if(z[i][j]>Max):
				Max = z[i][j]
				Max_x = i
				Max_y = j 
	#Mean =  Val_sum / (float)( m * n )

	initial_guess = (Max,0,Max_x,Max_y,2)
	popt, pcov = opt.curve_fit(luv_mix_2D_10, (x, y), z.ravel(), p0=initial_guess)
	data_fitted = luv_mix_2D_10((x, y), *popt)
	
	print 'Optimum fitting parameters [amplitude, offset, x0, y0, k]'
	print popt

	f, axarr = plt.subplots(3, sharex=True)
	axarr[0].imshow(z, cmap='Greys')
	axarr[1].imshow(data_fitted.reshape(m, n), cmap='Greys')
	axarr[2].imshow(z-data_fitted.reshape(m, n), cmap='Greys')
	plt.figure()
	plt.imshow(z-data_fitted.reshape(m, n), cmap='Greys')
	polar_image = plot_polar_image(z-data_fitted.reshape(m, n), (popt[2],popt[3]))
	row_mean, row_median, row_index = get_statistics(polar_image)
	plt.figure()
	plt.plot(row_index,row_mean,label='Mean Row Intensity')
	plt.plot(row_index,row_median,label='Median Row Intensity')
	plt.show()
