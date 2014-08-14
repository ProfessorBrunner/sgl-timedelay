import scipy.optimize as opt
import numpy as np
import pylab as plt
from scipy import misc
from scipy import misc
import math

'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Fitting for two dimentional gaussians and one dimentional 
sersics.
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
'''
def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()


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

	#Statistics for getting the initial estimate
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
	Mean =  Val_sum / (float)( m * n )


	#initial_guess = (Max,Max_x,Max_y,1,1,0,0)
	initial_guess = (100,4,Max_x,100,4,Max_y,0)


	#popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), z.ravel(), p0=initial_guess)
	#data_fitted = twoD_Gaussian((x, y), *popt)

	popt, pcov = opt.curve_fit(sersic, (x, y), z.ravel(), p0=initial_guess)
	data_fitted = sersic((x, y), *popt)

	plt.figure()
	#plt.imshow(z, cmap=plt.cm.jet, extent=(x.min(), x.max(), y.min(), y.max()))
	plt.imshow(z)

	plt.figure()
	plt.imshow(data_fitted.reshape(m, n))

	plt.figure()
	plt.imshow(z-data_fitted.reshape(m, n))

	plt.show()

