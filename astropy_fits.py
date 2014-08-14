from scipy import misc
import numpy as np
from astropy.modeling import models, fitting
import astropy.modeling.functional_models as M
import matplotlib.pyplot as plt
import scipy as sp
import scipy.ndimage

'''
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Fitting models using astropy framework (bad results).
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
'''


def index_coords(data, origin=None):
    """Creates x & y coords for the indicies in a numpy array "data".
    "origin" defaults to the center of the image. Specify origin=(0,0)
    to set the origin to the lower left corner of the image."""
    ny, nx = data.shape
    if origin is None:
        origin_x, origin_y = nx // 2, ny // 2
    else:
        origin_x, origin_y = origin
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    x -= origin_x
    y -= origin_y
    return x, y

def cart2polar(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta

def polar2cart(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

def reproject_image_into_polar(data, origin=None):
    """Reprojects a 3D numpy array ("data") into a polar coordinate system.
    "origin" is a tuple of (x0, y0) and defaults to the center of the image."""
    ny, nx = data.shape
    if origin is None:
        origin = (nx//2, ny//2)

    # Determine that the min and max r and theta coords will be...
    x, y = index_coords(data, origin=origin)
    r, theta = cart2polar(x, y)

    r_i = np.linspace(r.min(), r.max(), nx)
    # Make a regular (in polar space) grid based on the min and max r & theta
    theta_i = np.linspace(theta.min(), theta.max(), ny)
    theta_grid, r_grid = np.meshgrid(theta_i, r_i)

    # Project the r and theta grid back into pixel coordinates
    xi, yi = polar2cart(r_grid, theta_grid)
    xi += origin[0] # We need to shift the origin back to 
    yi += origin[1] # back to the lower-left corner...
    xi, yi = xi.flatten(), yi.flatten()
    coords = np.vstack((xi, yi)) # (map_coordinates requires a 2xn array)
    zi = sp.ndimage.map_coordinates(data, coords, order=1).reshape((nx, ny))

  
    return zi, r_i, theta_i

    # Reproject each band individually and the restack
    # (uses less memory than reprojection the 3-dimensional array in one step)
    #bands = []
    #for band in data.T:
    #    zi = sp.ndimage.map_coordinates(band, coords, order=1)
    #    bands.append(zi.reshape((nx, ny)))
    #output = np.dstack(bands)
    #return output, r_i, theta_i



z = misc.imread('arc.JPG',True)
m , n = z.shape

x, y = np.mgrid[:m, :n]

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

print Mean


# Fit the data using astropy.modeling
p_init = M.Gaussian2D(Max,Max_x,Max_y,1,10)


#p_init = M.Ring2D(210,Max_x,Max_y,0,30)

f = fitting.LinearLSQFitter()
p = f(p_init, x, y, z)

# Plot the data with the best-fit model
plt.figure(figsize=(8,2.5))
plt.subplot(1,4,1)
plt.imshow(z, interpolation='nearest',vmin = 0 , vmax = Max)
plt.title("Data")
plt.subplot(1,4,2)
plt.imshow(p(x, y), interpolation='nearest',vmin = 0 , vmax = Max)
plt.title("Model")
plt.subplot(1,4,3)
plt.imshow(z - p(x, y), interpolation='nearest',vmin = 0 , vmax = Max)
plt.title("Residual")
plt.subplot(1,4,4)
polar_grid, r, theta = reproject_image_into_polar(z - p(x, y), (Max_x,Max_y) )
plt.imshow( polar_grid, vmin = 0 , vmax = Max)
plt.title("Polar")
#plt.colorbar(orientation='horizontal')
plt.colorbar()
plt.show()
