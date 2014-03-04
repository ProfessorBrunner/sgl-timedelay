#! /usr/bin/env python
# Robert J. Brunner
# February 2, 2014
#
# This file is directly executable: 
# rb $ ./generateLC.py
#
# Generate and plot quasar lightcurves following B. Kelly 2009 paper.
# Need to review all code. 
#
# Kelly argues based on fitting models to real lightcurves that optical
# variations are due to thermal fluctuations (i.e., longest timescales).

import numpy as np
import matplotlib.pyplot as plt

# Function to calculate stochastic integral
# Stochastic Integral doesn't require integrands in code?

def stochasticIntegral(tau):
    steps = 5 # Default number of sampling points for integral evaluation
    
    # Normal parameters
    mean = 0.0
    sigma = 1.0
    
    # Our integral sample points
    sample = np.arange(steps)
    
    # Our integration is simply summing up the random dt and the function value

    values = np.exp(-sample/tau) * np.random.normal(mean, sigma, steps)

    return (np.sum(values))

# Mass is in units of a solar mass
# Radius is in terms of the Schwarzschild Radius

def initTimescales(mass, radius, alpha):
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

# tsteps    Number of time steps
def generateLC(mass, radius, alpha, b, sigma, tsteps):

    # Define parameters:
    ntimes = 3 # Number of timescales
    
    # First get timescales
    tau = initTimescales(mass, radius, alpha)
    
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
            ssi[j] = sigma * stochasticIntegral(tau[j])
            
        fluxes[:,i] = exp_tau * fluxes[:, i - 1] + bt * (1.0 - exp_tau) + ssi
 
    return fluxes

# Main function, here we initialize values to match Kelly paper.
tsteps = 512
mass = 1.0E8
radius = 100
sigma = 1.0
b = 0.
alpha = 0.01

# Now generate timesteps and flux values
time = np.linspace(0,2555,512) #time is in days

flux = generateLC(mass, radius, alpha, b, sigma, tsteps)

#generate light curve plots
f, ax = plt.subplots(3, sharex=True)

ax[0].plot(time,flux[0,:])
ax[0].set_title('Light Crossing Light Curve')

ax[1].plot(time,flux[1,:])
ax[1].set_title('Orbital Light Curve')

ax[2].plot(time,flux[2,:])
ax[2].set_title('Thermal Light Curve')

plt.xlabel('Time [days]')
ax[1].set_ylabel('Flux [arbitrary units]')
plt.savefig('LightCurves.pdf')
plt.close()


