# -*- coding: utf-8 -*-
"""
Example of simulated observer performing MLDS

"""
import matplotlib.pyplot as plt
import numpy as np
import mlds

#######
# 1. First, we set a stimulus dimension of interest
N = 11
stim = np.linspace(0, 1, N)

#######
# 2. Then, we need to set a sensory function that maps the sensory dimension
# with the perceived dimension. 
# This package includes some rudimentary functions
from mlds.sensoryfunctions import PowerSensoryFunc

# sets a quadratic sensory function with exponent 2:  x^2 + epsilon, where x
# is the stimulus dimension. epsilon is the gaussian-distributed noise
# added to the sensory representation
fn = PowerSensoryFunc()
fn.exponent = 3.0

# it sets a fixed noise value 
fn.sigmamax = 0.05
fn.sigmamin = 0.05


######
# 3. Now we simulate the observer performing the method of triads.
# The simulated results are stored in the file fname
# (Reduce nblocks to 1 if you just want to see the result format.)
fname = mlds.simulateobserver(fn, stim, nblocks=15)


#######
# 4. Finally, we want to estimate the scale from the simulated data.
obs = mlds.MLDSObject(fname, boot=True, standardscale=False, verbose=False)

# number of samples for the calculation of bootstrap confidence intervals.
# here only 100 just for demostration purposes and speed. Recommeded is 
# at least 1000, we've used 10000 just to be sure (but that adds much more 
# computation time)
obs.nsamples = 100 

# calling to run the analysis
obs.run()  # this takes a while as bootstrap is done
obs.printinfo()


# we can plot the scale
# the shape of the scale should coincide with a power function with exponent 2
mlds.plotscale(obs)
plt.xlim([-0.1, 1.1])
plt.show()

# and the noise estimated by MLDS should be the double of the noise
# introduced at the sensory level in the simulation
print("noise estimate from MLDS: %.3f" % (1 / obs.mns[-1]))
print("which must be the double of the sensory noise %.3f" % fn.sigmamax)


# We can also check the goofness of fit of the model. As we have simulated
# an observer that exactly behaves like MLDS assumes a observer decides, then
# the goodness of fit tests should all be OK. 
obs.rundiagnostics()

print("GoF measures:")
#print('AIC: %f, DAF: %f' % (obs.AIC, obs.DAF))
print('p-val: %f' % obs.prob) 
# p-value should be not signifficant, that is, we do not discard the MLDS model
# as a model that can explain the data.


# EOF
