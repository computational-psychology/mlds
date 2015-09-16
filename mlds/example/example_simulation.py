# -*- coding: utf-8 -*-
"""
Example of simulated observer performing MLDS

"""

import matplotlib.pyplot as plt
import numpy as np
# appends directory where mlds.py is, and imports
import sys
sys.path.append("../")
import mlds


#######
# 1. First, we set a stimulus dimension of interest
N = 11
stim = np.linspace(0, 1, N)

#######
# 2. Then, we need to set a sensory function that maps the sensory dimension
# with the perceived dimension. 
# This repository includes some rudimentary functions
from sensoryfunctions import PowerSensoryFunc

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
#    The simulated results are stored in the file fname
fname = mlds.simulateobserver(fn, stim, nblocks=15)


#######
# 4. Finally, we want to estimate the scale from the simulated data.
obs = mlds.MLDSObject(fname, boot=True, standardscale=False)
obs.load()  # this takes a while as bootstrap is done
obs.printinfo

# we can plot the scale
# the shape of the scale should coincide with a power function with 
# exponent 2
mlds.plotscale(obs)
plt.xlim([-0.1, 1.1])
plt.show()

# and the noise estimated by MLDS should be the double of the noise
# introduced at the sensory level in the simulation
print "noise estimate from MLDS: %.3f" % (1/obs.mns[-1])
print "which must be the double of the sensory noise %.3f" % fn.sigmamax


# finally, GoF measures should be OK, as we are simulating an observer
# that actually performs the decision model assumed by MLDS
# obs.rundiagnostics()
# print "GoF measures:"
# print 'AIC: %f, DAF: %f' % (obs.AIC, obs.DAF)
# print 'p-val: %f' % obs.prob

# EOF
