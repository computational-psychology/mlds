# -*- coding: utf-8 -*-
"""
Example use of MLDS class and the method of triads

data.csv was created by simulation.

@author: G. Aguilar.
"""

import matplotlib.pyplot as plt
import mlds

plt.style.use('seaborn-v0_8-notebook')

## creating MLDS object and running analysis
obs = mlds.MLDSObject('data_triads.csv',
                      standardscale=False, # standard scale normalizes range from 0 to 1
                      boot=True) # to calculate confidence intervals via bootstrap
obs.nsamples = 1000

print("Running analysis...")
# the wrapper now sends the commands to R, and R runs the fitting.
obs.run()  # always runs the analysis, overwrites the .MLDS file if exists.

# The wrapper makes R save everything (scale values, bootstrap CIs etc)
# in a .MLDS file. This file can be opened manually with R itself, or
# used by this wrapper to plot and get useful information about
# the fitted model (see below).

# Instead of run(), one can also call .load()
#obs.load()   # loads results from .MLDS file and does not run the analysis again.

# This option can be useful as bootstrapping takes long. If there is no change in our
# data, we can call run() only once, and later just call load() to read the data
# from the previously saved MLDS file. If the MLDS file does not exist,
#  load() will automatically call run().


## plotting the scale
print("plotting...")
fig, ax = plt.subplots(1, 1, figsize=(5, 4))
ax.errorbar(obs.stim, obs.scale, yerr=abs(obs.ci95 - obs.mns), color='#4C72B0',
             linewidth=2, elinewidth=2)

ax.set_xlabel('Stimulus')
ax.set_ylabel('Perceptual scale')
ax.spines[['right', 'top']].set_visible(False)
plt.tight_layout()
fig.savefig('example_triads.png')
plt.show()


## if you want, you can create an .R file containing all the R commands
# that are run in the background. This might be useful to debug problems.
# obs.saveRcommands()


## Goodness of fit metrics
print("GoF measures:")
print('AIC: %f, DAF: %f' % (obs.AIC, obs.DAF))

## Deviance analysis - takes long as it runs a simulation
print("..deviance analysis ...")
obs.rundiagnostics()
print('p-val: %f' % obs.prob)

# plotting deviance. Plot as in Knoblauch & Maloney's book.
fig = obs.plotdiags()
plt.show()


#######################################################################
## important values that are accesible directly
# stimuli:                      obs.stim
# scale:                        obs.scale
# bootstrapped scale mean:      obs.mns
# bootstrapped scale CI95:      obs.ci95

# bootstrap samples:            obs.scalesbt
# deviance residuals:           obs.residuals

# GOF measures:                 obs.AIC, obs.DAF, obs.prob


## EOF
