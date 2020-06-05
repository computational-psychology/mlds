# -*- coding: utf-8 -*-
"""
Example use of MLDS class

data.csv was created by simulation.

"""

import matplotlib.pyplot as plt

# appends directory where mlds.py is, and imports
import sys
sys.path.append("../")
import mlds


## creating MLDS object and running analysis
obs = mlds.MLDSObject('data.csv', standardscale = False,
                      boot=True, verbose=True)

print("Running analysis...")
# obs.run() # always runs the analysis, overwrites the .MLDS file if exists.
obs.load()  # loads results from .MLDS file; runs the analysis if doesn't exist


## plotting
print("plotting..")
plt.figure()
plt.errorbar(obs.stim, obs.mns, yerr=abs(obs.ci95 - obs.mns), color='#4C72B0',
             linewidth=2, elinewidth=2 )

plt.xlabel('Stimulus')
plt.ylabel('Difference Scale')
plt.xlim((-0.1, 1.1))
plt.show()


## if you want, you can create an .R file containing all the R commands
# that are run in the background
obs.saveRcommands()


## Goodness of fit
print("GoF measures:")
print('AIC: %f, DAF: %f' % (obs.AIC, obs.DAF))
# if you want to run deviance bootstrap analysis, uncomment the lines below
# print("..deviance analysis ...")
# obs.rundiagnostics()
# print('p-val: %f' % obs.prob)
# # plotting deviance
# fig = obs.plotdiags()
# fig.show()

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
