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
                       boot=True)

obs.load()  # load results from file if it exists; otherwise runs the analysis
# obs.run() # always runs the analysis, overwrites a results file if exists.


## plotting
plt.errorbar(obs.stim, obs.mns, yerr=abs(obs.ci95 - obs.mns), color='#4C72B0', 
                linewidth=2, elinewidth=2 )

plt.xlabel('stimulus')
plt.ylabel('difference scale')
plt.xlim((-0.1, 1.1))
plt.show()


## Goodness of fit
print 'AIC: %f, DAF: %f' % (obs.AIC, obs.DAF)

obs.rundiagnostics()
print 'p-val: %f' % obs.prob


## access to bootstrap samples
obs.scalesbt

## access to deviance residuals
obs.residuals
