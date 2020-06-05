# -*- coding: utf-8 -*-
"""
@author: guille
"""
import sys
sys.path.append('../')
import mlds
from mlds import plotscale
import mlds_gam
import matplotlib.pyplot as plt

def setplotproperties(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['left'].set_bounds(0, 1)
    ax.set_xlim(( -0.2, 1.2))
    ax.set_ylim(( -0.2, 1.2))
    ax.spines['bottom'].set_bounds(0, 1)
    ax.tick_params(labelsize=12)
    ax.set_xlabel('Stimulus')
    ax.set_ylabel('Difference scale', labelpad=5)
    return ax


st = True
files=['0.csv', '1.csv', '2.csv']
gammas = [1, 2, 1.5]

glms = []
for f in files:
    O = mlds.MLDSObject(f, boot= True, standardscale=st, save=True)
    #O.parallel=True
    #O.run()
    O.load()
    glms.append( O )

obsGAM = mlds_gam.MLDSGAMCompare(files, standardscale=st)
obsGAM.run()

pval = obsGAM.anovares[4]
print(pval)

fig = plt.figure(figsize=(8, 6))
ax = plt.gca()

l = plt.plot(obsGAM.stim, obsGAM.scales, linewidth=2)
for o, obs in enumerate(glms):
    plotscale(obs, '$\gamma = %.2f$' % gammas[o], l[o].get_color(), linewidth=0, marker='o')
    # plt.plot(obs.stim, obs.scale, linewidth=0, marker='o', markerfacecolor = l[o].get_color(), label='$\gamma = %.2f$' % gammas[o])

ax = setplotproperties(ax)
plt.legend(loc=2, frameon=False)
plt.title('Example using power fn',  {'fontsize': 10} )
fig.savefig('example.pdf')
plt.show()
        
    


             
