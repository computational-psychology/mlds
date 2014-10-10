# -*- coding: utf-8 -*-
"""
Functions to simulate and fit difference scales from a known power function.

@author: GA, May 2014
"""

import matplotlib.pyplot as plt
import sys, os
sys.path.append("/home/guille/git/slantfromtex/mlds")
import mlds # <- maxim. likelihood difference scaling tools
import csv
import numpy as np
import random
from math import sqrt
import uuid


def psi(s, e):
    """Simple power function"""
    return s**e
    

def simulateobserver(func, sigma, stim, nblocks=1, gamma=1):
    """ Simulate an observer responding according to a function
    
    """
    
    ## Simulation triads experiment
    triads, indxt, order = mlds.generate_triads(stim)   

    ## where to save the results
    filename = str(uuid.uuid4()) + '.csv'
    rfl = open( filename, 'wb')
    writer = csv.writer(rfl, delimiter=' ')
    writer.writerow(['Trial', 'Response', 's1', 's2', 's3', 'invord', 'i1', 'i2', 'i3'])
    
    # going through every trial
    for b in range(nblocks):
        for t in range( len(triads) ):
            
            # slants of current triad
            if order[t]==1:
                s1 = triads[t][2]
                s2 = triads[t][1]
                s3 = triads[t][0]
            else:
                s1 = triads[t][0]
                s2 = triads[t][1]
                s3 = triads[t][2]
            
            # get mean response 
            v1 = func(s1, gamma)
            v2 = func(s2, gamma)
            v3 = func(s3, gamma)
            
    
            # responses from cue
            # random.gauss accepts negative sigmas, it probably gets 
            # the absolute number. np.random.gauss, instead, doesnt, 
            # gives an error.. This could be problematic for unconstrained 
            # optimization routines that take sigma as parameter
  
            f1 =  random.gauss( v1, sigma) 
            f2 =  random.gauss( v2, sigma)
            # if we used the same value of f2, apparently introduces correlation and scale noise estimation (sigma) is lower     
            #f2b = f2
            # we should then sample twice for the value of stimulus 2, so variances are as predicted.
            f2b = random.gauss( v2, sigma)   
            f3 =  random.gauss( v3, sigma)    
             
            # decision variable
            delta = (f3 - f2b) - (f2-f1) 
            #delta = (v3 - v2)  - (v2 -v1) + random.gauss(0, sigma)
            
            if delta>0:
                response=1
            elif delta<0:
                response=0
            else:
                response=None
                print "None"
            
            if order[t]==1:
                response = int(not response)
            writer.writerow([t, response, "%.2f" % triads[t][0], "%.2f" % triads[t][1], "%.2f" % triads[t][2], order[t], indxt[t][0]+1, indxt[t][1]+1, indxt[t][2]+1])
             
    rfl.close()

    return filename



def runsimulation(sigmas, nblocks, nruns, plotflag=True, gamma=2, N=11):
    
    stim = np.linspace(0,1,N) # stimulus vector

    scales = np.zeros((len(sigmas), nruns, N))   # all simulation data in one variable (nsigmas x nruns x spacing)
    
    ## change of sigma (noise parameter), influence on standard and nonstandard scales
    if plotflag:
        plt.figure()
    
    for j,sigma in enumerate(sigmas):
           
        for i in range(nruns):
            f = simulateobserver(psi, sigma, stim, nblocks, gamma)      # simulate observer
            obs = mlds.MLDSObject( f , boot=False, standardscale=False) # initialize object
            obs.run()                                                   # runs Mlds in R
            scales[j][i] = obs.scale
            
            os.remove(f)  # removing simulation file
        
        if plotflag:
            plt.errorbar(stim, scales[j].mean(axis=0), yerr = scales[j].std(axis=0), label="$\sigma=%.3f$" % sigma, fmt='*')
          
        # from unconstrained scales, sigma by design is 1, leaving the maximum of the scale as 1/sigma. Then,
        # sigma can be extracted from an uncontrained scale by simply calculating 1/ the scale maximum.
        shat = 1/ scales[j][:,-1].mean()
        
        print "sigma2 estimated: %.4f" % shat**2
        print "sigma2 combined: %.4f" % (4*sigma**2)
        
        if plotflag:
            plt.plot(stim, psi(stim, gamma) * 1/sqrt(4*sigma**2),'k--')
        
    if plotflag:
        plt.xlabel("Stimulus")
        plt.ylabel("Difference scale")
        plt.xlim([0, 1.05])
        plt.legend(loc=2)
        plt.show()
        
    return scales





# displaying    
#s = np.linspace(0, 1, 100) # stimulus values
#gammas = np.linspace(1,5,11) # psi response 
#
#plt.figure()
#for i, g in enumerate(gammas):
#    plt.plot(s, psi(s, g), label="$\gamma$ = %.2f" % g)
#plt.legend(loc=2)
#plt.show()

    