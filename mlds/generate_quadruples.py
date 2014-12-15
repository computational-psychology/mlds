# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 14:45:45 2013

@author: guille
"""

import numpy as np
import itertools
import random

stim = np.append(np.linspace(0, 0.9, 10), 0.98)
#stim= list(np.append(np.linspace(0, 0.9, 10), 0.98))

numstim = len(stim)

# getting all possible combinations in a list, C(n,k)
allTrials = [list(s) for s in itertools.combinations(range(numstim), 4)]

# randomize order of trials
random.shuffle(allTrials)

# randomize order or pairs 
topbot = np.random.binomial(1, 0.5, len(allTrials) )

for t in range(len(allTrials)):
    if topbot[t]:
        allTrials[t] = [ allTrials[t][2], allTrials[t][3], allTrials[t][0], allTrials[t][1] ]
        
# create list of stimulus quadruples
quadruples = [stim[t]   for t in allTrials]        
        
        
        
        
        


