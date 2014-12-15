# -*- coding: utf-8 -*-
"""
Utilities for MLDS  experiments

@author: G. Aguilar, Nov 2013
"""

import numpy as np
import itertools
import random


def generate_quadruples(stim):
    """
    Generate quadruples for difference scaling experiment.
    
    Usage:   quads, indx = generate_quadruples( stim )
    
    
    Input:  stim:   numpy array of stimulus values
    
    Output: quads:  list of all possible non-overlapping quadruples (stimulus values)
            indx:    (idem)   but indices values
    
    Ported from routines runQuadExperiment.R   from R package "MLDS"
    
    Reference: 
        Knoblauch, K. and Maloney, L. T. (2008) MLDS: Maximum likelihood difference scaling in R.
        Journal of Statistical Software, 25:2, 1–26, http://www.jstatsoft.org/v25/i02.

    """
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

    # returns the quadruples, and the indices
    return (quadruples, allTrials, topbot)
    
    
    
def generate_triads(stim):
    """
    Generate triads for difference scaling experiment.
    
    Usage:   traids, indx = generate_quadruples( stim )
    
    
    Input:  stim:   numpy array of stimulus values
    
    Output: triads:  list of all possible non-overlapping triads (stimulus values)
            indx:    (idem)   but indices values


    Ported from routines runTriadExperiment.R   from R package "MLDS"
    
    Reference: 
        Knoblauch, K. and Maloney, L. T. (2008) MLDS: Maximum likelihood difference scaling in R.
        Journal of Statistical Software, 25:2, 1–26, http://www.jstatsoft.org/v25/i02.
        
    """
    numstim = len(stim)

    # getting all possible combinations in a list, C(n,k)
    allTrials = [list(s) for s in itertools.combinations(range(numstim), 3)]
    
    # randomize order of trials
    random.shuffle(allTrials)
    
    # randomize order or pairs 
    topbot = np.random.binomial(1, 0.5, len(allTrials) )
    
    for t in range(len(allTrials)):
        if topbot[t]:
            allTrials[t] = [ allTrials[t][2], allTrials[t][1], allTrials[t][0] ]
            # careful here to make non-overlapping triads.. when inverting, the stim 1 must stay in the middle.
     
    # create list of stimulus triads
    triads = [stim[t]   for t in allTrials]

    # returns the triads, and the indices
    return (triads, allTrials, topbot)
 
    
    
    
    
    
