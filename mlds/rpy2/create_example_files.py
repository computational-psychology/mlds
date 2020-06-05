# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 16:09:39 2014

@author: guille
"""

import sim_scales_powerfn
import os
import numpy as np


gammas = [1, 2, 1.5, 1.8]
sigma = 0.2

N = 11
nblocks = 15


stim = np.linspace(0,1,N)

for i in range(len(gammas)):
    f = sim_scales_powerfn.simulateobserver(sim_scales_powerfn.psi, sigma, stim, nblocks, gammas[i]) 
    os.rename(f, '%d.csv' % i ) 


