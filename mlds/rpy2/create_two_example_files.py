# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 16:09:39 2014

@author: guille
"""

import sim_scales_powerfn 
import os
import numpy as np


gamma1 = 2
gamma2 = 2.4
N = 11
nblocks = 15

sigma1 = 0.2
sigma2 = 0.2



stim = np.linspace(0,1,N)

f = sim_scales_powerfn.simulateobserver(sim_scales_powerfn.psi, sigma1, stim, nblocks, gamma1) 
os.rename(f, 'first.csv')

f = sim_scales_powerfn.simulateobserver(sim_scales_powerfn.psi, sigma2, stim, nblocks, gamma2) 
os.rename(f, 'second.csv')

