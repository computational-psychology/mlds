# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 15:51:47 2013

@author: guille
"""


import numpy as np
import mlds

## definition of stimulus dimension, and values (spacing)
stim = np.append(np.linspace(0, 0.9, 10), 0.98)


## generating quadruples
quads, indx, topbot = mlds.generate_quadruples(stim)


## generating triads
triads, indxt, topbot = mlds.generate_triads(stim)


