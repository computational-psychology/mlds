# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 15:51:47 2013

@author: guille
"""


import numpy as np
import mlds


## definition of stimulus values
stim = np.linspace(0.0, 1.0, 11)


## generating quadruples
quads, indxq, invordq = mlds.generate_quadruples(stim)


## generating triads
triads, indxt, invordt = mlds.generate_triads(stim)


# EOF
