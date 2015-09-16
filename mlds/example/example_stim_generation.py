# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 15:51:47 2013

@author: guille
"""


import numpy as np
# appends directory where mlds.py is, and imports
import sys
sys.path.append("../")
import mlds


## definition of stimulus values
stim = np.linspace(0.0, 1.0, 11)


## generating quadruples
quads, indx, topbot = mlds.generate_quadruples(stim)


## generating triads
triads, indxt, topbot = mlds.generate_triads(stim)


# EOF
