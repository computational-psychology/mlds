# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 15:56:35 2014

@author: guille
"""

import sys
sys.path.append('../')
import mlds

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
MLDS = importr("MLDS")


obs1 = mlds.MLDSObject('first.csv', standardscale =False)
obs2 = mlds.MLDSObject('second.csv', standardscale =False)

obs1.run()
obs2.run()


# loading R datafiles .MLDS that were previously saved
robjects.r['load']("%s" % obs1.Rdatafile)
obs1mlds = robjects.r['obs.mlds']

robjects.r['load']("%s" % obs2.Rdatafile)
obs2mlds = robjects.r['obs.mlds']


# 
