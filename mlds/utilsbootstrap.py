# -*- coding: utf-8 -*-
"""
Utilities for loading and saving python dictionaries

@author: G. Aguilar, Nov2014
update March 2020: to work on python3

"""

from scipy import stats
import numpy as np


def sumjacknife(retbt, e=1):

    N = len(retbt)
    T_i = np.zeros(N)

    for i in range(len(retbt)):
        T_i[i] = np.mean( np.delete(retbt, i) )

    T_dot = np.sum(T_i) / len(T_i)

    return np.sum((T_dot - T_i) ** e)
    
    
def getCI_percentile(array, conf=0.95):
    """
    
    Calculates confidence intervals from a boostrap sample array 
    based on 'percentiles'. 
    
    Input
    ------
    
    array : 1 dim np.array 
            Array containing the bootstrap samples of a statistic
    conf  : float
            Confidence desired. Default = 0.95    
    
    Output
    ------
    
    tuple: (low, high) C.I. values 
            Lower and higher boundary of confidence interval
    
    
    Reference
    ----------
    Efron & Tibshirani (1993). Introduction to the bootstrap. 
    Chapter 13-3
    
    """
    l = np.percentile(array, ((1-conf)/2.0)*100)
    u = np.percentile(array, (1-(1-conf)/2.0)*100)
    
    return l, u
	


def getCI_BCa(array, t, conf=0.95):
    """
    
    Calculates 'bias-corrected and accelerated' confidence intervals 
    from a boostrap sample array 
    
    Input
    ------
    
    array : 1 dim np.array 
            Array containing the bootstrap samples of a statistic
    t:  float
        estimated statistic from sample
        
    conf  : float
            Confidence desired. Default = 0.95    
    
    Output
    ------
    
    tuple: (low, high) C.I. values 
            Lower and higher boundary of confidence interval
    
    
    Reference
    ----------
    Efron & Tibshirani (1993). Introduction to the bootstrap. 
    Chapter 14-3. pp 185 - 186
    
    """
    a1 = (1-conf)/2.0
    a2 = 1-(1-conf)/2.0
    
    # bias part
    bias = stats.norm.ppf( np.sum(array<t)/float(array.size))
    
    if np.isinf(bias):
        return np.nan, np.nan
    
    E_l3 = sumjacknife(array, 3)
    E_l2 = sumjacknife(array, 2)
    
    # acceleration part
    acc =  E_l3 / (6* E_l2**(3.0/2.0))
    
    # corrected percentiles
    alpha1 = stats.norm.cdf( bias + ( stats.norm.ppf(a1) + bias ) / (1-acc*(stats.norm.ppf(a1) + bias )) )
    alpha2 = stats.norm.cdf( bias + ( stats.norm.ppf(a2) + bias ) / (1-acc*(stats.norm.ppf(a2) + bias )) )
    
    #print 'a1, a2 = (%.4f, %.4f)' % (alpha1, alpha2)
    
    l = np.percentile(array, alpha1*100)
    u = np.percentile(array, alpha2*100)
    
    
    return l, u
    
    
