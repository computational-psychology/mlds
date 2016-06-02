# -*- coding: utf-8 -*-
"""
Class definition of sensory represention functions to be used in observer
simulations.

@author: G. Aguilar, June 2015
"""

import random
import numpy as np
from scipy import interpolate


class PowerSensoryFunc:
    """
    Generic power sensory function with gaussian noise. Noise can be
    independent or to depend on the sensory function value in a linear fashion.
    Vectorial function, it can be called for a single value or a list of values

    """

    def __init__(self):

        # power function
        self.exponent = 3  # default value
        self.a = 1
        self.c = 0

        self.func = lambda x: self.a*x**self.exponent + self.c

        # definition of relationship of noise with sensory rep value
        # if sigmamax == sigmamin, the noise is constant
        # otherwise a linear relationship can be defined
        # i.e. linearly increasing or decreasing noise with sensory value
        self.sigmamax = 0.3  # default values
        self.sigmamin = 0.1

        # noise function
        self.sigmafunc = lambda x: (self.sigmamax - self.sigmamin) * self.func(x) + self.sigmamin

    def __call__(self, x):
        """  returns a sample at stimulus intensity x  """

        if isinstance(x, float) or isinstance(x, int):
            v = random.gauss(self.func(x), self.sigmafunc(x))
        else:
            try:
                v = [random.gauss(self.func(xx), self.sigmafunc(xx)) for xx in x]
            except:
                raise ValueError('provide either a float, an int, a list or a numpy vector')

        return v


class QuadraticSensoryFunc:
    """
    Generic quatratic sensory function with gaussian noise. Noise can be
    independent or to depend on the sensory function value in a linear fashion.
    Vectorial function, it can be called for a single value or a list of values

    """

    def __init__(self):

        # power function
        self.a = 2  # default value
        self.b = 2
        self.c = 0

        self.func = lambda x: self.a*x**2 + self.b*x + self.c

        # definition of relationship of noise with sensory rep value
        # if sigmamax == sigmamin, the noise is constant
        # otherwise a linear relationship can be defined
        # i.e. linearly increasing or decreasing noise with sensory value
        self.sigmamax = 0.3  # default values
        self.sigmamin = 0.1

        # noise function
        self.sigmafunc = lambda x: (self.sigmamax - self.sigmamin) * self.func(x) + self.sigmamin

    def __call__(self, x):
        """  returns a sample at stimulus intensity x  """
        
        if isinstance(x, float) or isinstance(x, int):
            v = random.gauss(self.func(x), self.sigmafunc(x))
        else:
            try:
                v = [random.gauss(self.func(xx), self.sigmafunc(xx)) for xx in x]
            except:
                raise ValueError('provide either a float, an int, a list or a numpy vector')

        return v



class Cue2DSensoryFunc:
    """
    Generic sensory function that follows a given 2D cue with gaussian noise.
    Noise can be independent or to depend on the sensory function value in a linear fashion.
    Not vectorial function.. to be called for a single value

    Input
    ----------
    cuefilename : string
                  filename (complete path) of .npy file containing stimulus
                  and cue values

    """

    def __init__(self, cuefilename, normalize=True):

        # load cue data
        try:
            self.rawcue = np.load(cuefilename)
        except:
            print "Error loading cue file"
            raise

        # normalize if asked
        if normalize:
            cmax = self.rawcue[1, :].max()
            cmin = self.rawcue[1, :].min()

            a = 1.0 / (cmax - cmin)
            b = -cmin / (cmax - cmin)
            self.rawcue[1, :] = a * self.rawcue[1, :] + b

        # defines cue function, as a linear interpolation object
        self.func = interpolate.interp1d(self.rawcue[0, :], self.rawcue[1, :], kind='linear')

        # define default noise parameters
        self.sigmamax = 0.3  # default values
        self.sigmamin = 0.1

        # noise function
        self.sigmafunc = lambda x: (self.sigmamax - self.sigmamin) * self.func(x) + self.sigmamin


    def __call__(self, x):
        """  returns a sample at stimulus intensity x  """

        v = random.gauss(self.func(x), self.sigmafunc(x))
        return v




class PowerSensoryFunc_corr:

    def __init__(self):

        # power function
        self.exponent = 2  # default values
        self.a = 1.0
        self.c = 0.0

        # gives the mean of the function
        self.func = lambda x: self.a*x**self.exponent + self.c

        self.sigmamin = 0.0
        self.sigmamax = 0.0

        # noise function
        self.sigmafunc = lambda x: (self.sigmamax - self.sigmamin) * self.func(x) + self.sigmamin
        
        # covariance matrix
        self.cov = np.zeros((4, 4))
        
        # nondiagonal elements
        self.sigma_nd = 0.00
        

    def setcov(self, x, corranchor=False):
        """ Set covariance matrix for a given stimulus vector """
               
        # non diagonal values are the correlation terms,, fixed for all combinations
        self.cov = np.ones((4, 4)) * (self.sigma_nd**2)  
                

        # diagonal indices depend on the stimulus value itself
        for i in range(4):
            self.cov[i,i] = self.sigmafunc(x[i])**2
            

        # covariance values between 1 and 2 are full correlated
        if corranchor:
            self.cov[1, 2] = 1.0
            self.cov[2, 1] = 1.0
            
        

    def __call__(self, x):
        """  returns a sample at stimulus intensities in the vector x,
        with covariance cov"""
        
        ##
        if len(x)<3:
            raise(ValueError, "stimulus vector should be at least have 3 values")
        elif len(x)==3:
            X = [x[0], x[1], x[1], x[2]]  # three values are converted to four values vector
            corranchor = True
        elif len(x)==4:
            X = x
            corranchor=False
        else:
            raise(ValueError, "stimulus vector should be 3 or 4 values in length")
               
            
        #print X
        # get means 
        mu_triad = [self.func(i) for i in X]
        #print mu_triad
        
        # and covariance
        self.setcov(X, corranchor)

        # get the draws from gauss distribution multivariate
        v = np.random.multivariate_normal(mu_triad, self.cov)
        return v
        

# how to calculate correlation matrix from covariance matrix     
# D = np.sqrt(np.diag(fn.cov))
# (1.0/D)*fn.cov*(1.0/D)




if __name__ == "__main__":
    import matplotlib.pyplot as plt

    cuefilename= "cue.npy"
    obj = Cue2DSensoryFunc(cuefilename)

    stim = np.array(range(0, 70, 5))
    vals = np.array([obj.func(x) for x in stim])

    plt.plot(stim, vals, 'bo')
    plt.plot(obj.rawcue[0, :], obj.rawcue[1, :], 'b')

    obj.sigmamax = 0.2
    obj.sigmamin = 0.005

    randomvals = np.array([obj(x) for x in stim])
    plt.plot(stim, randomvals, 'ro')
    plt.show()
