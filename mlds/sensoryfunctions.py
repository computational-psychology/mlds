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

        if isinstance(x, list):
            v = [random.gauss(self.func(xx), self.sigmafunc(xx)) for xx in x]
        elif isinstance(x, float) or isinstance(x, int):
            v = random.gauss(self.func(x), self.sigmafunc(x))
        else:
            raise ValueError('provide either a float or int, or a list')

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
        
        if isinstance(x, list):
            v = [random.gauss(self.func(xx), self.sigmafunc(xx)) for xx in x]
        elif isinstance(x, float) or isinstance(x, int):
            v = random.gauss(self.func(x), self.sigmafunc(x))
        else:
            raise ValueError('provide either a float or int, or a list')

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
