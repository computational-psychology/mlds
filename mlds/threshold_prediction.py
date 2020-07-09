# -*- coding: utf-8 -*-
"""
Functions to derive predictions from MLDS data in terms of SDT.
Use function scale_to_2afc  to feed an MLDSObject, the desired standards and
d-prime around each standard.

@author: G. Aguilar, Nov2014, Jan2015
update March 2020: to work on python3

"""
import json
import multiprocessing
import os.path
import sys

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import optimize
from scipy.interpolate import UnivariateSpline

from .utilsbootstrap import getCI_BCa, getCI_percentile


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def findroots_spline(sp, offset=0):
    """
    Find the x values on which sp(x) = offset.
    In other words, find the root(s) of a spline when moved vertically
    an amount offset.

    """

    f = lambda x: sp(x) - offset

    # find root via fsolve, starting at values 0 and 45 deg
    r = [optimize.fsolve(f, 0, xtol=0.0001),
         optimize.fsolve(f, 45, xtol=0.0001),
         optimize.fsolve(f, 90, xtol=0.0001)]
    r = np.array(r).round(2)
    return np.unique(r)


def findroots_spline_new(sp, offset=0):
    """
    Find the x values on which sp(x) = offset.
    In other words, find the root(s) of a spline when moved vertically
    an amount offset.

    Find the roots numericaly by fsolve, only at the spline knots 
    
    """

    f = lambda x: sp(x) - offset

    # find root via fsolve, starting at the spline knots values
    r = []
    for k in sp.get_knots():
        r.append(optimize.fsolve(f, k))
    r = np.array(r)

    return np.unique(r)


def findroots_spline_zerocrossings(xs, ys, offset=0):
    """
    Find the x values on which sp(x) = offset.
    In other words, find the root(s) of a spline when moved vertically
    an amount offset. 
    
    Find the roots numericaly by calculating zero crossings

    """
    zero_crossings = np.where(np.diff(np.sign(ys - offset)))[0]
    r = xs[zero_crossings]
    r = r.round(2)
    return np.unique(r)


#####
def getvalue(xs, sp, st, d, tol):
    """
    Given a spline (sp) and its x-vector (xs), finds the value corresponding
    to a d-prime difference (d) around the standard (st).
    Returns nan if value cannot be found inside a certain tolerance (tol)

    """

    # Workflow:
    # 0. get spline values y
    ys = sp(xs)

    # 1. set standard st, get id and value of discrete vector xs
    iv, val = find_nearest(xs, st)
    assert (abs(st - val) < tol)  # check that st value is found

    # 2. get target value at y = st + d
    target = ys[iv] + d

    # 3. find value of x corresponding to that target

    # profiling diff algorithms
    # X = findroots_spline(sp, target)  # **1
    X = findroots_spline_new(sp, target)  # **2 the fastest
    # X = findroots_spline_zerocrossings(xs, sp(xs), target) # **3

    X = X[X > 0]

    # 4. check that interpolated value in spline is actually close to desired value,
    # discard if findroots does find a value that does not correspond to original spline fit
    if len(X) == 0:
        ret = [np.nan]
    else:
        ret = []
        for i in range(len(X)):
            xi, xv = find_nearest(xs, X[i])
            err = abs(target - ys[xi])
            ret.append(np.nan if err > tol else xv)

    return ret


def getdprimefromspline(xs, sp, st, d, sp_bt=None, citype="percentile", tol=0.1, warn=True):
    """
    Obtain stimulus value d-prime units (d) around a standard (st), from the
    mean fitted data on a spline (sp), and c.i. around this mean from the
    bootstrapped data (sp_bt).
    C.I. are either 'percentile', or 'BCa' (bias corrected and accelerated)
    
    tol: estimation error tolerance in the scale dimension (y axis)

    """

    # get values for mean scale,
    ret = getvalue(xs, sp, st, d, tol)

    if len(ret) > 1:
        ret = np.max(ret)
        # if d > 0:
        #    # takes the
        #    ret = ret[np.argmax(np.array(ret))]
        # elif d < 0 :
        #    ret = ret[np.argmin(np.array(ret))]

    else:
        ret = ret[0]

    # if there's no value of threshold return, dont bother to do the same with the
    # bootstrap samples, immediately return 
    if np.isnan(ret):
        retl, retu, retm = np.nan, np.nan, np.nan
        retbt_r = np.array([np.nan])
        sys.stdout.flush()  # flushing print output
        return ret, retm, retl, retu, retbt_r
  

    if sp_bt is not None:  # we get a spline with CI
        # variability estimation: bootstrap
        # get values for all bootstrap runs
        print("thresholds from bootstrap samples for st: %f, d'=%.1f" % (st, d))
        retbt = np.zeros((len(sp_bt)))
        n = 0

        for i in range(len(sp_bt)):
            val = np.array(getvalue(xs, sp_bt[i], st, d, tol))

            val = val[np.logical_not(np.isnan(val))]

            if len(val) == 0:
                retbt[i] = np.nan
            elif len(val) == 1:
                retbt[i] = val[0]
            elif len(val) > 1:
                # takes the highest value
                retbt[i] = val[np.argmax(np.array(val))]
                n += 1
                # print i
                # print val

                # retbt[i] = val[np.argmin(np.array(val))]
                # if d > 0:
                #    retbt[i] = val[np.argmax(np.array(val))]
                # elif d < 0 :
                #    retbt[i] = val[np.argmin(np.array(val))]
                # else:
                #    print "something is very wrong with your choice of d'"

        # all values are nan, warning and returning nans
        if np.all(np.isnan(retbt)):
            print("WARNING: there are NaNs in all bootstrap samples, thresholds cannot be derived")
            retl, retu, retm = np.nan, np.nan, np.nan
            retbt_r =  np.array([np.nan])
            
            
        else:

            # print "%d times two values were recovered" % n
            retbt_r = np.copy(retbt)
    
            # checks if there are at least 100 different values in the histogram
            if len(np.unique(retbt_r)) < 100:
                print("WARNING: there are only %d different histogram values" % len(np.unique(retbt_r)))
                print("for the bootstrap thresholds.. increase resolution of x dimension")
    
                if warn:
                    import matplotlib.pyplot as plt
                    plt.hist(retbt, 100)
                    plt.show()

            ########################### CI
            retm, retl, retu = calculateCI(ret, retbt, len(sp_bt), citype)

    else:
        retl, retu, retm = np.nan, np.nan, np.nan
        retbt_r = np.array([np.nan])

    sys.stdout.flush()  # flushing print output

    return ret, retm, retl, retu, retbt_r


#########################
def calculateCI(ret, retbt, Nsim, citype='percentile'):
    if (len(retbt) > Nsim * 0.1) & (ret is not np.nan):

        retbt = np.array(retbt)[~np.isnan(retbt)]

        ## CI calculation *** See note below
        if citype == "percentile":
            retl, retu = getCI_percentile(retbt)

        elif citype == "BCa":
            retl, retu = getCI_BCa(retbt, ret)

        else:
            raise NameError("Type of confidence interval unknown")

        retm = np.percentile(retbt, 50.0)  # median

    else:
        retl, retu, retm = np.nan, np.nan, np.nan

    return retm, retl, retu


def getinparallel(st, xs, sp, dp, sp_bt, citype, tol, warn, bootdata):
    data = {}
    retbt = {}

    for d in dp:
        if bootdata:
            tup = getdprimefromspline(xs, sp, st, d, sp_bt, citype, tol=tol, warn=warn)
            data[d] = tup[0:4]
            retbt[d] = list(tup[4])
        else:
            data[d] = getdprimefromspline(xs, sp, st, d, tol=tol, warn=warn)

    return data, retbt


####################
def predict_thresholds(obsGLM, sts, dp, k=3, factor=2, citype="percentile",
                       tol=0.1, res=0.01, rangex=None, warn=True, save=True):
    """
    Takes a MLDS Class object, and calculates a prediction for performance on
    a 2AFC task at different standards.

    Arguments
    -----------

    obsGLM: MLDS Object

    sts: list, ndarray
         Contains the value of the standards that will be used as anchors

    dp: list
        Contains the value of d-prime to be calculated around every standard

    k: int
        Degree of polynomial for spline fitting

    factor: double
        Factor of re-parametrization of the scale. Default: 2

    citype: str {'percentile', 'BCa'}
        Confidence interval type

    res: float
         resolution in the x dimension. Default: 0.1
    
    tol: float
         error tolerance for the interpolation in the y dimension (d' dimension)
            Default: 0.01

    rangex: list
    warn: bool
    save: bool

    
    
    Return
    ----------

    data: dict
        Dictionary containin estimated thresholds for all sts and d' asked

    xs, ys: stimulus vector (xs) and spline fit for the scale (ys) in d'
    
    scaleGLM: scale in d'


    Reference
    ----------
    Devinck et al. (JoV 2014)
    Aguilar et al. (JoV 2017)

    """

    if obsGLM.mns is not None:
        scaleGLM = obsGLM.mns * factor
        # scaleGLM_ci= obsGLM.ci95*2
        bootdata = True
    else:
        scaleGLM = obsGLM.scale * factor
        bootdata = False

    # k degree parameter, <=5
    # s = 1
    sp = UnivariateSpline(obsGLM.stim, scaleGLM, k=k)
    # sp = interpolate.interp1d(obsGLM.stim, scaleGLM, kind=5)
    # sp= np.poly1d( np.polyfit(obsGLM.stim, scaleGLM,5) )

    if rangex is None:
        xs = np.arange(obsGLM.stim[0], obsGLM.stim[-1], res)  # resolution of spline fit
    else:
        xs = np.arange(rangex[0], rangex[1], res)  # resolution of spline fit

    # evaluate spline with xs as a vector
    ys = sp(xs)

    ##
    if bootdata:

        # if file exists, then it does not calculate once again
        fname_th = "%s_s2afc_%d_%s_%s.csv" % (obsGLM.Rdatafile.split(".")[0], k, factor, citype)
        fname_bt = "%s_retbt_%d_%s.json" % (obsGLM.Rdatafile.split(".")[0], k, factor)

        if os.path.isfile(fname_th):
            print("loading previously calculated c.i.")
            df = pd.read_csv(fname_th)

        elif os.path.isfile(fname_bt):
            print("loading previously bootstrapped thresholds")
            print("and calculating CI")
            with open(fname_bt, 'r') as fp:
                retbt = json.load(fp)
            data = []
            for st in sts:
                for d in dp:
                    ret = getvalue(xs, sp, st, d, tol)
                    if len(ret) > 1:
                        ret = np.max(ret)
                    else:
                        ret = ret[0]
                    tup = calculateCI(ret, retbt[str(st)][str(d)], obsGLM.scalesbt.shape[1], citype)
                    data.append([st, d, ret, tup[0], tup[1], tup[2]])
            df = pd.DataFrame(data,
                              columns=['st', 'dprime', 'point_estimate', 'mean', 'CIl', 'CIh'])

            if save:
                df.to_csv(fname_th, index=False)

        else:
            ### proceed to calculations
            nboot = obsGLM.scalesbt.shape[1]
            scaleGLM_bt = obsGLM.scalesbt * factor

            sp_bt = [UnivariateSpline(obsGLM.stim, scaleGLM_bt[:, i], k=k) for i in range(nboot)]
            ys_bt = np.zeros((len(xs), nboot))

            for i in range(nboot):
                ys_bt[:, i] = sp_bt[i](xs)

            ## slow performance here, commented cos it's just for visualization
            # ys_l = np.percentile(ys_bt, 2.5, axis=1)
            # ys_h = np.percentile(ys_bt, 97.5, axis=1)
            # ys = np.vstack((ys, ys_l,  ys_h))
            # plt.plot(obsGLM.stim, scaleGLM, 'ko')
            # plt.plot(xs, ys[0,:])
            # plt.plot(xs, ys[1,:],'r')
            # plt.plot(xs, ys[2,:],'g')
            # plt.show()

            #### get prediction for all standards and d.primes
            data = []
            retbt = {}

            ncpus = multiprocessing.cpu_count()
            parallelizer = Parallel(n_jobs=ncpus)

            tasks_iterator = (
            delayed(getinparallel)(st, xs, sp, dp, sp_bt, citype, tol, warn, bootdata) for st in
            sts)

            result = parallelizer(tasks_iterator)

            for i, st in enumerate(sts):
                retbt[st] = result[i][1]
                for dp, th in result[i][0].items():
                    data.append([st, dp, th[0], th[1], th[2], th[3]])

            df = pd.DataFrame(data,
                              columns=['st', 'dprime', 'point_estimate', 'mean', 'CIl', 'CIh'])

            if save:
                df.to_csv(fname_th, index=False)
                with open(fname_bt, 'w') as fp:
                    json.dump(retbt, fp)

    else:
        data = []
        for st in sts:
            for d in dp:
                data.append([st, d] + list(getdprimefromspline(xs, sp, st, d, tol=tol, warn=warn)))
        df = pd.DataFrame(data,
                          columns=['st', 'dprime', 'point_estimate', 'mean', 'CIl', 'CIh', 'bt_r'])

        print("variability not estimated: it needs bootstrapped data")

    return df, xs, ys, scaleGLM

##############################################################################
## Note about BCa CI calculation
# as in pypsignifit, CI are calculated as 'bias-corrected and accelerated' CI,
# as described in Efron (1993), pp 185 - 186.
# This is the same procedure used in psignifit, as seen in its sourcecode,
# specifically in file bootstrap.cc,  functions determineBCa(),
# which computes the bias and acceleration parameters,, and in file
# psignidata.py with computes the whole CI range considering the two terms.

# By passing citype="", the function here calculated percentile CI, which is
# not well suited cos it has been shown that CI in this way are way too small
# Wichmann (2001)b, see also Rasmussen (1987, 1988)
