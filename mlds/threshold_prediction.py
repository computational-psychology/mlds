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
from scipy.interpolate import UnivariateSpline

from .utilsbootstrap import getCI_BCa, getCI_percentile


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def getvalue(xs, sp, st, d, tol):
    """
    Given a *monotonically increasing* spline (sp) and its x-vector (xs),
    finds the value corresponding to a d-prime difference (d) around the standard (st).
    Returns nan if value cannot be found inside a certain tolerance (tol).
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
    iv, val = find_nearest(ys, target)
    if abs(val - target) > tol:
        return np.nan
    else:
        return xs[iv]


def calculateCI(ret, retbt, Nsim, citype='percentile'):
    """
    TODO
    """
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
        retm, retl, retu = np.nan, np.nan, np.nan

    return retm, retl, retu


def getdprimefromspline(xs, sp, st, d, sp_bt=None, citype="percentile", tol=0.1, warn=True):
    """
    Obtain stimulus value d-prime units (d) around a standard (st), from the
    mean fitted data on a spline (sp), and c.i. around this mean from the
    bootstrapped data (sp_bt).
    C.I. are either 'percentile', or 'BCa' (bias corrected and accelerated)

    tol: estimation error tolerance in the scale dimension (y axis)

    """

    # get values for mean scale
    ret = getvalue(xs, sp, st, d, tol)

    if sp_bt is not None:
        # we get a spline with CI
        # variability estimation: bootstrap
        # get values for all bootstrap runs
        print("thresholds from bootstrap samples for st: %f, d'=%.1f" % (st, d))
        retbt = [getvalue(xs, s, st, d, tol) for s in sp_bt]

        # checks if there are at least 100 different values in the histogram
        if len(np.unique(retbt)) < 100:
            print("WARNING: "
                  "there are only %d different histogram values for the bootstrap thresholds... "
                  "increase resolution of x dimension" % len(np.unique(retbt)))

            if warn:
                import matplotlib.pyplot as plt
                plt.hist(retbt, 100)
                plt.show()

        # CI
        retm, retl, retu = calculateCI(ret, retbt, len(sp_bt), citype)

    else:
        retl, retu, retm = np.nan, np.nan, np.nan
        retbt = []

    sys.stdout.flush()  # flushing print output

    return ret, retm, retl, retu, retbt


def getalldprime(st, xs, sp, dp, sp_bt=None, citype="percentile", tol=0.1, warn=True):
    """
    Calculates getdprimefromspline() for all d'-values given as dp, for a given standard st.
    """
    data = {}
    ret_bt = {}

    for d in dp:
        ret, retm, retl, retu, retbt = getdprimefromspline(xs, sp, st, d, sp_bt, citype, tol, warn)
        data[d] = (ret, retm, retl, retu)
        ret_bt[d] = retbt

    return data, ret_bt


####################
def predict_thresholds(obsGLM, sts, dp, k=3, factor=2, citype="percentile",
                       tol=0.1, res=0.01, rangex=None, warn=True, save=True, debug=False):
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
    debug: bool

    
    
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

    df_columns = ['st', 'dprime', 'point_estimate', 'mean', 'CIl', 'CIh']

    if obsGLM.mns is not None:
        scaleGLM = obsGLM.mns * factor
        # scaleGLM_ci = obsGLM.ci95*2
        bootdata = True
    else:
        scaleGLM = obsGLM.scale * factor
        bootdata = False

    # k degree parameter, <=5
    sp = UnivariateSpline(obsGLM.stim, scaleGLM, k=k)
    # sp = interpolate.interp1d(obsGLM.stim, scaleGLM, kind=5)
    # sp = np.poly1d( np.polyfit(obsGLM.stim, scaleGLM,5) )

    if rangex is None:
        xs = np.arange(obsGLM.stim[0], obsGLM.stim[-1], res)  # resolution of spline fit
    else:
        xs = np.arange(rangex[0], rangex[1], res)  # resolution of spline fit

    # evaluate spline with xs as a vector
    ys = sp(xs)

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
                ret_bt = json.load(fp)
            results = []
            for st in sts:
                for d in dp:
                    ret = getvalue(xs, sp, st, d, tol)
                    retbt = ret_bt[str(st)][str(d)]
                    nboot = obsGLM.scalesbt.shape[1]

                    retm, retl, retu = calculateCI(ret, retbt, nboot, citype)
                    results.append([st, d, ret, retm, retl, retu])
            df = pd.DataFrame(results, columns=df_columns)

            if save:
                df.to_csv(fname_th, index=False)

        else:
            # proceed to calculations
            nboot = obsGLM.scalesbt.shape[1]
            scaleGLM_bt = obsGLM.scalesbt * factor

            sp_bt = [UnivariateSpline(obsGLM.stim, scaleGLM_bt[:, i], k=k) for i in range(nboot)]

            ## slow performance here, commented cos it's just for visualization
            # ys_bt = np.zeros((len(xs), nboot))
            # for i in range(nboot):
            #     ys_bt[:, i] = sp_bt[i](xs)
            # ys_l = np.percentile(ys_bt, 2.5, axis=1)
            # ys_h = np.percentile(ys_bt, 97.5, axis=1)
            # ys = np.vstack((ys, ys_l,  ys_h))
            # plt.plot(obsGLM.stim, scaleGLM, 'ko')
            # plt.plot(xs, ys[0,:])
            # plt.plot(xs, ys[1,:],'r')
            # plt.plot(xs, ys[2,:],'g')
            # plt.show()

            # get prediction for all standards and d.primes
            results = []
            ret_bt = {}

            params = (xs, sp, dp, sp_bt, citype, tol, warn)
            if debug:
                results_per_st = [getalldprime(st, *params) for st in sts]
            else:
                tasks_iterator = (delayed(getalldprime)(st, *params) for st in sts)
                parallelizer = Parallel(n_jobs=(multiprocessing.cpu_count()))
                results_per_st = parallelizer(tasks_iterator)

            for i, st in enumerate(sts):
                data, retbt = results_per_st[i]
                for dp, (ret, retm, retl, retu) in data.items():
                    results.append([st, dp, ret, retm, retl, retu])
                ret_bt[st] = retbt

            df = pd.DataFrame(results, columns=df_columns)

            if save:
                df.to_csv(fname_th, index=False)
                with open(fname_bt, 'w') as fp:
                    json.dump(ret_bt, fp)

    else:
        results = []
        for st in sts:
            for d in dp:
                ret, retm, retl, retu, _ = getdprimefromspline(xs, sp, st, d, tol=tol, warn=warn)
                results.append([st, d, ret, retm, retl, retu])
        df = pd.DataFrame(results, columns=df_columns)

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
