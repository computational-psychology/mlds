#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runs MLDS using the mlds wrapper in python (which calls R in the background)
and predicts the 2AFC thresholds at given standard values.
Uses bootstrap samples from MLDS to also provide thresholds with bootstrap CIs

@author: G. Aguilar, April 2020
"""

import mlds
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context('talk')
sns.set_style("ticks", {'xtick.direction': 'in', 'ytick.direction': 'in',
                        "xtick.major.size": 2, "ytick.major.size": 2})

name = 'data.csv'

# parameters for prediction form a scale
dp = [-2, -1, -0.5, 0.5, 1, 2]  # d' around st to be predicted from scale
k = 3  # spline degree, 3 for cubic spline
factor = 2.0  # factor correction
citype = 'BCa'  # bias-corrected and accelerated CIs
# citype='percentile'  # raw percentile CI

nsamples = 1000  # number of bootstrap samples 
res = 0.000015

# standards
sts = np.array([0.0303])  # , 0.0466, 0.0680])

# %% Estimating scale calling MLDS
ml = mlds.MLDSObject(name, boot=True, standardscale=False, verbose=False)
ml.nsamples = nsamples
ml.load()

# plotting
plt.figure(figsize=(5, 5))
plt.plot(ml.stim, ml.scale, c='k')
yerr = abs(ml.ci95 - ml.mns)
plt.errorbar(ml.stim, ml.mns, yerr=yerr, fmt="none", c='k')
sns.despine()
plt.xlabel('Stimulus')
plt.ylabel('Perceptual scale')
plt.show()

# runs GoF diagnostics
# ml.rundiagnostics()
# ml.plotdiags()

# %% Deriving threshold prediction from scale
print('estimating thresholds from scales')

# predict thresholds from scale at different d' around st
# mldsthrs is a pandas DataFrame, also saved as a CSV in the same folder
mldsthrs, xs, ys, scaleGLM = mlds.predict_thresholds(ml, sts, dp, k=k, factor=factor,
                                                     citype=citype, tol=0.01, res=res,
                                                     warn=True)

plt.figure(figsize=(5, 5))
plt.plot(ml.stim, ml.scale * factor, c='k')
yerr = abs(ml.ci95 - ml.mns)
plt.errorbar(ml.stim, ml.mns * factor, yerr=yerr, fmt="none", c='k')
sns.despine()
plt.xlabel('Stimulus')
plt.ylabel('Perceptual scale')
plt.plot(xs, ys, label='Interpolated scale')
plt.title('Interpolated scale')
plt.show()

# EOF
