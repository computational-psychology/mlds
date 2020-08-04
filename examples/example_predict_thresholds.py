#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Runs MLDS using the mlds wrapper in python (which calls R in the background)
and predicts the 2AFC thresholds at given standard values.
Uses bootstrap samples from MLDS to also provide thresholds with bootstrap CIs

It saves threhsold data in a JSON. 

Check out the following references for more details.

[1] Aguilar, Wichmann & Maertens (2017) Journal of Vision.
[2] Devinck & Knoblauch (2012) Journal of Vision.

@author: G. Aguilar, April 2020
"""

import mlds
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context('talk')
sns.set_style("ticks", {'xtick.direction': 'in', 'ytick.direction': 'in',
                        "xtick.major.size": 2, "ytick.major.size": 2})

# raw csv file of MLDS data
name = 'data.csv'

#### parameters for prediction from a perceptual scale

# interpolation of the perceptual scale
k = 3  # spline degree, 3 for cubic spline

# factor multiplication to  transform perceptual scale in units of d'
# assuming an internal representation that is a gaussian random variable with
# variance sigma, the decision variable in the MLDS model will also be 
# be a gaussian random variable with variance 4*sigma. MLDS estimates the overall 
# variance \hat{sigma} of the decision variable, and the inverse is the 
# the scale's height. Thus, to transform the height to units of the *internal*
# representation, we multiply by 2. (more detailed derivation in Appendix of [1])
factor = 2.0  

#
nsamples = 1000  # number of bootstrap samples 
res = 0.000015  # interpolation resolution. decrease if needed.

# Standards: at which stimulus levels you want to predict thresholds?
sts = np.array([0.25, 0.5, 0.75, 0.99])
# d' : thresholds will be derived to the following d' performance 
dp = [-2, -1, -0.5, 0.5, 1, 2]  # d' around the  standard at which we will read out thresholds
# (negative values indicate comparisons below the standard, positive values for 
# comparisons above the standard)


# type of confidence interval. 
citype = 'BCa'  # bias-corrected and accelerated CIs
# citype = 'percentile'  # raw percentile CI

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
                                                     warn=True, debug=False)

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

print(mldsthrs)

# EOF
