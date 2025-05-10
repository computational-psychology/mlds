# Python wrapper for MLDS R package 

[![Tests](https://github.com/computational-psychology/mlds/actions/workflows/ci-tests.yml/badge.svg)](https://github.com/computational-psychology/mlds/actions/workflows/ci-tests.yml)
[![DOI](https://zenodo.org/badge/42587765.svg)](https://zenodo.org/doi/10.5281/zenodo.12658147)

This python package contains:

- a python implementation that wraps the R package MLDS. This wrapper makes easier to analyse the data obtained in MLDS experiments. It also provides the extra functionality to using multi-thread, making the bootstrap calculation much faster.

- utilities for designing MLDS experiments (method of triads and quadruples)

- functions to simulate an observer performing an MLDS experiment (so far only for the method of triads).


## Requirements

- Python >= 3.8 with numpy, subprocess, multiprocessing, and rpy2

- R, with the *MLDS*, *psyphy* and *snow* packages

Python module dependencies are installed automatically.

R packages must be installed manually, either from CRAN (see below)
or using the files provided in this repository (mlds/CRAN).


## Installation

- Install R
 
- Install the required R packages. You can do it from inside R with

```R
install.packages(c("MLDS", "psyphy", "snow"))
```

or from the command line with

```bash
R -e 'install.packages(c("MLDS", "psyphy", "snow"))'
```

- Install the mlds wrapper (this package). In the console run

`pip install https://github.com/computational-psychology/mlds/tarball/master`. 

The python dependencies will be installed automatically.


## Quick start

Given a CSV file containing the data from a triads experiment, called 'data_triads.csv'
the code

```python
import matplotlib.pyplot as plt
import mlds

obs = mlds.MLDSObject('data_triads.csv', standardscale=False, 
                      boot=True, verbose=False)
```

creates an object `obs` that will talk to R, pass the data, do the fit and
return the results to python. Until now the object is just initalized; 
to run the actual analysis we do

```python
obs.run()  # the wrapper now sends the commands to R, and R runs the fitting. 
```

Now we can get the stimulus values and estimated scales with

```python
print(obs.stim)
print(obs.scale)
```

and plot the perceptual scale with

```python
plt.figure()
plt.errorbar(obs.stim, obs.scale)
plt.xlabel('Stimulus')
plt.ylabel('Difference Scale')
plt.show()
```

A more detailed usage example can be found in *examples/example.py*.



## Usage examples

- *examples/example.py* gives usage example for data from a triads experiment.
- *examples/example_stim_generation.py* gives usage example for designing the triads or quadruples given the stimulus intensities.
- *examples/example_simulation.py* shows how to simulate an observer performing the method of triads.
- *examples/example_predict_thresholds.py* shows how to predict discrimination thresholds from a perceptual scale, using signal detection theory assumptions. See Aguilar, Wichmann & Maertens (2017) for details.



Contact
=======
Questions? Feedback? Don't hesitate to ask Guillermo Aguilar (guillermo.aguilar@mail.tu-berlin.de)

This implementation has a testing suite that gets run on every release. Testing is done in python >=3.8 and Ubuntu.

