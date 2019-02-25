Python wrapper for MLDS R package [![Build Status](https://travis-ci.org/TUBvision/mlds.svg?branch=master)](https://travis-ci.org/TUBvision/mlds)


Contents
========

It contains:

- a python implementation that wraps the R package MLDS. This wrapper makes easier to analyse the data obtained in MLDS experiments. It also provides the extra functionality to using multi-thread, making the bootstrap calculation much faster.

- utilities for designing MLDS experiments (method of triads and quadruples)
- functions to simulate an observer performing an MLDS experiment (so far only for the method of triads).


Requirements
============

- Python modules: numpy, subprocess, multiprocessing and rpy2 (>=2.3.10)

- R (>=3.0), with the *MLDS* and *psyphy* packages already installed.
- Optional but recommended: *snow* package in R.

To install python modules, use the package manager of your sytem (in Debian, apt-get)
R packages can be installed from CRAN or using the files provided in this repository (mlds/mlds/CRAN)



Installation
============

- Install all requirements.
- Clone the repository from github  (*git clone https://github.com/TUBvision/mlds.git*) 
- For adding multi-thread functionality, also install the package *snow* in R.
- Run the tests: go to the directory *mlds/mlds/test* and execute: *python -m unittest discover*. All tests should pass.



Usage examples
==============

- *mlds/example/example.py*  gives usage example for MLDS analysis
- *mlds/example/example_stim_generation.py*   gives usage example for designing the triads or quadruples.
- *mlds/example/example_simulation.py*   gives usage example for simulating an observer performing the method of triads.



Contact
=======
Questions? Feedback? Don't hesitate to fork, pull request, or 
contact me (guillermo@bccn-berlin.de)

This repository has so far only been tested in Linux (Debian Wheezy and Jessie, Ubuntu) 
