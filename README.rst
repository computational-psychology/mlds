Python wrapper for MLDS R package.


Contents
========

It contains a python implementation that wraps the
R package MLDS. This wrapper makes easier to analyse the data
obtained in MLDS experiments. It also provides the extra
functionality to using multi-thread, making the bootstrap calculation
much faster.


It also contains:

- utilities for designing MLDS experiments (method of triads and quadruples)
- functions to simulate an observer performing an MLDS experiment.


Requirements
============

- Python modules: numpy, rpy2, subprocess, multiprocessing
- R, with *MLDS* package already installed. 
- Optional but recommended: *snow* package in R.


Installation
============

- Install all requirements. 
- Clone the repository from github.
- For multi-thread: ensure that you can log in via ssh to your own localhost. For that, check *instructions_snowpackage.txt*
- Run the tests: in directory *mlds/mlds/test*, execute *python -m unittest discover*. All tests should pass.



Usage examples
==============

- *mlds/example/example.py*  gives usage example for MLDS analysis
- *mlds/example/example_stim_generation.py*   gives usage example for designing the triads or quadruples.
- *mlds/example/example_simulation.py*   gives usage example for simulating an observer performing the method of triads.



Contact
=======
Questions? Feedback? Don't hesitate to fork, pull request, or 
contact me (guillermo@bccn-berlin.de)

This repository has so far only been tested in Linux. 
