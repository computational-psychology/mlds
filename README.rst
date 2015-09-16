Python wrapper for MLDS R package.


Contents
========

It contains a python implementation that wraps the
R package MLDS. This wrapper makes easier to analyse the data
obtained in MLDS experiments. It also provides the extra
functionality to using multi-thread for the bootstrap.


It also contains:

- utilities for designing MLDS experiments (method of triads and quadruples)
- functions to simulate an observer performing an MLDS experiment.


Requirements
============

- Python modules: rpy2, subprocess, multiprocessing
- R, with *MLDS* package already installed. 
- Optional but recommended: *snow* package


Installation
============

- Install all requirements. 
- Clone the repository from github.
- For multi-thread: ensure that you can log in via ssh to your own localhost. For that, check *instructions_snowpackage.txt*
- Run the tests: in directory *mlds/mlds/test*, execute *python -m unittest discover*. All tests should pass.



Examples
========

- *mlds/example/example.py*  gives usage example for MLDS analysis
- *mlds/example/example_stim_generation.py*   gives usage example for designing the triads or quadruples.



Contact
=======
Questions, feedback? Don't hesitate to fork, pull request, or 
contact me (guillermo@bccn-berlin.de)
