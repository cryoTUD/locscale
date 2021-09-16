Electron Microscopy Data Analytical Tool (EMDA)
===============================================

EMDA is an importable Python library for Electron Microscopy map and model manipulations written in Python3.
EMDA comes under MPL-2.0 license. 

A selected set of EMDA’s core functionalities include:

- Map input/output and format conversion
- Resolution related functionalities
- Map statistics calculation in Image and Fourier space
- Likelihood based map fitting and Bayesian map averaging
- Likelihood based magnification refinement
- Local correlation map calculation for map and model validation
- FSC based map-model validation


Installation
-------------  
To install EMDA, just type following line in the terminal

.. code-block:: 

   $ pip install emda

All necessary files will be installed under pythonx.x/site-packages/

e.g.: /Users/ranganaw/anaconda3/lib/python3.6/site-packages

To update EMDA, just type

.. code-block:: 

   $ pip install -U emda

To test the success of installation, please run the following command in the terminal

.. code-block:: 

   $ emda_test
   
All tests will pass, if the installation is successful.
   

Dependencies
-------------
EMDA has several dependencies in addition to general python modules (e.g. Numpy, Scipy).
They are Pandas, Gemmi, Mrcfile and Matplotlib.
All these dependencies will be automatically checked and installed during
EMDA installation.


How to use EMDA
---------------
There are two ways to access EMDA’s underlying functionalities. EMDA is available as a standalone command line tool to ensure its easy usage. Each functionality is callable with a keyword followed by a set of arguments. Also, EMDA has a simple to use, yet powerful, API for advanced usage documented under emda_methods.

.. toctree::
   :maxdepth: 4
   :caption: Contents:

   emda_methods
   emda_cmd_caller
