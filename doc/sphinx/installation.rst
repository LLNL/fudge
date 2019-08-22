.. _installation:

Installation
============

Fudge is designed to be simple to install and use. For basic use the package only requires 
Python (versions 2.4-2.7, not version 3.0 and higher).  Python is likely already present 
on all unix-based systems (if in doubt, try entering "python" on the command-line),
and is available for all platforms at http://www.python.org/getit/.

To install the package:
-----------------------

1) Obtain the current version (available online at https://ndclx4.bnl.gov/gf/project/gnd/ or ftp://gdo-nuclear.ucllnl.org/), 
   and place it on your system.

2) To unzip the package, type

    ::

        $tar -xvf fudge-4.0.0.tar.gz
    
3) The unzipped directory is named 'fudge-4.0.0'. Navigate into the new directory, and type 

    ::

        $make
    
   this builds extra tools for better performance. If you encounter errors at this step, you can still use the ENDF
   translation tool but advanced features will be unavailable.

4) You are now set to use Fudge! For advanced use, see the optional components below.


[Optional] Extra python packages
--------------------------------
  
For more advanced use, Fudge depends on these additional, optional packages:

* gnuplot and Gnuplot.py (for plotting, http://gnuplot-py.sourceforge.net)
* matplotlib (alternate style of plotting, http://matplotlib.sourceforge.net)
* numpy (for some advanced features including checking and manipulating covariance matrices, http://numpy.scipy.org)


[Optional] Setting Environment Variables
----------------------------------------

For general use of the fudge package, some changes should be made to your
computer's environment. The following lines make the required change. 
Note that <path_to_fudge> indicates the path to
the directory containing the Fudge README.txt file:

* On unix with bash, ksh, etc, put this line in the .bashrc or equivalent:

    ::

        $export PYTHONPATH=$PYTHONPATH:<path_to_fudge>

* On csh, tcsh or similar, do:
    
    ::
    
        $setenv PYTHONPATH $PYTHONPATH:<path_to_fudge>

* On Windows, the environment variable should be added to the registry 
  (see for example http://www.support.tabs3.com/main/R10463.htm)


[Optional] Building Extensions
------------------------------
  
Fudge is primarily written in python, which is easily portable across many different
systems. Some advanced features of fudge (computationally intensive tasks) are however 
implemented in C/C++ or fortran, and must be compiled before use. To build these extensions, 
use the Makefile included in the main Fudge directory:

::

    $ make

If you encounter trouble with installing, building extensions or setting environment 
variables, please let us know (email: mattoon1@llnl.gov).

