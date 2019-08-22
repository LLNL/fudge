Welcome to FUDGE (For Updating Data and Generating Evaluations), a nuclear data platform built in python that permits
reading, viewing, modifying, writing and processing of nuclear data. FUDGE is a product of the Computational Nuclear
Physics Group (CNP) at Lawrence Livermore National Lab (LLNL).

The current release of Fudge focuses on introducing and extending the new Generalized Nuclear Data (GND) format,
intended to modernize storage of nuclear data. The package includes tools to translate other formats to and from GND,
plus tools for testing and visualizing GND-formatted data.

If you are encountering Fudge and GND for the first time, we strongly
recommend you look at the html documentation and tutorial available inside the
package at fudge-4.0.0/doc/html/index.html, where you will find much more
detail about installing and using the software.

Contents of this README:
    - 1) Getting Started
    - 2) Basic Use (translating ENDF files)
    - 3) Getting Help

1) Getting Started:
  Fudge is designed to be simple to install and use. For basic use the package only requires python (version 2.4-2.7).
Python is likely already present on all unix-based systems (if in doubt, try entering "python" on the command-line),
and is available for all platforms at www.python.org/getit

  To install the package:
    1) obtain the current version (available online at nuclear.llnl.gov), and place it on your system.
    2) to unzip the package, type
        >tar -xvf fudge-4.0.0.tar.gz	(or latest version)
    3) The unzipped directory is named 'fudge-4.0.0' (or latest version). Navigate into the new directory, and type 
        >make
       this builds extra tools for better performance. If you encounter errors at this step, you can still use the ENDF
       translation tool but advanced features will be unavailable.
    3) You are now set to use Fudge! For advanced use, see the optional components below.

  [Optional] Extra python packages: For more advanced use, Fudge depends on these additional, optional packages:
    gnuplot and Gnuplot.py (for plotting, http://gnuplot-py.sourceforge.net)
    matplotlib (alternate style of plotting, http://matplotlib.sourceforge.net)
    numpy (for some advanced features including checking and manipulating covariance matrices, http://numpy.scipy.org)

  [Optional] Setting Environment Variables: for general use of the fudge package, some changes should be made to your
computer's environment. The following lines make the required change. Note that <path_to_fudge> indicates the path to
the directory containing this README file:

  on unix with bash, ksh, etc, put this line in the .bashrc or equivalent:
    export PYTHONPATH=$PYTHONPATH:<path_to_fudge>

  on csh, tcsh or similar, do:
    setenv PYTHONPATH $PYTHONPATH:<path_to_fudge>

  on Windows, the environment variable should be added to the registry (see for example
       http://www.support.tabs3.com/main/R10463.htm)

  [Optional] Building Extensions: Fudge is primarily written in python, which is easily portable across many different
systems. Some advanced features of fudge (computationally intensive tasks) are however implemented in C, and
must be compiled before use. To build these extensions, use the Makefile included in this directory:
    >make

  If you encounter trouble with installing, building extensions or setting environment variables, please let us know
     (email: mattoon1@llnl.gov)


2) Basic Use:
  After installing Fudge, most users will be interested in translating ENDF files to the new format and in measuring how
well the translation works. Two tools for translating ENDF are included in the 'bin' directory:

  rePrint.py: translate one ENDF-formatted file to GND, then translate back to ENDF for comparison with the original.
Usage:
    <path_to_fudge>/bin/rePrint.py filename.endf
  Several files are produced:
    -test.endf6.xml	(the gnd-formatted file)
    -test.endf6-covar.xml  (covariances in gnd, produced only if the original file has covariances)
    -test.endf6.noLineNumbers  (the gnd file translated back to ENDF format)
    -test.endf6.orig.noLineNumbers  (the original ENDF file with line numbers stripped for easy comparison)

  After running rePrint.py, you may wish to compare the original/new ENDF files using diff,kompare, etc:
    >diff test.endf6.orig.noLineNumbers test.endf6.noLineNumbers

  See Section 5 below for a discussion of the most common differences between original/translated ENDF files.

  Another tool included in the release is 'rePrintSample.py'. This is very similar to rePrint, except it randomly picks 
an ENDF-formatted file from the included samples.

3) Getting Help:
  Documentation in html is packaged along with Fudge. To view the documentation, open the file "doc/html/index.html"
inside fudge. This offers extensive documentation of all the python code included in the package, and is meant to be
used as a reference.
For other questions about the package, please feel free to contact the maintainer at mattoon1@llnl.gov

