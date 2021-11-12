# FUDGE (For Updating Data and Generating Evaluations)
FUDGE is a nuclear data platform built in Python that permits reading, viewing, modifying,
writing and processing of nuclear data. FUDGE is a product of the Computational Nuclear
Physics Group (CNP) at Lawrence Livermore National Lab (LLNL).

The current release of Fudge focuses on extending support for v2.0 of the Generalized Nuclear Database Structure (GNDS) format,
intended to modernize storage of nuclear data. The package includes tools to translate other formats to and from GNDS,
plus tools for testing, visualizing, and processing GNDS-formatted data.

# Contents of this README:
    - 1) Getting Started
    - 2) Basic Use (translating ENDF files)
    - 3) Getting Help

## Getting Started:
  For basic use the package only requires Python (version 3.6 or higher) and NumPy (version 1.15 or higher).
  The following installations options are currently supported:

    - Installation via `pip install`: This is ideal for use in virtual environments and the following steps are recommended:

        (i)  Ensure that NumPy (version 1.15 or higher) and wheels are installed in your Python environment

        (ii) Run the command `pip install git+ssh://git@czgitlab.llnl.gov:7999/nuclear/fudge/fudge.git@fudge4.3-rc5`

    - Installation via Makefiles: This is the typical development mode for active FUDGE code maintenance and improvements.
      The following steps are recommended:

        (i)   Ensure that NumPy (version 1.15 or higher) is installed

        (ii)  Use the following command to clone FUDGE to the current folder: 
              `git clone --recurse-submodules ssh://git@czgitlab.llnl.gov:7999/nuclear/fudge/fudge.git`

        (iii) Build FUDGE with the following command: `cd fudge; make -s`

        (iv)  [Optional] Setting Environment Variables: for general use of the FUDGE package, some changes should be made to your
              computer's environment. The following lines make the required change. Note that <path_to_FUDGE> indicates the path to
              the directory containing this README file:

              on unix with bash, ksh, etc, put this line in the .bashrc or equivalent:
                export PYTHONPATH=$PYTHONPATH:<path_to_FUDGE>

              on csh, tcsh or similar, do:
                setenv PYTHONPATH $PYTHONPATH:<path_to_FUDGE>

              on Windows, the environment variable should be added to the registry (see for example
                   http://www.support.tabs3.com/main/R10463.htm)

## Requirements for interactive plotting
   With the change from Python2 to Python3, FUDGE is replacing its Gnuplot-based interactive plotting with a PyQt5 based
   one. While PyQt5 is not a requirement for the general use of FUDGE, it is needed to use the new FUDGE interactive 
   plotting module  

## Basic Use:
  After installing FUDGE, most users will be interested in translating ENDF files to the new format and in measuring how
well the translation works. Two tools for translating ENDF are included in the 'brownies/bin' directory:

  rePrint.py: translate one ENDF-formatted file to GNDS, then translate back to ENDF for comparison with the original.
Usage:
    <path_to_FUDGE>/brownies/bin/rePrint.py filename.endf -v
  Several files are produced:
    -test.endf6.xml	(the gnds-formatted file)
    -test.endf6-covar.xml  (covariances in gnds, produced only if the original file has covariances)
    -test.endf6.noLineNumbers  (the gnds file translated back to ENDF format)
    -test.endf6.orig.noLineNumbers  (the original ENDF file with line numbers stripped for easy comparison)

  After running rePrint.py, you may wish to compare the original/new ENDF files using diff,kompare, etc:
    >diff test.endf6.orig.noLineNumbers test.endf6.noLineNumbers

  endf2gnds.py: one-way translation from ENDF to GNDS. Useful for translating multiple ENDF-6 files since the GNDS versions will
  follow the same naming convention. For example, n-001_H_001.endf produces n-001_H_001.endf.gnds.xml (and n-001_H_001.endf.gnds-covar.xml)
Usage:
  <path_to_FUDGE>/brownies/bin/endf2gnds.py filename.endf -v

## Getting Help:
  The FUDGE documentation is currently undergoing an update but the historical documentation in html is packaged along with FUDGE for 
deployments via the Makefile option. To view the documentation, open the file "doc/html/index.html" inside FUDGE. This offers extensive 
documentation of all the Python code included in the package, and is meant to be used as a reference.

For other questions about the package, please feel free to contact the maintainer at mattoon1@llnl.gov

# License
FUDGE is distributed under the terms of the 3-clause BSD license.

All new contributions must be made under the 3-clause BSD license.

See LICENSE, COPYRIGHT and NOTICE for details.

SPDX-Licence-Identifier: BSD-3-Clause

This package includes several components: \
LLNL-CODE-683960	(FUDGE)

LLNL-CODE-770134	(numericalFunctions \
LLNL-CODE-771182	(statusMessageReporting) \
LLNL-CODE-725546	(Merced)

