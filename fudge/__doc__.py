# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Fudge is a Python package which allows one to view, plot and modify GNDS nuclear data, and also process GNDS
into forms suitable for Monte Carlo and deterministic transport applications.
Processed data are stored as additional data 'forms' inside a GNDS file. The GIDIplus C++ API can then be used
to read and sample from the processed GNDS files.

Getting at fudge's python scripts.
==================================
To use fudge, one must add the location of the fudge scripts to the PYTHONPATH environment variable. (PYTHONPATH is the environment variable
used by python when searching for imported python modules (files). On LLNL's computing system (LC) this would look something like,

export PYTHONPATH=$PYTHONPATH:/usr/gapps/fudge4/current  OR  PYTHONPATH=$PYTHONPATH:/usr/gapps/fudge4/development

for the bash shell.

Alternatively, one can add the following lines near the top of a python script (or type them at the prompt)

    >>> import sys
    >>> sys.path.append( "/usr/gapps/fudge4/current" )

>>> import fudgeDefaults
>>> fudgeDefaults.NDF_DATABASE_DIR = /my/personal/database/processed

Reading fudge's documentation.
==============================
In general, one should first create a reactionSuite class object (also sometimes called a Protare for
Projectile + Target + Evaluation). For example the beginning of a fudge session may look like,

    >>> from fudge import reactionSuite
    >>> RS = reactionSuite.readXML("path/to/GNDS/reactionSuite.xml")

It is therefore important to read the documentation on the reactionSuite class.
For an overview of creating, visualizing and processing GNDS reactionSuites we recommend starting with
the documentation in doc/html/index.html (in particular the GNDS tutorial).
"""
