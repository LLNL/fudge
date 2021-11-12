######################
    Fudge and GNDS
######################

:Release: |version|
:Date: |today|

Welcome to FUDGE (For Updating Data and GEnerating Evaluations), a nuclear data platform built in Python that permits
reading, viewing, modifying, writing and processing of nuclear data. FUDGE is a product of the Computational Nuclear
Physics Group (CNP) at Lawrence Livermore National Lab (LLNL).

This release of Fudge expands support for processing and visualizing data in the Generalized Nuclear Data (GNDS) format.
GNDS is intended to modernize how nuclear data are stored and used, eventually replacing the legacy ENDF-6 format.
The package includes tools to translate other formats to and from GNDS, plus tools
for testing, visualizing and processing GNDS-formatted data.

As of this release, Fudge is compatible with Python 3 (version 3.6 and higher) as well as Python 2.7.
See the installation notes for more details on working with multiple versions of Python.

New users (and returning users) are encouraged to check out :doc:`tutorial/index`.

If you just want to see what a GNDS file looks like, visit the ENDF/B-VIII library pages. The library
is available in both ENDF-6 and `GNDS format <http://www.nndc.bnl.gov/endf/b8.0/gndsfiles.html>`_.


Contents
========

.. toctree::
   :maxdepth: 1

   The Fudge Package <fudge/index>
   externals
   gnds
   gidi
   Tutorial <tutorial/index>
   downloading
   installation
   reporting-bugs
   todo
   thanks
   history

.. toctree::
   :hidden:
   
   glossary 

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`glossary`
* :ref:`search`

