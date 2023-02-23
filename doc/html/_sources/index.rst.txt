######################
    Fudge and GNDS
######################

:Release: |version|
:Date: |today|

Welcome to FUDGE (For Updating Data and GEnerating Evaluations), a nuclear data platform built with Python that permits
reading, viewing, modifying, writing and processing of nuclear data. FUDGE is a product of the Nuclear Data and
Theory group (NDT) at Lawrence Livermore National Laboratory (LLNL).

For installation instruction be review the REAME.md file are visit https://github.com/llnl/fudge.

This release of Fudge expands support for processing and visualizing data in the Generalized Nuclear Data (GNDS) format.
GNDS is intended to modernize how nuclear data are stored and used, eventually replacing the legacy ENDF-6 format.
This version of FUDGE supports the latest version of the GNDS standard, version 2.0, as well as older versions '2.0.LLNL_4' and '1.10'.
The package includes tools to translate other formats to and from GNDS, plus tools
for testing, visualizing and processing GNDS-formatted data.
Python code within the *fudge* directory only support GNDS formatted data. Python code to support lecagy formats 
(e.g., ENDF-6 and ACE) reside in the brownies directory.

Fudge requires Python 3 (version 3.7 and higher).
See the installation notes for more details on working with multiple versions of Python.

New users (and returning users) are encouraged to check out :doc:`tutorial/index`.

If you just want to see what a GNDS file looks like, visit the ENDF/B-VIII library pages. The library
is available in both ENDF-6 and `GNDS format <http://www.nndc.bnl.gov/endf/b8.0/gndsfiles.html>`_.

Contents
========

.. toctree::
   :maxdepth: 1

   downloading
   installation
   gnds
   Tutorial <tutorial/index>
   externals
   reporting-bugs
   todo
   thanks
   The pqu package <pqu/modules.rst>

.. toctree::
   :hidden:
   
   The Fudge Package <fudge/index>
   glossary 

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`glossary`
* :ref:`search`
