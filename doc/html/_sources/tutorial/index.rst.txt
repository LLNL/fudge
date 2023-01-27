######################
  The Fudge Tutorial
######################

:Release: |version|
:Date: |today|

The Generalized Nuclear Database Structure (GNDS) has been accepted by the international
Nuclear Data community as the successor to the venerable
Evaluated Nuclear Data File (ENDF) format.  Currently, GNDS can be written out in either
the XML or HDF5 format (or in a hybrid mixture of the two).
The ``fudge`` Python package provides tools to read, write,
manipulate and process nuclear reaction evaluations in either the ENDF format or GNDS
formats.  You can find more about ``fudge`` from the `github site
<https://github.com/LLNL/fudge/>`_ and about GNDS from `C. Mattoon, et al.,
Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145--3171 
<http://dx.doi.org/10.1016/j.nds.2012.11.008>`_.
 
In this tutorial, we will work through the development of several simple
fudge scripts that will: 

    * translate to/from the GNDS and ENDF nuclear reaction data formats;
    * check the physics content of an evaluation in either GNDS or ENDF format;
    * dig deeper into the GNDS hierarchy and the ``fudge`` code
    * outline a larger script for processing an evaluation for use 
      in a transport code;

I will assume throughout this tutorial that you have correctly installed fudge 
(see :ref:`installation` if you haven't), the Python scripting language (version 3.7 and later)
and the Python extensions numpy and matplotlib.  I will also assume that you have working knowledge
of the Python scripting language and that you are involved in nuclear data enterprise in 
some fashion (so I don't have to define things like "evaluation").

Let's begin with the most important thing: if you installed FUDGE using the Makefile, make sure fudge is in your PYTHONPATH!
You can do this a few different ways, depending on what shell you use.

    * For bash (or ksh, zsh)::
        $ export PYTHONPATH="/full/path/to/directory/containing/fudge"  # only changes the current terminal session
        or add that line to .bashrc (or similar)

    * For csh (or tcsh)::
        $ setenv PYTHONPATH "/full/path/to/directory/containing/fudge"
        or add that line to .cshrc or .tcshrc
    
    * Alternatively, you can amend the ``sys.path`` variable to point to the directory 
      containing the ``fudge`` python module *inside* any python script you write.::
      $ import sys
      $ sys.path.append("/full/path/to/directory/containing/fudge")


You can skip this step if you installed using pip, although you may need to activate a Python environment.
With that settled, you can start with the section on translating between ENDF and GNDS.

.. note::  If you find yourself lost, the Python interpretor has excellent built-in help 
           (just type ``help( object )`` to learn about an object.  Then, of course, there 
           is the Fudge Documentation itself which should have been included with your 
           Fudge installation in ``doc/html/index.html``.  

**Contents:**

.. toctree::
   :maxdepth: 1
   :numbered:
   
   translating
   checking
   navigating
   diggingDeeper
   processing
   other-examples

.. commented out stuff
    Indices and tables
    ==================
    
    * :ref:`genindex`
    * :ref:`search`

