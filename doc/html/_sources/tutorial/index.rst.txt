######################
  The Fudge Tutorial
######################

:Release: |version|
:Date: |today|

The Generalized Nuclear Data (GND) language is a proposed successor to the venerable 
Evaluated Nuclear Data File (ENDF) format.  Currently, GND can be expressed in either 
the XML or HDF5 formats.  The ``fudge`` Python package provides tools to read, write, 
manipulate and process nuclear reaction evaluations in either the ENDF format or GND 
formats.  You can find more about ``fudge`` from the `GND Project Page 
<https://ndclx4.bnl.gov/gf/project/gnd/>`_ and about GND from `C. Mattoon, et al., 
Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145--3171 
<http://dx.doi.org/10.1016/j.nds.2012.11.008>`_.
 
In this tutorial, we will work through the development of several simple
fudge scripts that will: 

    * translate to/from the GND and ENDF nuclear reaction data formats;
    * check the physics content of an evaluation in either GND or ENDF format;
    * dig deeper into the GND hierarchy and the ``fudge`` code
    * outline a larger script for processing an evaluation for use 
      in a transport code;

I will assume throughout this tutorial that you have correctly installed fudge 
(see :ref:`installation` if you haven't), the Python scripting language (version 2.7.x) and the Python
extensions numpy and matplotlib.  I will also assume that you have working knowledge 
of the Python scripting language and that you are involved in nuclear data enterprise in 
some fashion (so I don't have to define things like "evaluation").

Let's begin with the most important thing.  Make sure fudge is in your PYTHONPATH!  
You can do this a few different ways:
    
    * Using the bash shell, at the command line do::
    
        $ export PYTHONPATH="/full/path/to/directory/containing/fudge"
        
    * Add this export directive to your ``.bashrc`` file.
    * Alternatively, you can amend the ``sys.path`` variable to point to the directory 
      containing the ``fudge`` python module *inside* any python script you write.    

OK, with that settled, you can start with the section on translating between ENDF and GND.

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

