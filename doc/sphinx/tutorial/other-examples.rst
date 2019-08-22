Other Fudge Examples
====================

After installing Fudge, most users will be interested in translating ENDF files to the new format and in measuring how
well the translation works. Two tools for translating ENDF are included in the ``bin`` directory.

rePrint.py 
----------
This is a tool translate one ENDF-formatted file to GNDS, then translate back to ENDF for comparison with the original.

:Usage: 
        ``<path_to_fudge>/bin/rePrint.py filename.endf``

Several files are produced:
    - ``test.endf6.xml``	(the gnds-formatted file)
    - ``test.endf6-covar.xml``  (covariances in gnds, produced only if the original file has covariances)
    - ``test.endf6.noLineNumbers``  (the gnds file translated back to ENDF format)
    - ``test.endf6.orig.noLineNumbers``  (the original ENDF file with line numbers stripped for easy comparison)

After running rePrint.py, you may wish to compare the original/new ENDF files using diff,kompare, etc::
    
    $diff test.endf6.orig.noLineNumbers test.endf6.noLineNumbers

Another tool included in the release is ``rePrintSample.py``. This is very similar to ``rePrint.py``, except it randomly picks 
an ENDF-formatted file from the included samples.

Making Plots
------------

After translating an ENDF file to GNDS, fudge can be used to make plots from the data. Plotting in fudge may require
installing extra components (``gnuplot`` and/or ``matplotlib``). A sample plotting script is included in the ``examples/``
directory::
  
    $ python <path_to_fudge>/examples/plotCrossSection.py <path_to_endf_file> <list_of_mt_numbers>
   
Another example can be used to compare a single cross section from two different evaluations. Evaluations may be in ENDF or GNDS format::

    $ python <path_to_fudge>/examples/compareCrossSections.py <mt> <first_file> <second_file>
