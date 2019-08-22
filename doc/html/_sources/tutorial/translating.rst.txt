Basic translating between ENDF and GNDS
======================================

Fudge comes with a tool (fudge/bin/rePrint.py) for translating ENDF-6 evaluations to and from GNDS.
For the purpose of demonstrating how Fudge works, in this tutorial we will reproduce the capability of rePrint.py
by writing two codes:
    
    * ``endf2gnds.py`` : a code to translate ENDF formatted data into the GNDS format
    * ``gnds2endf.py`` : a code to translate GNDS formatted data into ENDF

In both scripts (and all the others we will write), we will use the Python `argparse` module to 
manage the command line options of the scripts.  Before we get to that, lets just start by learning how to translate between ENDF and GNDS.

Translating from ENDF into GNDS
------------------------------
Let's experiment with fudge in the Python shell.

::

    $ python
    Python 2.7.1 (r271:86832, Jun 16 2011, 16:59:05) 
    [GCC 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2335.15.00)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.

    
First, import the ``endfFileToGNDS`` function from ``fudge.legacy.converting``:

    >>> from fudge.legacy.converting.endfFileToGNDS import endfFileToGNDS
    
I assume that you have an ENDF formatted file somewhere and for the sake of argument, I'll 
assume that it is ``n-001_H_001.endf`` and that it is in the current directory, ``'.'``.  Given this, 
simply do:

    >>> resultsDict = endfFileToGNDS( 'n-001_H_001.endf', toStdOut=True, skipBadData=True )
    
Before hitting return, note a few things.  First, the ``endfFileToGNDS()`` function takes one mandatory argument,
the ENDF file name, and returns a Python dictionary.
The other arguments ``toStdOut`` and ``skipBadData`` do not need to be specified.
The ``skipBadData`` option tells ``fudge`` to keep going when reading an ENDF file, even if it is "broken".
``Fudge`` is picky and will otherwise throw an exception when it encounters bad or malformed data
(an all-too-common occurrance sadly).  The ``toStdOut`` option tells ``fudge`` to redirect all warnings and errors
to ``stdout`` rather than ``stderr``.

So, let's hit return and see what happens:

    >>> resultsDict = endfFileToGNDS( 'n-001_H_001.endf', toStdOut=True, skipBadData=True )
      2 [3, 4, 33] : MF=4, LTT = 1
    102 [3, 6, 33] : MF=6 : ZAP=0, LAW=2, LANG=0 : ZAP=1002, LAW=4
      1 [3, 33]
    Reading resonances (MF=2 MT=151)
    Reading covariances (MFs 33)

Voila, your ENDF file is translated into GNDS.  The first number listed is the ENDF reaction designator
MT (2, 102 and 1 are present in this file).  The list following the MT is the list of observables
for which there is data in this reaction, i.e. the list of MF's.  The final set of values denote any special formats for those MF's.

The dictionary returned by endfFileToGNDS contains several items:

    >>> resultsDict.keys()
['info', 'reactionSuite', 'errors', 'covarianceSuite']

resultsDict['reactionSuite'] is an instance of the Fudge ``reactionSuite`` class containing the entire `GNDS` file.
In other words, endfFileToGNDS reads and translates the data at the same time.
The dictionary also contains a ``covarianceSuite`` instance with covariances from the evaluation.
If the original ENDF file has no covariances, then the value will be set to ``None``.
The dictionary also contains

To write translated results to a file, we need to use the ``toXMLList()`` member functions
of the ``reactionSuite`` and ``covarianceSuite`` classes:
    
    >>> open( gnds, mode='w' ).writelines(
            resultsDict['reactionSuite'].toXMLList(verbosity=0) )
    >>> if myCov is not None:
    >>>     open( "n-001_H_001.gndsCov.xml", mode='w' ).writelines( resultsDict['covarianceSuite'].toXMLList(verbosity=0) )

In both cases, we may hand the ``toXMLList()`` member function extra processing directives
(in this case, we simply set 'verbosity' to 0 to indicate that we want the output to be quiet).

Translating GNDS files back into ENDF
------------------------------------

Now let us reconstruct the original ENDF file.  First we need to import an additional module that supports writing
GNDS back to ENDF:

    >>> from site_packages.legacy.toENDF6 import toENDF6

After this import, all ``reactionSuite`` instances contain a member function ``toENDF6()``.  Let's try it out:

    >>> evalStyle = resultsDict['info'].style
    >>> myENDF = resultsDict['reactionSuite'].toENDF6( evalStyle, {'verbosity':0},
            covarianceSuite=resultsDict['covarianceSuite'], lineNumbers = False )
    >>> open( "junk.endf", mode='w' ).write( myENDF )

A few comments on the previous lines: the 'evalStyle' is needed when going back to ENDF-6 because a GNDS file
can contain more than one style of data, including one or more evaluations along with processed data.
Only one style (indicated by a string label) is translated back to ENDF-6.  The ``{'verbosity':0}`` argument
contains extra flags used when writing back to ENDF-6. We also supply a ``covarianceSuite`` so covariances are
written back, and tell Fudge not to write line numbers.

Now you should quit Python (using ^d), and check out what you made.

::

    $ diff junk.endf n-001_H_001.endf 
    1c1
    <                                                                      1 0  0
    ---
    >  $Rev:: 532      $  $Date:: 2011-12-05#$                             1 0  0    0
    93c93
    <                                 1        451        101          0 125 1451   92
    ---
    >                                 1        451        101          5 125 1451   92
    95,102c95,102
    <                                 3          1         35          0 125 1451   94
    <                                 3          2         35          0 125 1451   95
    <                                 3        102         35          0 125 1451   96
    <                                 4          2        196          0 125 1451   97
    <                                 6        102        201          0 125 1451   98
    <                                33          1          5          0 125 1451   99
    <                                33          2         21          0 125 1451  100
    <                                33        102         21          0 125 1451  101
    ---
    >                                 3          1         35          4 125 1451   94
    >                                 3          2         35          4 125 1451   95
    >                                 3        102         35          5 125 1451   96
    >                                 4          2        196          4 125 1451   97
    >                                 6        102        201          4 125 1451   98
    >                                33          1          5          5 125 1451   99
    >                                33          2         21          5 125 1451  100
    >                                33        102         21          5 125 1451  101
    113c113
    <          30          5         96          2          0          0 125 3  1    3
    ---
    >          30          5         96          2                       125 3  1    3
    149c149
    <          96          2          0          0          0          0 125 3  2    3
    ---
    >          96          2                                             125 3  2    3
    185c185
    <          30          5         96          2          0          0 125 3102    3
    ---
    >          30          5         96          2                       125 3102    3
    223c223
    <          96          2          0          0          0          0 125 4  2    4
    ---
    >          96          2                                             125 4  2    4
    420c420
    <           2          2          0          0          0          0 125 6102    3
    ---
    >           2          2                                             125 6102    3
    423c423
    <          96          2          0          0          0          0 125 6102    6
    ---
    >          96          2                                             125 6102    6
    617c617
    <           2          2          0          0          0          0 125 6102  200
    ---
    >           2          2                                             125 6102  200    

Not bad...  There are obviously several differences.  Let's examine them:

**Line 1:**
      The ``$Rev::$`` and ``$Date::`` fields are put in by the NNDC on the 
      very first line of every ENDF file simply to enable subversion version control
      keyword substitutions.  This line is not part of the ENDF standard and may be 
      safely ignored.
**Lines 92-101:**
      These lines are the ENDF dictionary in the end of the free text discriptive
      section (MF1/MT451).  The only difference here is that the ENDF section version numbers 
      were are set to 0.  In this case, this messes up the versioning of ``n-001_H_001.endf``, 
      however we note that few evaluators remember to set these values in practice.
**Remainder of lines:**
      In each case, the original ENDF file did not quite follow the ENDF format
      strictly and entered empty strings where the integer ``0`` should have been used.

When translating from ENDF, you may notice some substantial differences between the original and re-translated file.
Some differences are due to sections that are not yet translated to the new format (for example, delayed gammas from ENDF
MF 1 MT 460 are not yet translated). Other differences include:

    - masses, which appear many times in ENDF and are often inconsistent. In a GNDS file, the mass is stored only once,
      so upon translation back to ENDF inconsistent masses are overwritten.

    - duplicate points: ENDF files sometimes contain two or more duplicate X-Y pairs in a cross section, multiplicity
      or distribution. Unless these appear at the boundary between interpolation regions, the ENDF-to-GNDS translator
      drops the second point as unnecessary, leading to differences when comparing to the original ENDF file.

    - interpolation regions: ENDF files permit using different interpolation (lin-lin, log-lin, etc) in different
      regions. GNDS also supports this, but where possible we have merged two or more regions into a single region (for
      example, 'flat' interpolation regions can be merged with lin-lin regions with no loss of precision). Also, ENDF
      files may contain discontinuous functions within a single interpolation region. Upon translating to GNDS, these are
      converted into multiple regions.

    - Improperly-formatted ENDF files: the GNDS translation tool strictly interprets the ENDF format as defined in the
      June 2010 version of the ENDF manual (available at https://ndclx4.bnl.gov/gf/project/endf6man). Some differences come
      from files in the ENDF library that do not strictly follow the format. As a common example, some ENDF files contain
      non-zero data in a reserved field. After translation, the entry is reset to '0'.



Reading GNDS XML files
---------------------

If I didn't have pre-made instances of ``reactionSuite`` and ``covarianceSuite``, how would I read in the XML files?
For this purpose, both the ``fudge.gnds.reactionSuite`` and ``fudge.gnds.covariances`` have the factory function ``readXML()``.
To use them do:

    >>> from fudge.gnds import reactionSuite
    >>> from fudge.gnds.covariances import covarianceSuite
    >>> myOtherEval = reactionSuite.readXML( "n-001_H_001.gnds.xml" )

This reads in the evaluation itself.  To read in the covariances, we need to tell the
To use them do:

    >>> from fudge.gnds import reactionSuite
    >>> from fudge.gnds.covariances import covarianceSuite
    >>> myOtherEval = reactionSuite.readXML( "n-001_H_001.gnds.xml" )

This reads in the evaluation itself.  To read in the covariances, we need to tell the
To use them do:

    >>> from fudge.gnds import reactionSuite
    >>> from fudge.gnds.covariances import covarianceSuite
    >>> myOtherEval = reactionSuite.readXML( "n-001_H_001.gnds.xml" )

This reads in the evaluation itself.  To read in the covariances, we need to tell the `covariances.readXML()` function
where the evaluation is so that it can set up the hyperlinks correctly:

    >>> myOtherCov = covarianceSuite.readXML( "n-001_H_001.gndsCov.xml", reactionSuite=myOtherEval )

Setting up the translator scripts
---------------------------------

In this final section of the first tutorial, we'll actually make the two scripts ``endf2gnds.py`` and ``gnds2endf.py``.
Let's start by making the files and then editing the first:
::

    $ touch endf2gnds.py gnds2endf.py
    $ chmod u+x endf2gnds.py gnds2endf.py
    $ vim endf2gnds.py
    
For ``endf2gnds.py``, we want to read one ENDF file and write the GNDS evaluation file and (if present) the GNDS covariance file. Since there are two output files, we want them to have the same prefix for bookkeeping purposes.  So, here is my version of ``endf2gnds.py`` (download it :download:`here <endf2gnds.py>`):
::

    #! /usr/bin/env python
    import argparse
    from fudge.legacy.converting.endfFileToGNDS import endfFileToGNDS
    
    # Process command line options
    parser = argparse.ArgumentParser(description='Translate ENDF into GNDS')
    parser.add_argument('inFile', type=str, help='The ENDF file you want to translate.' )
    parser.add_argument('-o', dest='outFilePrefix', default=None, help='''Specify the output file's prefix to be ``outFilePrefix``.  The outputted files have extensions ".gnds.xml" and ".gndsCov.xml"vfor the GNDS main evaluations and covariance files and ".endf" for ENDF files.''' )
    args = parser.parse_args()
    
    # Compute output file names
    if args.outFilePrefix != None:
        outEvalFile = args.outFilePrefix + '.gnds.xml'
        outCovFile = args.outFilePrefix + '.gndsCov.xml'
    else:
        outEvalFile = args.inFile.replace( '.endf', '.gnds.xml' )
        outCovFile = args.inFile.replace( '.endf', '.gndsCov.xml' )
        
    # Now translate
    results = endfFileToGNDS( args.inFile, toStdOut=True, skipBadData=True )
    myEval = results['reactionSuite']
    myCov = results['covarianceSuite']
    open( outEvalFile, mode='w' ).writelines( line+'\n' for line in myEval.toXMLList( {'verbosity':0} ) )
    if myCov is not None:
         open( outCovFile, mode='w' ).writelines( line+'\n' for line in myCov.toXMLList( {'verbosity':0} ) )

I urge you to try it out.  If you are unsure how to use it, type ``./endf2gnds.py --help``.

``gnds2endf.py`` is similar.  However, we need to specify an input file prefix, the style to translate,
and optionally the  output file name.  This is my version of ``gnds2endf.py`` (download it :download:`here <gnds2endf.py>`):
::

    #! /usr/bin/env python
    import argparse, os
    from fudge.gnds import reactionSuite
    from fudge.gnds.covariances import covarianceSuite
    
    # Process command line options
    parser = argparse.ArgumentParser(description='Translate GNDS into ENDF')
    parser.add_argument('inFilePrefix', type=str, help='The prefix of the GNDS files you want to translate.' )
    parser.add_argument('--style', type=str, help='Data style to translate back to ENDF-6', default='eval' )
    parser.add_argument('-o', dest='outFile', default=None, help='Specify the output file' )
    args = parser.parse_args()
    
    # Compute input file names
    inEvalFile = args.inFilePrefix + '.gnds.xml'
    inCovFile = args.inFilePrefix + '.gndsCov.xml'
    
    # Compute the output file name
    if args.outFile == None: outFile = args.inFilePrefix + '.endf'
    else:                    outFile = args.outFile
        
    # Read in XML files
    myEval = reactionSuite.readXML( inEvalFile )
    if os.path.exists( inCovFile ): myCov = covarianceSuite.readXML( inCovFile, reactionSuite=myEval )
    else:                           myCov = None
    
    # Now translate
    open( outFile, mode='w' ).write( myEval.toENDF6( args.style, {'verbosity':0}, covarianceSuite=myCov ) )

