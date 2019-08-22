Checking ENDF or GND files
==========================

Let's face it.  ENDF evaluations often have mistakes and bugs.  Fortunately `fudge` comes 
armed with a large number of pre-made checks.  Running them is as simple as running the 
``check()`` member function on any fudge class instance.  The quickest way to test an 
evaluation and its covariance is just to run the main ``check()`` function of either an 
instance of ``reactionSuite`` or ``covarianceSuite``.  ``Fudge`` will then recursively descend
through the instances and check every data member.

Checking ``reactionSuite`` instances
------------------------------------

Using the ``myEval`` and ``myCov`` instances from the previous tutorial, we find:

    >>> myEval.check()
    n + H1
    H2 + photon
    <fudge.gnd.warning.context object at 0x101cb6190>
    
Hmm... This doesn't look so helpful.  This is because ``fudge`` understand the hierarchy 
and categorizes the errors and warnings accordingly.  Try this:

    >>> print myEval.check()
    ReactionSuite: n + H1
    WARNING: Wick's limit too low by 3.975% at 1000000.0
    WARNING: Wick's limit too low by 11.792% at 1200000.0
    WARNING: Wick's limit too low by 17.839% at 1400000.0
    WARNING: Wick's limit too low by 22.677% at 1600000.0
    WARNING: Wick's limit too low by 26.645% at 1800000.0
    WARNING: Wick's limit too low by 29.961% at 2000000.0
    WARNING: Wick's limit too low by 32.778% at 2200000.0
    WARNING: Wick's limit too low by 35.199% at 2400000.0
    WARNING: Wick's limit too low by 37.305% at 2600000.0
    WARNING: Wick's limit too low by 39.151% at 2800000.0
    WARNING: Wick's limit too low by 40.784% at 3000000.0
    WARNING: Wick's limit too low by 42.238% at 3200000.0
    WARNING: Wick's limit too low by 43.540% at 3400000.0
    WARNING: Wick's limit too low by 44.713% at 3600000.0
    WARNING: Wick's limit too low by 45.774% at 3800000.0
    WARNING: Wick's limit too low by 46.738% at 4000000.0
    WARNING: Wick's limit too low by 47.618% at 4200000.0
    WARNING: Wick's limit too low by 48.424% at 4400000.0
    WARNING: Wick's limit too low by 49.164% at 4600000.0
    WARNING: Wick's limit too low by 49.845% at 4800000.0
    WARNING: Wick's limit too low by 50.475% at 5000000.0
    WARNING: Wick's limit too low by 51.855% at 5500000.0
    WARNING: Wick's limit too low by 53.007% at 6000000.0
    WARNING: Wick's limit too low by 53.979% at 6500000.0
    WARNING: Wick's limit too low by 54.807% at 7000000.0
    WARNING: Wick's limit too low by 55.514% at 7500000.0
    WARNING: Wick's limit too low by 56.123% at 8000000.0
    WARNING: Wick's limit too low by 56.647% at 8500000.0
    WARNING: Wick's limit too low by 57.099% at 9000000.0
    WARNING: Wick's limit too low by 57.489% at 9500000.0
    WARNING: Wick's limit too low by 57.824% at 10000000.0
    WARNING: Wick's limit too low by 58.111% at 10500000.0
    WARNING: Wick's limit too low by 58.354% at 11000000.0
    WARNING: Wick's limit too low by 58.559% at 11500000.0
    WARNING: Wick's limit too low by 58.728% at 12000000.0
    WARNING: Wick's limit too low by 58.865% at 12500000.0
    WARNING: Wick's limit too low by 58.973% at 13000000.0
    WARNING: Wick's limit too low by 59.053% at 13500000.0
    WARNING: Wick's limit too low by 59.107% at 14000000.0
    WARNING: Wick's limit too low by 59.138% at 14500000.0
    WARNING: Wick's limit too low by 59.146% at 15000000.0
    WARNING: Wick's limit too low by 59.133% at 15500000.0
    WARNING: Wick's limit too low by 59.101% at 16000000.0
    WARNING: Wick's limit too low by 59.049% at 16500000.0
    WARNING: Wick's limit too low by 58.980% at 17000000.0
    WARNING: Wick's limit too low by 58.894% at 17500000.0
    WARNING: Wick's limit too low by 58.792% at 18000000.0
    WARNING: Wick's limit too low by 58.674% at 18500000.0
    WARNING: Wick's limit too low by 58.542% at 19000000.0
    WARNING: Wick's limit too low by 58.396% at 19500000.0
    WARNING: Wick's limit too low by 58.236% at 20000000.0
    reaction label 1: H2 + gamma
        WARNING: Calculated and tabulated Q-values disagree: 2735357.48694634 eV vs 2224631 eV!

Much more useful.  In fact, you can loop over the warnings to look for specific ones:

    >>> for w in myEval.check():
    ...     if 'Q-value' in str( w ): print w
    ...
    reaction label 1: H2 + gamma
        WARNING: Calculated and tabulated Q-values disagree: 2735357.48694634 eV vs 2224631 eV!

In other words, the tabulated Q-value for the capture reaction (MT 102) is inconsistent with the
value calculated from incident and outgoing particle masses.
Note, you don't just get the warning or error, you get the entire context in which it occurred.  
This is very useful for bug tracking. In the example above, we can tell that the problem is in the capture
reaction by looking at the list of outgoing particles: 'H2 + gamma'.

Checking ``covarianceSuite`` instances
--------------------------------------

The covariances have ``check()`` functions too and they find all sorts of stuff:

    >>> myCov.check()
    CovarianceSuite: n + H1
        Section: H2 + gamma
            Form covarianceMatrix:
                WARNING: Ratio of smallest/largest eigenvalue (2.359558e-10) is too small

Setting up a checker script
---------------------------

We will set up a script to use ``fudge`` to check ENDF files.  I'll leave to you to figure out how
to do the same with GND files (if you do, remember to load both the evaluation and the covariance!).

This is what I came up with (download it :download:`here <checkendf.py>`):

::

    #! /usr/bin/env python
    import argparse
    from fudge.legacy.converting.endfFileToGND import endfFileToGND
    
    # Process command line options
    parser = argparse.ArgumentParser(description='Check an ENDF file')
    parser.add_argument('inFile', type=str, help='The ENDF file you want to translate and check.' )
    args = parser.parse_args()
    
    # Now translate
    results = endfFileToGND( args.inFile, toStdOut=True, skipBadData=True )
    myEval = results['reactionSuite']
    myCov = results['covarianceSuite']
    print '\n\n'
    
    # Check the evaluation
    print "Checking evaluation for "+args.inFile
    print "------------------------------------------------"
    print myEval.check()
    
    print '\n'
    
    # Check the covariance
    print "Checking covariances for "+args.inFile
    print "------------------------------------------------"
    print myCov.check()

Try it out!
