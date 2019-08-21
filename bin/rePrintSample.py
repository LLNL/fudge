#! /usr/bin/env python

# <<BEGIN-copyright>>
# <<END-copyright>>

import sys, os, glob, random

binDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert( 0, os.path.dirname( binDir ) )
from fudge.core.utilities.brb import Pause
from fudge.core.utilities import subprocessing

endfPath = "../examples"
if( len( sys.argv ) > 1 ) : endfPath = sys.argv[1]

files = glob.glob( os.path.join(endfPath,"*.endf") )
if( not files ) :
    print "No ENDF files were found in %s. Please give the path to a directory containing ENDF files." % endfPath
    sys.exit( 1 )

random.shuffle( files )
for index, file in enumerate( files ) :
    returncode, stdout, stderr = subprocessing.executeCommand( [ 'python', os.path.join( binDir, 'rePrint.py' ), file ], raiseOnError = False )
    print file
    try:
        returncode, stdout, stderr = subprocessing.executeCommand( [ 'xdiff', '--geometry', '1500x1150+10+10', 
            'test.endf6.orig.noLineNumbers.cleanAndFixed', 'test.endf6.noLineNumbers' ], raiseOnError = False )
    except:
        print """I couldn't find the xdiff command on your system! If you have another tool (like kompare), try the following:
> kompare test.endf6.orig.noLineNumbers.cleanAndFixed test.endf6.noLineNumbers
        """
    response = Pause( 'Enter <RET> to continue or "q" to quit [%d of %d sampled] ' % ( index + 1, len( files ) ) )
    if( response == 'q' ) : break
