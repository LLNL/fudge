#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
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
