#! /bin/env python

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

import sys
sys.path.insert( 0, '../../../../../lib' )

import os, copy
import pointwiseXY_C

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__

def checkSlicing( XYs_base, i1, i2, doPrint = False ) :

    XYs = copy.copy( XYs_base )
    pXYs = pointwiseXY_C.pointwiseXY_C( data = XYs, initialSize = 20 )
    sXYs = XYs[i1:i2]
    if( doPrint ) :
        print
        print i1, i2
        for xy in XYs : print xy
        print
        print sXYs
    spXYs = pXYs[i1:i2]
    if( doPrint ) : spXYs
    if( len( sXYs ) != len( spXYs ) ) : raise Exception( "%s: len( sXYs ) = %d != len( spXYs ) = %d: index1 = %d, i2 = %d" % \
        ( __file__, len( sXYs ), len( spXYs ), i1, i2 ) )
    for i, xy in enumerate( sXYs ) :
        if( xy[0] != spXYs[i][0] ) : raise Exception( "%s: difference at index = %d: %e %e" % ( __file__, i, xy[0], spXYs[i][0] ) )
    
XYs = [ [ float( x ), float( x )**2 ] for x in xrange( 12 ) ]

checkSlicing( XYs, 4, 8 )
checkSlicing( XYs, -4, 8 )
checkSlicing( XYs, -4, -2 )
checkSlicing( XYs, 4, -2 )
checkSlicing( XYs, 4, -8 )
checkSlicing( XYs, -4 - 7 * len( XYs ), 8 )
checkSlicing( XYs, 3, -4 - 7 * len( XYs ) )
