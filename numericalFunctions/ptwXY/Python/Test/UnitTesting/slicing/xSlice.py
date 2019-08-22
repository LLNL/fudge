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

verbose = False
if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__
    if( '-v' in options ) : verbose = True

def printXYIfVerbose( XYs, xMin = None, xMax = None, fill = None ) :

    if( verbose ) :
        if( not( xMin is None ) ) :
            print '# xMin =', xMin
            print '# xMax =', xMax
            print '# fill =', fill
        print '# length =', len( XYs )
        print XYs
        print

def checkSlicing( XYs_base, xMin, xMax, fill, length ) :

    XYs = XYs_base.xSlice( xMin = xMin, xMax = xMax, fill = fill )
    printXYIfVerbose( XYs, xMin, xMax, fill )
    if( len( XYs ) != length ) : raise Exception( "%s: len( XYs ) = %d != length = %d: %g %d %d" % \
        ( __file__, len( XYs ), length, xMin, xMax, fill ) )
    xMin_ = max( xMin, XYs_base.xMin( ) )
    xMax_ = min( xMax, XYs_base.xMax( ) )
    if( fill ) :
        if( xMin_ != XYs.xMin( ) ) : raise Exception( "%s: xMin_ = %g != XYs.xMin( ) = %g: %g %d %d" % ( __file__, xMin_, XYs.xMin( ), xMin, xMax, fill ) )
        if( xMax_ != XYs.xMax( ) ) : raise Exception( "%s: xMax_ = %g != XYs.xMax( ) = %g: %g %d %d" % ( __file__, xMax_, XYs.xMax( ), xMax, xMax, fill ) )
    else :
        if( xMin_ > XYs.xMin( ) ) : raise Exception( "%s: xMin_ = %g > XYs.xMin( ) = %g: %g %d %d" % ( __file__, xMin_, XYs.xMin( ), xMin, xMax, fill ) )
        if( xMax_ < XYs.xMax( ) ) : raise Exception( "%s: xMax_ = %g < XYs.xMax( ) = %g: %g %d %d" % ( __file__, xMax_, XYs.xMax( ), xMin, xMax, fill ) )

if( verbose ) :
    doc = pointwiseXY_C.pointwiseXY_C.xSlice.__doc__.split( '\n' )
    for d in doc : print "#", d

XYs = pointwiseXY_C.gaussian( .1, -10, xMax = 8, sigma = 2.2 )
printXYIfVerbose( XYs )

checkSlicing( XYs, -100, 100, 0, 61 )
checkSlicing( XYs, -100, 100, 1, 61 )

checkSlicing( XYs, -100, 3.3, 0, 43 )
checkSlicing( XYs, -100, 3.3, 1, 44 )

checkSlicing( XYs, -7, 3.3, 0, 19 )
checkSlicing( XYs, -7, 3.3, 1, 21 )

checkSlicing( XYs, 3.3, 100, 0, 18 )
checkSlicing( XYs, 3.3, 100, 1, 19 )
