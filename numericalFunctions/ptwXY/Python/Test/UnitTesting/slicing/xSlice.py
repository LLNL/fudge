# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import copy

from numericalFunctions import pointwiseXY_C

verbose = False
if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print( __file__ )
    if( '-v' in options ) : verbose = True

def printXYIfVerbose( XYs, xMin = None, xMax = None, fill = None ) :

    if( verbose ) :
        if( not( xMin is None ) ) :
            print( '# xMin =', xMin )
            print( '# xMax =', xMax )
            print( '# fill =', fill )
        print( '# length =', len( XYs ) )
        print( XYs )
        print( )

def checkSlicing( XYs_base, xMin, xMax, fill, length ) :

    XYs = XYs_base.domainSlice( domainMin = xMin, domainMax = xMax, fill = fill )
    printXYIfVerbose( XYs, xMin, xMax, fill )
    if( len( XYs ) != length ) : raise Exception( "%s: len( XYs ) = %d != length = %d: %g %d %d" % \
        ( __file__, len( XYs ), length, xMin, xMax, fill ) )
    xMin_ = max( xMin, XYs_base.domainMin( ) )
    xMax_ = min( xMax, XYs_base.domainMax( ) )
    if( fill ) :
        if( xMin_ != XYs.domainMin( ) ) : raise Exception( "%s: xMin_ = %g != XYs.domainMin( ) = %g: %g %d %d" % ( __file__, xMin_, XYs.domainMin( ), xMin, xMax, fill ) )
        if( xMax_ != XYs.domainMax( ) ) : raise Exception( "%s: xMax_ = %g != XYs.domainMax( ) = %g: %g %d %d" % ( __file__, xMax_, XYs.domainMax( ), xMax, xMax, fill ) )
    else :
        if( xMin_ > XYs.domainMin( ) ) : raise Exception( "%s: xMin_ = %g > XYs.domainMin( ) = %g: %g %d %d" % ( __file__, xMin_, XYs.domainMin( ), xMin, xMax, fill ) )
        if( xMax_ < XYs.domainMax( ) ) : raise Exception( "%s: xMax_ = %g < XYs.domainMax( ) = %g: %g %d %d" % ( __file__, xMax_, XYs.domainMax( ), xMin, xMax, fill ) )

if( verbose ) :
    doc = pointwiseXY_C.pointwiseXY_C.domainSlice.__doc__.split( '\n' )
    for d in doc : print( "#", d )

XYs = pointwiseXY_C.gaussian( .1, -10, domainMax = 8, sigma = 2.2 )
printXYIfVerbose( XYs )

checkSlicing( XYs, -100, 100, 0, 61 )
checkSlicing( XYs, -100, 100, 1, 61 )

checkSlicing( XYs, -100, 3.3, 0, 43 )
checkSlicing( XYs, -100, 3.3, 1, 44 )

checkSlicing( XYs, -7, 3.3, 0, 19 )
checkSlicing( XYs, -7, 3.3, 1, 21 )

checkSlicing( XYs, 3.3, 100, 0, 18 )
checkSlicing( XYs, 3.3, 100, 1, 19 )
