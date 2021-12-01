# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import sys

from numericalFunctions import pointwiseXY_C

verbose = False

options = []
if( 'CHECKOPTIONS' in os.environ ) : options = os.environ['CHECKOPTIONS'].split( )
for argv in sys.argv[1:] : options += [ argv ]

if( '-e' in options ) : print( __file__ )
if( '-v' in options ) : verbose = True

def printIfVerbose( ) :

    if( not( verbose ) ) : return
    print( )
    print( '# length =' , len( data ), '  allocated size =', data.allocatedSize( ) )
    print( '# overflow length =', data.overflowLength( ), '  overflow allocated size =', data.overflowAllocatedSize( ) )

if( verbose ) :
    doc = pointwiseXY_C.pointwiseXY_C.reallocatePoints.__doc__.split( '\n' )
    for d in doc : print( "#", d )
    print( "#" )
    doc = pointwiseXY_C.pointwiseXY_C.reallocateOverflowPoints.__doc__.split( '\n' )
    for d in doc : print( "#", d )

data = [ [ i, i * 10 + 3 ] for i in range( 51 ) ]
data = pointwiseXY_C.pointwiseXY_C( data, initialSize = 10, overflowSize = 10 )

printIfVerbose( )

data.reallocatePoints( 1000 )
data.reallocateOverflowPoints( 25 )
if( data.allocatedSize( ) != 1000 ) : raise Exception( '%s: data.allocatedSize( ) = %s != 1000' % ( __file__, data.allocatedSize( ) ) )
if( data.overflowAllocatedSize( ) != 25 ) : raise Exception( '%s: data.overflowAllocatedSize( ) = %s != 25' % i\
    ( __file__, data.overflowAllocatedSize( ) ) )

printIfVerbose( )

for i in range( 103 ) : data.setValue( i + 0.5, i * 10 + 3.5 )
printIfVerbose( )
