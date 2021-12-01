#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Some simple stuff that should repeat.
"""

from numericalFunctions import pointwiseXY_C

f = pointwiseXY_C.pointwiseXY_C( initialSize = 40, overflowSize = 10 )
print( 'len( f ) =', len( f ) )
print( 'f.allocatedSize( ) =', f.allocatedSize( ) )
print( 'f.overflowAllocatedSize( ) =', f.overflowAllocatedSize( ) )
print( 'f.overflowLength( ) =', f.overflowLength( ) )

data = [ [ 1, 1. ], [ 3., 3.14 ], [ 4.5, -5 ], [ 6.6, 2 ], [ 7.5, 4], [ 8.6, 2 ], [ 11.3, -5 ], [ 12.6, 3 ] ]
f.setData( data )
print( 'len( f ) =', len( f ) )

print( 'f[2] = ', f[2] )

print( 'f.evaluate( 2. ) =', f.evaluate( 2. ) )
print( 'f.evaluate( 5.5 ) =', f.evaluate( 5.5 ) )
print( 'f.evaluate( 7. ) =', f.evaluate( 7. ) )

print( 'f.domainMin( ) = ', f.domainMin( ) )
print( 'f.domainMax( ) = ', f.domainMax( ) )

print( )
print( "printing each element of f" )
for i in f : print( "  ", i )

print( f.toString.__doc__ )
print( )
print( "printing f" )
print( f )

print( 'printing f.toString( pairsPerLine = 3 )' )
print( f.toString( pairsPerLine = 3 ) )
print( """printing f.toString( format = " %6.3f %10.3e%%", pairsPerLine = 2, pairSeparator = ',' )""" )
print( f.toString( format = " %6.3f %10.3e%%", pairsPerLine = 2, pairSeparator = ',' ) )
print( """printing f.toString( format = " %6.3e %10.3E", pairSeparator = ',' )""" )
print( f.toString( format = " %6.3e %10.3E", pairSeparator = ',' ) )
