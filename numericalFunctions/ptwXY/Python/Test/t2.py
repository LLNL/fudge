#! /bin/env python

# <<BEGIN-copyright>>
# <<END-copyright>>

"""
Some simple stuff that should repeat.
"""
import sys
sys.path.append( '../../../lib' )

import pointwiseXY_C

data = [ [ 1, 1. ], [ 3., 3.14 ], [ 4.5, -5 ], [ 6.6, 2 ] ]
f = pointwiseXY_C.pointwiseXY_C( data, initialSize = 40, overflowSize = 10 )
print 'len( f ) =', len( f )
print 'f.allocatedSize( ) =', f.allocatedSize( )
print 'f.overflowAllocatedSize( ) =', f.overflowAllocatedSize( )
print 'f.overflowLength( ) =', f.overflowLength( )

print 'f[2] = ', f[2]

print 'f.getValue( 2. ) =', f.getValue( 2. )
print 'f.getValue( 5.5 ) =', f.getValue( 5.5 )
print 'f.getValue( 7. ) =', f.getValue( 7. )

print 'f.xMin( ) = ', f.xMin( )
print 'f.xMax( ) = ', f.xMax( )

print
print "printing each element of f"
for i in f : print "  ", i

print
print "printing f"
print f
