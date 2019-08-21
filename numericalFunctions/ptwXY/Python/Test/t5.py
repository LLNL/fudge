#! /bin/env python

# <<BEGIN-copyright>>
# <<END-copyright>>

"""
Some simple stuff that is not repeatable as pointers are printed.
"""
import sys
sys.path.append( '../../../lib' )

import random
import pointwiseXY_C

f = pointwiseXY_C.pointwiseXY_C( initialSize = 7, overflowSize = 3 )

xMin, xMax = -100, 10
yMin, yMax = -1., 0.
n = 100
i = 0
random.seed( 314159 )
while( i < n ) :
    r = random.random()
    x = xMin * r + xMax * ( 1. - r )
    r = random.random()
    y = yMin * r + yMax * ( 1. - r )
    print 'f.setValue( %.18e, %.18e )' % ( x, y )
    f.setValue( x, y )
    i += 1
    xm = xMin - 1
    for x, y in f :
        if( x < xm ) : raise Exception( 'x = %16e < xm = %16e: i = %d' % ( x, xm, i ) )
        xm = x

length = len( f )
f.showInteralStructure( printPointersAsNull = True )
