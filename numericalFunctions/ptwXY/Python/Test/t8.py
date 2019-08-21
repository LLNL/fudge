#! /bin/env python

# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This stuff should be repeatable. Simple test of operators '+', '-', '*', 'neg', 'abs', 'pow' and union method.
"""
import sys
sys.path.append( '../../../lib' )

import random
import pointwiseXY_C

d = pointwiseXY_C.pointwiseXY_C( initialSize = 100, overflowSize = 31 )

d.setValue( 1, 11. )
d.setValue( 3.14, 3.12 )
d.setValue( 7.66, -1.5 )
d.setValue( 3.66, -1.6 )
d.setValue( 2.66,  1.7 )

print '---- d ----'
print d

e = d
del e[2]

print '---- del d[2] ----'
print d

d = pointwiseXY_C.pointwiseXY_C( initialSize = 100, overflowSize = 31 )

d.setValue( 1, 11. )
d.setValue( 3.14, 3.12 )
d.setValue( 7.66, -1.5 )
d.setValue( 3.66, -1.6 )
d.setValue( 2.66,  1.7 )
e[2] = d[2]

print '---- e ----'
print e

e = pointwiseXY_C.pointwiseXY_C( initialSize = 5, overflowSize = 3 )
xMin, xMax = 0.5, 8
yMin, yMax = -1.2, 10.
random.seed( 314159 )
for i in xrange( 7 ) :
    r = random.random( )
    x = xMin * r + xMax * ( 1. - r )
    r = random.random( )
    y = yMin * r + yMax * ( 1. - r )
    e.setValue( x, y )

print '---- d ----'
print d

print '---- neg of d ----'
print -d

print '---- abs of d ----'
print abs( d )

print '---- d.pow( 3 ) ----'
print d.__pow__( 3 )

print '---- e ----'
print e

print '---- d.union( e ) ----'
u = d.union( e )
print u

print '---- union and map e onto d ----'
u = d.union( e, True )
print u

print '---- union and map d onto e ----'
u = e.union( d, True )
print u

print '---- union and map and trim of e onto d ----'
u = d.union( e, True, True )
print u

e[ 0] = [ e[ 0][0], 0. ]        # needed to make e's domain mutual with d
e[-1] = [ e[-1][0], 0. ]

print '---- e ----'
print e

print '---- d + e ----'
print d + e

print '---- e + d ----'
print e + d

print '---- d + 1.41 ----'
print d + 1.41

print '----  1.41 + d ----'
print 1.41 + d

print '---- 1.41 - d ----'
print 1.41 - d

print '---- d - 1.41 ----'
print d - 1.41

print '---- d * 1.41 ----'
print d * 1.41

print '---- 1.41 * d ----'
print 1.41 * d

print '---- d + d ----'
print d + d

print '---- d - d ----'
print d - d

print '---- d * d ----'
print d * d

print '---- d * d * d * d ----'
print d * d * d * d

print '---- d**4 ----'
print d**4

print '---- d**-5 ----'
print d**-5
