#! /bin/env python

# <<BEGIN-copyright>>
# <<END-copyright>>

"""
Some simple stuff that is not repeatable as pointers are printed.
"""
import sys
sys.path.append( '../../../lib' )

import pointwiseXY_C

def setShow( f, x, y ) :
    print
    print 'Setting x, y to', x, y
    f.setValue( x, y )
    f.showInteralStructure( printPointersAsNull = True )

f = pointwiseXY_C.pointwiseXY_C( initialSize = 6, overflowSize = 3 )
f.showInteralStructure( printPointersAsNull = True )

setShow( f, 1., -1 )
setShow( f, 0.1, -2 )
setShow( f, 0.091, -1.2 )
setShow( f, 0.91, -1.3 )
setShow( f, 91, -1.6 )
setShow( f, 1e-3, -1.05 )
setShow( f, 1e-9, -1.15 )
setShow( f, 1e-6, -1.25 )
