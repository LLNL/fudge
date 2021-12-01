#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Some simple stuff that is not repeatable as pointers are printed.
"""

from numericalFunctions import pointwiseXY_C

def setShow( f, x, y ) :
    print( )
    print( 'Setting x, y to', x, y )
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
