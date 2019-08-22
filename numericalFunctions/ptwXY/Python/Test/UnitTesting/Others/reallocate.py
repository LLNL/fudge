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

import os
import pointwiseXY_C

verbose = False

options = []
if( 'CHECKOPTIONS' in os.environ ) : options = os.environ['CHECKOPTIONS'].split( )
for argv in sys.argv[1:] : options += [ argv ]

if( '-e' in options ) : print __file__
if( '-v' in options ) : verbose = True

def printIfVerbose( ) :

    if( not( verbose ) ) : return
    print
    print '# length =' , len( data ), '  allocated size =', data.allocatedSize( )
    print '# overflow length =', data.overflowLength( ), '  overflow allocated size =', data.overflowAllocatedSize( )

if( verbose ) :
    doc = pointwiseXY_C.pointwiseXY_C.reallocatePoints.__doc__.split( '\n' )
    for d in doc : print "#", d
    print "#"
    doc = pointwiseXY_C.pointwiseXY_C.reallocateOverflowPoints.__doc__.split( '\n' )
    for d in doc : print "#", d

data = [ [ i, i * 10 + 3 ] for i in xrange( 51 ) ]
data = pointwiseXY_C.pointwiseXY_C( data, initialSize = 10, overflowSize = 10 )

printIfVerbose( )

data.reallocatePoints( 1000 )
data.reallocateOverflowPoints( 25 )
if( data.allocatedSize( ) != 1000 ) : raise Exception( '%s: data.allocatedSize( ) = %s != 1000' % ( __file__, data.allocatedSize( ) ) )
if( data.overflowAllocatedSize( ) != 25 ) : raise Exception( '%s: data.overflowAllocatedSize( ) = %s != 25' % i\
    ( __file__, data.overflowAllocatedSize( ) ) )

printIfVerbose( )

for i in xrange( 103 ) : data.setValue( i + 0.5, i * 10 + 3.5 )
printIfVerbose( )
