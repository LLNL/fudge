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

import sys, os
sys.path.insert( 0, '../../../../../lib' )

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__

import pointwiseXY_C

value = 2.4
orig = pointwiseXY_C.pointwiseXY_C( [ [ 0, 1 ], [ 1, 3 ], [ 3, -1 ] ] )
errs = []

def check( operator, other ) :

    global orig
    m1 = eval( 'orig %s other' % operator )
    m2 = orig.copy( )
    m3 = m2
    exec( 'm2 %s= other' % operator )
    d1 = m2 - m1
    if( d1.yMin( ) != 0. ) : errs.append( 'failure for operator "%s": d1.yMin( ) = %e != 0.' % ( operator, d1.yMin( ) ) )
    if( d1.yMax( ) != 0. ) : errs.append( 'failure for operator "%s": d1.yMax( ) = %e != 0.' % ( operator, d1.yMax( ) ) )
    if( m2 is not m3 ) : errs.append( 'failure for operator "%s": m2 is not m3' % operator )

check( '+', value )
check( '-', value )
check( '*', value )
check( '/', value )

other = pointwiseXY_C.pointwiseXY_C( [ [ 0, 1 ], [ 1, 4 ], [ 3, 3 ] ] )
check( '+', other )
check( '-', other )
check( '*', other )
check( '/', other )

if( len( errs ) > 0 ) :
    for err in errs : print '   ', err
    raise Exception( '%s: in place binary operator failed' % __file__ )
