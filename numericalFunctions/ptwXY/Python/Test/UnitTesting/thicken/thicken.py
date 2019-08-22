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

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__

CPATH = '../../../../Test/UnitTesting/thicken'

os.system( 'cd %s; thicken -v > v' % CPATH )
f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

def getData( ls, hasLabel ) :

    i = 0
    for l in ls :
        if( l.strip( ) != '' ) : break
        i = i + 1
    ls = ls[i:]
    if( len( ls ) == 0 ) : return( None, None, None )
    label = None
    if( hasLabel ) : label, ls = ls[0].strip( ), ls[1:]
    length, ls = ls[0], ls[1:]
    if( '# length = ' != length[:11] ) : raise Exception( 'Line does not contain length info: "%s"' % ls[0].strip( ) )
    length = int( length.split( '=' )[1] )
    data = [ map( float, ls[i].split( )[:2] ) for i in xrange( length ) ]
    return( ls[length:], label, pointwiseXY_C.pointwiseXY_C( data, initialSize = 10, overflowSize = 10 ) )

def compareValues( label, i, v1, v2 ) :

    sv1, sv2 = '%.7g' % v1, '%.7g' % v2
    if( sv1 != sv2 ) : raise Exception( 'Values %s %s diff at %d for label = %s' % ( v1, v2, i, label ) )

def thicken( label, original, data ) :

    values = label.split( ':' )[1].split( '=' )
    sectionSubdivideMax = int( values[1].split( )[0] )
    dxMax = float( values[2].split( )[0] )
    fxMax = float( values[3].split( )[0] )
    thick = original.thicken( sectionSubdivideMax = sectionSubdivideMax, dxMax = dxMax, fxMax = fxMax )
    if( len( data ) != len( thick ) ) : raise Exception( 'len( data ) = %d != len( thick ) = %d for label = "%s"' % \
        ( len( data ), len( thick ), label ) )
    if( 'log-log' in label ) : return
    for i, xy in enumerate( data ) :
        xc, yc = xy
        xp, yp = thick[i]
        compareValues( label, i, xc, xp )
        compareValues( label, i, yc, yp )

hasLabel = False
while( 1 ) :
    ls, label, data = getData( ls, hasLabel )
    if( ls is None ) : break
    if( hasLabel ) :
        thicken( label, original, data )
    else :
        original = data
    hasLabel = True
