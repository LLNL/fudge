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

CPATH = '../../../../Test/UnitTesting/convolution'

os.system( 'cd %s; ./gaussian -v > v' % CPATH )
f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

def getDoubleValue( name, ls ) :

    s = "# %s = " % name
    n = len( s )
    if( ls[0][:n] != s ) : raise Exception( '%s: line does not contain %s info: "%s"' % ( __file__, name, ls[0][:-1] ) )
    value = float( ls[0].split( '=' )[1] )
    return( ls[1:], value )

def compareValues( label, i, v1, v2 ) :

    sv1, sv2 = '%.12g' % v1, '%.12g' % v2
    sv1, sv2 = '%.8g' % float( sv1 ), '%.8g' % float( sv2 )
    if( sv1 != sv2 ) : print '<%s> <%s>' % ( sv1, sv2 )
    if( sv1 != sv2 ) : raise Exception( '%s: values %s %s diff by %g at %d for label = %s' % ( __file__, v1, v2, v2 - v1, i, label ) )

def getData( ls ) :

    i = 0
    for l in ls :
        if( l.strip( ) != '' ) : break
        i = i + 1
    ls = ls[i:]
    if( len( ls ) == 0 ) : return( ls )
    if( ls[0][:9] == "# area = " ) : ls = ls[1:]
    if( len( ls ) == 0 ) : return( ls )

    label = ls[0][2:-1]
    if( label not in [ 'ptwXY_createGaussianCenteredSigma1', 'ptwXY_createGaussian' ] ) :
        raise Exception( '%s: invalid label = "%s"' % ( __file__, ls[0][:-1] ) )
    ls = ls[1:]
    ls, accuracy = getDoubleValue( 'accuracy', ls )
    ls, offset = getDoubleValue( 'offset', ls )
    ls, sigma = getDoubleValue( 'sigma', ls )
    ls, amplitude = getDoubleValue( 'amplitude', ls )
    ls, xMin = getDoubleValue( 'xMin', ls )
    ls, xMax = getDoubleValue( 'xMax', ls )
    ls, dullEps = getDoubleValue( 'dullEps', ls )

    length, ls = ls[0], ls[1:]
    if( '# length = ' != length[:11] ) : raise Exception( '%s: line does not contain length info: "%s"' % ( __file__, length.strip( ) ) )
    length = int( length.split( '=' )[1] )

    data = [ map( float, ls[i].split( )[:2] ) for i in xrange( length ) ]
    data = pointwiseXY_C.pointwiseXY_C( data, initialSize = len( data ), overflowSize = 10 )
    ls = ls[length:]
    if( label == 'ptwXY_createGaussianCenteredSigma1' ) :
        gaussian = pointwiseXY_C.basicGaussian( accuracy )
    else :
        gaussian = pointwiseXY_C.gaussian( accuracy = accuracy, xMin = xMin, xMax = xMax, offset = offset, sigma = sigma, 
            amplitude = amplitude, dullEps = dullEps )
    if( len( data ) != len( gaussian ) ) : raise Exception( '%s: len( data ) = %d != len( gaussian ) = %d for label "%s"' %
        ( __file__, len( data ), len( gaussian ), label ) )
    for i , dXY in enumerate( data ) :
        gXY = gaussian[i]
        compareValues( label, i, dXY[0], gXY[0] )
        compareValues( label, i, dXY[1], gXY[1] )
    return( ls )

while( len( ls ) ) : ls = getData( ls )
