# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

from numericalFunctions import pointwiseXY_C

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print( __file__ )

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
    if( sv1 != sv2 ) : print( '<%s> <%s>' % ( sv1, sv2 ) )
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

    data = [ list( map( float, ls[i].split( )[:2] ) ) for i in range( length ) ]
    data = pointwiseXY_C.pointwiseXY_C( data, initialSize = len( data ), overflowSize = 10 )
    ls = ls[length:]
    if( label == 'ptwXY_createGaussianCenteredSigma1' ) :
        gaussian = pointwiseXY_C.basicGaussian( accuracy )
    else :
        gaussian = pointwiseXY_C.gaussian( accuracy = accuracy, domainMin = xMin, domainMax = xMax, offset = offset, sigma = sigma, 
            amplitude = amplitude, dullEps = dullEps )
    if( len( data ) != len( gaussian ) ) : raise Exception( '%s: len( data ) = %d != len( gaussian ) = %d for label "%s"' %
        ( __file__, len( data ), len( gaussian ), label ) )
    for i , dXY in enumerate( data ) :
        gXY = gaussian[i]
        compareValues( label, i, dXY[0], gXY[0] )
        compareValues( label, i, dXY[1], gXY[1] )
    return( ls )

while( len( ls ) ) : ls = getData( ls )
