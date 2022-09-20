# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

from numericalFunctions import pointwiseXY_C

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print( __file__ )

CPATH = '../../../Test/initial'

os.system( 'cd %s; make -s clean; ./testExp_1' % CPATH )

def skipBlankLines( ls ) :

    i = 0
    for i, l in enumerate( ls ) :
        if( l.strip( ) != '' ) : break
    ls = ls[i:]
    if( ( len( ls ) == 1 ) and ( ls[0].strip( ) == '' ) ) : ls = []
    return( ls )

def getIntegerValue( name, ls ) :

    s = "# %s = " % name
    n = len( s )
    if( ls[0][:n] != s ) : raise Exception( '%s: missing %s info: "%s"' % ( __file__, name, ls[0][:-1] ) )
    value = int( ls[0].split( '=' )[1] )
    return( ls[1:], value )

def getDoubleValue( name, ls ) :

    s = "# %s = " % name
    n = len( s )
    if( ls[0][:n] != s ) : raise Exception( '%s: missing %s info: "%s"' % ( __file__, name, ls[0][:-1] ) )
    value = float( ls[0].split( '=' )[1] )
    return( ls[1:], value )

def compareValues( label, v1, v2 ) :

    sv1, sv2 = '%.12e' % v1, '%.12e' % v2
    sv1, sv2 = '%.7e' % float( sv1 ), '%.7e' % float( sv2 )
    if( sv1 != sv2 ) : print( '<%s> <%s>' % ( sv1, sv2 ) )
    if( sv1 != sv2 ) : raise Exception( '%s: values %e and %e diff by %e for label = %s' % ( __file__, v1, v2, v2 - v1, label ) )

def getXYData( fileName, accuracy, biSectionMax ) :

    f = open( os.path.join( CPATH, fileName  + '.dat' ) )
    ls = f.readlines( )
    ls, length = getIntegerValue( 'length', ls )
    data = [ list( map( float, ls[i].split( ) ) ) for i in range( length ) ]
    data = pointwiseXY_C.pointwiseXY_C( data, initialSize = len( data ), overflowSize = 10, accuracy = accuracy, biSectionMax = biSectionMax )
    return( data )

def checkDatas( label, d1, d2 ) :

    if( len( d1 ) != len( d1 ) ) : raise Exception( '%s: for %s len( d1 ) = %d != len( d2 ) = %d' % ( __file__, label, len( d1 ), len( d2 ) ) )
    for i, xy in enumerate( d1 ) :
        compareValues( '%s: x-values at index %d' % ( label, i ), xy[0], d2[i][0] )
        compareValues( '%s: y-values at index %d' % ( label, i ), xy[1], d2[i][1] )

f = open( os.path.join( CPATH, 'info.dat' ) )
ls = f.readlines( )
f.close( )
ls, accuracy = getDoubleValue( 'accuracy', ls )
ls, biSectionMax = getIntegerValue( 'biSectionMax', ls )

curve_u = getXYData( 'curve_u', accuracy, biSectionMax )
exactExp_u = getXYData( 'exactExp_u', accuracy, biSectionMax )
exp_u = getXYData( 'exp_u', accuracy, biSectionMax )
exp_u_times_exp_minusUC = getXYData( 'exp_u_times_exp_minusU', accuracy, biSectionMax )
exp = curve_u.exp( 1. )

checkDatas( 'exp', exp_u, exp )

exp_u_times_exp_minusU = curve_u.exp( 1 ) * curve_u.exp( -1 )
checkDatas( 'exp(y) * exp(-y)', exp_u_times_exp_minusU, exp_u_times_exp_minusUC ) 
