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
sys.path.insert( 0, '../../../../lib' )

import os
import pointwiseXY_C

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__

CPATH = '../../../Test/initial'

os.system( 'cd %s; make -s clean; ./testInterpolation_1' % CPATH )

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
    if( sv1 != sv2 ) : print '<%s> <%s>' % ( sv1, sv2 )
    if( sv1 != sv2 ) : raise Exception( '%s: values %e and %e diff by %e for label = %s' % ( __file__, v1, v2, v2 - v1, label ) )

def getXYData( fileName, accuracy, biSectionMax, interpolation ) :

    f = open( os.path.join( CPATH, fileName  + '.dat' ) )
    ls = f.readlines( )
    ls, length = getIntegerValue( 'length', ls )
    data = [ map( float, ls[i].split( ) ) for i in xrange( length ) ]
    data = pointwiseXY_C.pointwiseXY_C( data, initialSize = len( data ), overflowSize = 10, accuracy = accuracy, biSectionMax = biSectionMax, interpolation = interpolation )
    return( data )

def checkDatas( label, d1, d2 ) :

    if( len( d1 ) != len( d1 ) ) : raise Exception( '%s: for %s len( d1 ) = %d != len( d2 ) = %d' % ( __file__, label, len( d1 ), len( d2 ) ) )
    for i, xy in enumerate( d1 ) :
        compareValues( '%s: x-values at index %d' % ( label, i ), xy[0], d2[i][0] )
        compareValues( '%s: y-values at index %d' % ( label, i ), xy[1], d2[i][1] )

accuracy = 1e-3
biSectionMax = 10
curve_u_loglog = getXYData( 'curve_u_loglog', accuracy, biSectionMax, 'log,log' )
curve_u_interpolatedToLinearC = getXYData( 'curve_u_interpolatedToLinear', accuracy, biSectionMax, 'linear,linear' )
curve_u_loglog_dense = getXYData( 'curve_u_loglog_dense', accuracy, biSectionMax, 'log,log' )

curve_u_interpolatedToLinear = curve_u_loglog.changeInterpolation( 'linear,linear' )
checkDatas( 'curve_u_loglog', curve_u_interpolatedToLinearC, curve_u_interpolatedToLinear )
