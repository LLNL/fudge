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

CPATH = '../../../../Test/UnitTesting/integrate'

os.system( 'cd %s; make -s clean; ./groupTestAll1 -v' % CPATH )

def getIntegerValue( name, ls ) :

    s = "# %s = " % name
    n = len( s )
    if( ls[0][:n] != s ) : raise Exception( '%s: %s does not contain %s info: "%s"' % ( __file__, name, ls[0][:-1] ) )
    value = int( ls[0].split( '=' )[1] )
    return( ls[1:], value )

def compareValues( label, i, v1, v2 ) :

    sv1, sv2 = '%.12e' % v1, '%.12e' % v2
    sv1, sv2 = '%.7e' % float( sv1 ), '%.7e' % float( sv2 )
    if( sv1 != sv2 ) : print '<%s> <%s>' % ( sv1, sv2 )
    if( sv1 != sv2 ) : raise Exception( '%s: values %e and %e diff by %e at %d for label = %s' % ( __file__, v1, v2, v2 - v1, i, label ) )

def compareGroups( fileName, norm, g1 ) :

    label = fileName + '_' + norm
    g2 = getXData( label )
    if( len( g1 ) != len( g2 ) ) : raise Exception( '%s: for %s len( g1 ) = %d != len( g2 ) = %d' %( __file__, label, len( g1 ), len( g2 ) ) )
    for i , g1X in enumerate( g1 ) : compareValues( label, i, g1X, g2[i] )

def getXData( fileName ) :

    fileName_ = os.path.join( CPATH, fileName + '.dat' )
    f = open( fileName_ )
    ls = f.readlines( )
    f.close( )
    ls, length = getIntegerValue( 'length', ls )
    if( len( ls ) != length ) : raise Exception( '%s: len( ls ) = %s != length = %d for file %s' % ( len( ls ), length, fileName ) )
    data = [ float( l ) for l in ls ]
    return( data )
    
def getXYData( fileName ) :

    fileName_ = os.path.join( CPATH, fileName )
    f = open( fileName_ )
    ls = f.readlines( )
    f.close( )
    ls, length = getIntegerValue( 'length', ls )
    data = [ map( float, l.split( ) ) for l in ls ]
    data = pointwiseXY_C.pointwiseXY_C( data, initialSize = len( data ), overflowSize = 10 )
    return( data )

def checkOneFunctionGrouping( fileName, groupBoundaries ) :

    flux = getXYData( fileName + '.dat' )
    flux_None = flux.groupOneFunction( groupBoundaries )
    compareGroups( fileName, 'None', flux_None )
    flux_dx = flux.groupOneFunction( groupBoundaries, norm = 'dx' )
    compareGroups( fileName, 'dx', flux_dx )
    flux_norm = flux.groupOneFunction( groupBoundaries, norm = flux_None )
    compareGroups( fileName, 'norm', flux_norm )
    return( flux, flux_None )

def checkTwoFunctionGrouping( fileName, groupBoundaries, flux, flux_None ) :

    crossSection = getXYData( fileName + '.dat' )
    crossSection_None = crossSection.groupTwoFunctions( flux, groupBoundaries )
    compareGroups( fileName, 'None', crossSection_None )
    crossSection_dx = crossSection.groupTwoFunctions( flux, groupBoundaries, norm = 'dx' )
    compareGroups( fileName, 'dx', crossSection_dx )
    crossSection_norm = crossSection.groupTwoFunctions( flux, groupBoundaries, norm = flux_None )
    compareGroups( fileName, 'norm', crossSection_norm )
    return( crossSection )

def checkThreeFunctionGrouping( fileName, groupBoundaries, flux, crossSection, flux_None ) :

    multiplicity = getXYData( fileName + '.dat' )
    multiplicity_None = multiplicity.groupThreeFunctions( flux, crossSection, groupBoundaries )
    compareGroups( fileName, 'None', multiplicity_None )
    multiplicity_dx = multiplicity.groupThreeFunctions( flux, crossSection, groupBoundaries, norm = 'dx' )
    compareGroups( fileName, 'dx', multiplicity_dx )
    multiplicity_norm = multiplicity.groupThreeFunctions( flux, crossSection, groupBoundaries, norm = flux_None )
    compareGroups( fileName, 'norm', multiplicity_norm )

groupBoundaries = getXData( 'groupBoundaries' )
flux, flux_None = checkOneFunctionGrouping( 'flux', groupBoundaries )
crossSection = checkTwoFunctionGrouping( 'crossSection', groupBoundaries, flux, flux_None )
checkThreeFunctionGrouping( 'multiplicity', groupBoundaries, flux, crossSection, flux_None )
