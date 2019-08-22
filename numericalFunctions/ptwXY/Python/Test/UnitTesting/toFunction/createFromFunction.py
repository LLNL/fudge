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
sys.path.insert( 0, '../../Utilities' )

import os, math
import pointwiseXY_C
import utilities

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__

CPATH = '../../../../Test/UnitTesting/toFunction'

os.system( 'cd %s; createFromFunction -v > v' % CPATH )
f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

def getData( ls ) :

    ls = utilities.skipBlankLines( ls )
    if( len( ls ) == 0 ) : return( None, None, None, None, None, None )
    label, accuracy = None, None
    label, ls = ls[0].strip( ), ls[1:]
    if( label != '# Xs' ) : raise Exception( '%s: missing line "# Xs": found "%s"' % ( __file__, label ) )
    ls, Xs = utilities.getXData( ls, '=' )
    ls, accuracy = utilities.getDoubleValue( 'accuracy', ls )
    ls, biSectionMax = utilities.getDoubleValue( 'biSectionMax', ls )
    ls, data = utilities.getXYData( ls )
    return( ls, label, accuracy, biSectionMax, Xs, pointwiseXY_C.pointwiseXY_C( data, accuracy = accuracy ) )

def compareValues( label, i, v1, v2 ) :

    sv1, sv2 = '%.7g' % v1, '%.7g' % v2
    if( sv1 != sv2 ) : raise Exception( '%s: values %s %s diff at %d for label = %s' % ( __file__, v1, v2, i, label ) )

def f_x_Sin_xx( x, args ) :

       return( x * math.sin( x * x ) )

def createFromFunction( label, accuracy, biSectionMax, Xs, original ) :

    data = pointwiseXY_C.createFromFunction( Xs, f_x_Sin_xx, None, accuracy, biSectionMax, checkForRoots = 1 )
    if( len( data ) != len( original ) ) : raise Exception( '%s: len( data ) = %d != len( original ) = %d for label = "%s"' % \
        ( __file__, len( data ), len( original ), label ) )

while( 1 ) :
    ls, label, accuracy, biSectionMax, Xs, original = getData( ls )
    if( ls is None ) : break
    createFromFunction( label, accuracy, biSectionMax, Xs, original )
