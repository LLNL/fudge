# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import random

from numericalFunctions import listOfDoubles_C

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print( __file__ )

def compareLists( Msg, mul, list1, CList ) :

    CList2 = [ item for item in CList ]
    Pylist = mul * list1
    if( Pylist != CList2 ) :
        print( mul, len( list1 ), len( Pylist ), len( CList ) )
        print( list1 )
        print( Pylist )
        print( CList2 )
        raise Exception( 'List differ for %s' % Msg )

def check( mul, list1 ) :

    CList1 = listOfDoubles_C.listOfDoubles_C( list1 )
    compareLists( 'mul * CList1', mul, list1, mul * CList1 )
    compareLists( 'mul * CList1', mul, list1, CList1 * mul )
    CList1 *= mul
    compareLists( 'CList1 *= mul', mul, list1, CList1 )

values1 = list( range( -3, 5 ) )
values2 = list( range( -30, 500 ) )
for i1 in range( 30 ) :
    check( random.choice( values1 ), [ float( random.choice( values1 ) ) for j1 in range( random.choice( values1 ) ) ] )
    check( random.choice( values2 ), [ float( random.choice( values2 ) ) for j1 in range( random.choice( values2 ) ) ] )
