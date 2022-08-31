# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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

def compareLists( Msg, list1, list2, CList ) :

    CList2 = [ item for item in CList ]
    Pylist = list1 + list2
    if( Pylist != CList2 ) :
        print( len( list1 ), len( list2 ) )
        print( list1 )
        print( list2 )
        print( Pylist )
        print( CList.toString( format = "%.0f", valuesPerLine = 25, valueSeparator = ', ' ) )
        print( CList2 )
        raise Exception( 'List differ for %s' % Msg )

def check( list1, list2 ) :

    CList1 = listOfDoubles_C.listOfDoubles_C( list1 )
    CList2 = listOfDoubles_C.listOfDoubles_C( list2 )
    compareLists( 'CList1 + list1',  list1, list1, CList1 + list1 )
    compareLists( 'CList1 + CList1', list1, list1, CList1 + CList1 )
    compareLists( 'CList1 + list2',  list1, list2, CList1 + list2 )
    compareLists( 'CList1 + CList2', list1, list2, CList1 + CList2 )
    compareLists( 'CList2 + list1',  list2, list1, CList2 + list1 )
    compareLists( 'CList2 + CList1', list2, list1, CList2 + CList1 )
    compareLists( 'CList2 + list2',  list2, list2, CList2 + list2 )
    compareLists( 'CList2 + CList2', list2, list2, CList2 + CList2 )
    CList1 += list1
    compareLists( 'CList1 += list1',  list1, list1, CList1 )

values1 = list( range( -3, 5 ) )
values2 = list( range( -30, 5000 ) )
for i1 in range( 30 ) :
    check( [ float( random.choice( values1 ) ) for j1 in range( random.choice( values1 ) ) ], [ float( random.choice( values1 ) ) for j1 in range( random.choice( values1 ) ) ] )
    check( [ float( random.choice( values2 ) ) for j1 in range( random.choice( values2 ) ) ], [ float( random.choice( values2 ) ) for j1 in range( random.choice( values2 ) ) ] )
