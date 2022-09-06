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

def compareLists( list1, CList ) :

    CList2 = [ item for item in CList ]
    Pylist = [ -x  for x in list1 ]
    if( Pylist != CList2 ) :
        print( list1 )
        print( Pylist )
        print( CList2 )
        raise Exception( 'List differ' )

def check( list1 ) :

    CList1 = listOfDoubles_C.listOfDoubles_C( list1 )
    CList2 = -CList1
    compareLists( list1, CList2 )

values1 = list( range( -3, 5 ) )
values2 = list( range( -300, 5000 ) )
for i1 in range( 100 ) :
    check( [ float( random.choice( values1 ) ) for j1 in range( random.choice( values1 ) ) ] )
    check( [ float( random.choice( values2 ) ) for j1 in range( random.choice( values2 ) ) ] )
