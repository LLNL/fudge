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

def compareLists( comparing, list1, list2, CList1, CList2 ) :

    e1, e2 = eval( "list1 %s list2" % comparing ), eval( "CList1 %s CList2" % comparing )
    if( e1 != e2 ) :
        print( comparing, eval( "list1 %s list2" % comparing ), eval( "CList1 %s CList2" % comparing ) )
        print( list1 )
        print( list2 )
        raise Exception( 'Comparison failed for operator "%s"' % comparing )

def check( list1, list2 ) :

    CList1 = listOfDoubles_C.listOfDoubles_C( list1 )
    CList2 = listOfDoubles_C.listOfDoubles_C( list2 )
    compareLists( '<',  list1, list2, CList1, CList2 )
    compareLists( '<=', list1, list2, CList1, CList2 )
    compareLists( '==', list1, list2, CList1, CList2 )
    compareLists( '!=', list1, list2, CList1, CList2 )
    compareLists( '>=', list1, list2, CList1, CList2 )
    compareLists( '>',  list1, list2, CList1, CList2 )

check( [], [] )
values1 = list( range( -10, 200 ) )
for i1 in range( 100 ) :
    check( [ float( random.choice( values1 ) ) for j1 in range( random.choice( values1 ) ) ],
            [ float( random.choice( values1 ) ) for j1 in range( random.choice( values1 ) ) ] )
