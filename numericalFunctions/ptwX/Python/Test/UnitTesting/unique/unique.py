# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import random
from numericalFunctions import listOfDoubles_C

debugging = False

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print( __file__ )

def compareLists( Msg, list_, uniquePy, uniqueC ) :

    uniqueC2 = [ item for item in uniqueC ]
    if( uniquePy != uniqueC2 ) :
        print( len( list_ ) )
        print( list_ )
        print( uniquePy )
        print( uniqueC.toString( format = "%.0f", valuesPerLine = 25, valueSeparator = ', ' ) )
        raise Exception( 'List differ for %s' % Msg )

def check( list_ ) :

    CList = listOfDoubles_C.listOfDoubles_C( list_ )
    uniquePy = []
    for item in list_ :
        if( item not in uniquePy ) : uniquePy.append( item )

    uniqueC = CList.unique( )
    compareLists( 'default', list_, uniquePy, uniqueC )

    uniqueC = CList.unique( order = -1 )
    compareLists( 'descending', list_, sorted( uniquePy, reverse = True ), uniqueC )

    uniqueC = CList.unique( order = 0 )
    compareLists( 'default 2', list_, uniquePy, uniqueC )

    uniqueC = CList.unique( order = 1 )
    compareLists( 'ascending', list_, sorted( uniquePy ), uniqueC )
    if( debugging ) :
        print( CList.toString( format = "%.0f", valuesPerLine = 25, valueSeparator = ', ' ) )
        print( CList.unique( ).toString( format = "%.0f", valuesPerLine = 25, valueSeparator = ', ' ) )
        print( CList.unique( order = -1 ).toString( format = "%.0f", valuesPerLine = 25, valueSeparator = ', ' ) )
        print( CList.unique( order =  0 ).toString( format = "%.0f", valuesPerLine = 25, valueSeparator = ', ' ) )
        print( CList.unique( order =  1 ).toString( format = "%.0f", valuesPerLine = 25, valueSeparator = ', ' ) )

values = list( range( -3, 5 ) )
check( [ random.choice( values ) for j1 in range( 10 ) ] )
check( [ random.choice( values ) for j1 in range( 100 ) ] )
values = list( range( -30, 500 ) )
check( [ random.choice( values ) for j1 in range( 100000 ) ] )
