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

def compareLists( CList1, sep ) :

    for prefixWidth in range( 5 ) :
        for suffixWidth in range( 5 ) : compareLists2( CList1, sep, prefixWidth, suffixWidth )

def compareLists2( CList1, sep, prefixWidth, suffixWidth ) :

    str1 = toString( CList1, sep, prefixWidth, suffixWidth )
    CList2, extraCharacters = listOfDoubles_C.createFromString( str1, sep = sep )
    if( CList1 == CList2 ) : return
    print( CList1 == CList2 )
    print( 'sep = <%s>' % sep, prefixWidth, suffixWidth )
    print( CList1 )
    print( CList2 )
    raise Exception( 'Comparison failed' )

def toString( values, sep, prefixWidth, suffixWidth ) :

    strList = []
    for value in values : 
        prefix = random.choice( range( prefixWidth + 1 ) ) * ' '
        suffix = random.choice( range( suffixWidth + 1 ) ) * ' '
        strList.append( "%s%.8g%s" % ( prefix, value, suffix ) )
    return( sep.join( strList ) )

def check( values ) :

    CList1 = listOfDoubles_C.listOfDoubles_C( values )
    compareLists( CList1, ' ' )
    compareLists( CList1, ',' )
    compareLists( CList1, ':' )

vMin, vMax = -1e3, 3.14e6
for i1 in range( 44 ) :
    values = [ float( "%.6g" % random.uniform( vMin, vMax ) ) for j1 in range( i1 ) ]
    check( values )
