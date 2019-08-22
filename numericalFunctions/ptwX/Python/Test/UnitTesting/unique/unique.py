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

import listOfDoubles_C
import os, random

debugging = False

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__

def compareLists( Msg, list_, uniquePy, uniqueC ) :

    uniqueC2 = [ item for item in uniqueC ]
    if( uniquePy != uniqueC2 ) :
        print len( list_ )
        print list_
        print uniquePy
        print uniqueC.toString( format = "%.0f", valuesPerLine = 25, valueSeparator = ', ' )
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
        print CList.toString( format = "%.0f", valuesPerLine = 25, valueSeparator = ', ' )
        print CList.unique( ).toString( format = "%.0f", valuesPerLine = 25, valueSeparator = ', ' )
        print CList.unique( order = -1 ).toString( format = "%.0f", valuesPerLine = 25, valueSeparator = ', ' )
        print CList.unique( order =  0 ).toString( format = "%.0f", valuesPerLine = 25, valueSeparator = ', ' )
        print CList.unique( order =  1 ).toString( format = "%.0f", valuesPerLine = 25, valueSeparator = ', ' )

values = xrange( -3, 5 )
check( [ random.choice( values ) for j1 in xrange( 10 ) ] )
check( [ random.choice( values ) for j1 in xrange( 100 ) ] )
values = xrange( -30, 500 )
check( [ random.choice( values ) for j1 in xrange( 100000 ) ] )
