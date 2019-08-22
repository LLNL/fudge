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

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__

def compareLists( list1, CList ) :

    CList2 = [ item for item in CList ]
    Pylist = [ -x  for x in list1 ]
    if( Pylist != CList2 ) :
        print list1
        print Pylist
        print CList2
        raise Exception( 'List differ' )

def check( list1 ) :

    CList1 = listOfDoubles_C.listOfDoubles_C( list1 )
    CList2 = -CList1
    compareLists( list1, CList2 )

values1 = xrange( -3, 5 )
values2 = xrange( -300, 5000 )
for i1 in xrange( 100 ) :
    check( [ float( random.choice( values1 ) ) for j1 in xrange( random.choice( values1 ) ) ] )
    check( [ float( random.choice( values2 ) ) for j1 in xrange( random.choice( values2 ) ) ] )
