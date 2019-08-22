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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../../../../lib' )

import listOfDoubles_C
import os, random

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__

def compareLists( CList1, sep ) :

    for prefixWidth in range( 5 ) :
        for suffixWidth in range( 5 ) : compareLists2( CList1, sep, prefixWidth, suffixWidth )

def compareLists2( CList1, sep, prefixWidth, suffixWidth ) :

    str1 = toString( CList1, sep, prefixWidth, suffixWidth )
    CList2, extraCharacters = listOfDoubles_C.createFromString( str1, sep = sep )
    if( CList1 == CList2 ) : return
    print CList1 == CList2
    print 'sep = <%s>' % sep, prefixWidth, suffixWidth
    print CList1
    print CList2
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
