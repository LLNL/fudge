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

import os
import pointwiseXY_C

verbose = False

options = []
if( 'CHECKOPTIONS' in os.environ ) : options = os.environ['CHECKOPTIONS'].split( )
for argv in sys.argv[1:] : options += [ argv ]

if( '-e' in options ) : print __file__
if( '-v' in options ) : verbose = True

CPATH = '../../../../Test/UnitTesting/Others'

os.system( 'cd %s; make -s clean; ./mutualify -v > v' % CPATH )

def skipBlankLines( ls ) :

    i = 0
    for i, l in enumerate( ls ) :
        if( l.strip( ) != '' ) : break
    ls = ls[i:]
    if( ( len( ls ) == 1 ) and ( ls[0].strip( ) == '' ) ) : ls = []
    return( ls )

def getIntegerValue( name, ls, verbose ) :

    s = "# %s = " % name
    n = len( s )
    if( ls[0][:n] != s ) : raise Exception( '%s: missing %s info: "%s"' % ( __file__, name, ls[0][:-1] ) )
    value = int( ls[0].split( '=' )[1] )
    if( verbose ) : print ls[0][:-1]
    return( ls[1:], value )

def getDoubleValue( name, ls ) :

    s = "# %s = " % name
    n = len( s )
    if( ls[0][:n] != s ) : raise Exception( '%s: missing %s info: "%s"' % ( __file__, name, ls[0][:-1] ) )
    value = float( ls[0].split( '=' )[1] )
    if( verbose ) : print ls[0][:-1]
    return( ls[1:], value )

def compareValues( label, i, v1, v2 ) :

    sv1, sv2 = '%.12e' % v1, '%.12e' % v2
    sv1, sv2 = '%.7e' % float( sv1 ), '%.7e' % float( sv2 )
    if( sv1 != sv2 ) : print '<%s> <%s>' % ( sv1, sv2 )
    if( sv1 != sv2 ) : raise Exception( '%s: values %e and %e diff by %e at %d for label = %s' % ( __file__, v1, v2, v2 - v1, i, label ) )

def getXYData( ls, verbose ) :

    ls, length = getIntegerValue( 'length', ls, False )
    data = [ map( float, ls[i].split( ) ) for i in xrange( length ) ]
    data = pointwiseXY_C.pointwiseXY_C( data, initialSize = len( data ), overflowSize = 10 )
    ls = ls[length:]
    ls = skipBlankLines( ls )
    if( verbose ) : printDataIfVerbose( data )
    return( ls, data )

def printDataIfVerbose( data ) :

    if( verbose ) :
        print "# length = %d" % len( data )
        print data
        print

def checkMutualifing( count, ls ) :

    ls = skipBlankLines( ls )
    ls, i1 = getIntegerValue( 'i1', ls, verbose )
    ls, i2 = getIntegerValue( 'i2', ls, verbose )
    ls, lowerEps = getDoubleValue( 'lowerEps', ls )
    ls, upperEps = getDoubleValue( 'upperEps', ls )
    ls, positiveXOnly = getIntegerValue( 'positiveXOnly', ls, verbose )
    ls, d1 = getXYData( ls, verbose )
    ls, d2 = getXYData( ls, verbose )
    ls, mutual1C = getXYData( ls, False )
    ls, mutual2C = getXYData( ls, False )
    d1m, d2m = d1.mutualify( lowerEps, upperEps, positiveXOnly, d2, lowerEps, upperEps, positiveXOnly )
    printDataIfVerbose( d1m )
    printDataIfVerbose( d2m )

    if( len( mutual1C ) != len( d1m ) ) : raise Exception( '%s: at %d len( mutual1C ) = %d != len( d1m ) = %d' % \
        ( __file__, count, len( mutual1C ), len( d1m ) ) )
    for i, xy in enumerate( d1m ) :
        compareValues( "d1m x mutualified", count, xy[0], mutual1C[i][0] )
        compareValues( "d1m y mutualified", count, xy[1], mutual1C[i][1] )

    if( len( mutual2C ) != len( d2m ) ) : raise Exception( '%s: at %d len( mutual2C ) = %d != len( d2m ) = %d' % \
        ( __file__, count, len( mutual2C ), len( d2m ) ) )
    for i, xy in enumerate( d2m ) :
        compareValues( "d2m x mutualified", count, xy[0], mutual2C[i][0] )
        compareValues( "d2m y mutualified", count, xy[1], mutual2C[i][1] )
    return( ls )

if( verbose ) :
    doc = pointwiseXY_C.pointwiseXY_C.mutualify.__doc__.split( '\n' )
    for d in doc : print "#", d

f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

count = 0
while( len( ls ) ) :
    count += 1
    ls = checkMutualifing( count, ls )
