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

import os
import pointwiseXY_C

verbose = False

options = []
if( 'CHECKOPTIONS' in os.environ ) : options = os.environ['CHECKOPTIONS'].split( )
for argv in sys.argv[1:] : options += [ argv ]

if( '-e' in options ) : print __file__
if( '-v' in options ) : verbose = True

CPATH = '../../../../Test/UnitTesting/Others'

os.system( 'cd %s; make -s clean; ./areDomainsMutual -v > v' % CPATH )

def skipBlankLines( ls ) :

    i = 0
    for i, l in enumerate( ls ) :
        if( l.strip( ) != '' ) : break
    ls = ls[i:]
    if( ( len( ls ) == 1 ) and ( ls[0].strip( ) == '' ) ) : ls = []
    return( ls )

def getIntegerValue( name, ls, verbose ) :

    ls = skipBlankLines( ls )
    s = "# %s = " % name
    n = len( s )
    if( ls[0][:n] != s ) : raise Exception( '%s: missing %s info: "%s"' % ( __file__, name, ls[0][:-1] ) )
    value = int( ls[0].split( '=' )[1] )
    if( verbose ) : print ls[0][:-1]
    return( ls[1:], value )

def getDoubleValue( name, ls ) :

    ls = skipBlankLines( ls )
    s = "# %s = " % name
    n = len( s )
    if( ls[0][:n] != s ) : raise Exception( '%s: missing %s info: "%s"' % ( __file__, name, ls[0][:-1] ) )
    value = float( ls[0].split( '=' )[1] )
    if( verbose ) : print ls[0][:-1]
    return( ls[1:], value )

def getXYData( ls, verbose ) :

    ls = skipBlankLines( ls )
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

def checkAreMutual( count, ls ) :

    ls, areMutual = getIntegerValue( 'areMutual', ls, verbose )
    ls, d1 = getXYData( ls, verbose )
    ls, d2 = getXYData( ls, verbose )
    pyResults = d1.areDomainsMutual( d2 )
    if( ( pyResults and not( areMutual ) ) or ( not( pyResults ) and areMutual ) ) :
        print '%s: for count = %d, pyResults = %s and areMutual = %s' % ( __file__, count, pyResults, areMutual )

    return( ls )

if( verbose ) :
    doc = pointwiseXY_C.pointwiseXY_C.areDomainsMutual.__doc__.split( '\n' )
    for d in doc : print "#", d

f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

count = 0
while( len( ls ) ) :
    count += 1
    ls = checkAreMutual( count, ls )
