#! /bin/env python

# <<BEGIN-copyright>>
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
