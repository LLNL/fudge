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
