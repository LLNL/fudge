#! /bin/env python

# <<BEGIN-copyright>>
# <<END-copyright>>

import sys
sys.path.insert( 0, '../../../../../lib' )

import os
import pointwiseXY_C

if( 'CHECKOPTIONS' in os.environ ) :
    options = os.environ['CHECKOPTIONS'].split( )
    if( '-e' in options ) : print __file__

CPATH = '../../../../Test/UnitTesting/binaryMath'

os.system( 'cd %s; make -s clean; ./add -v > v' % CPATH )

def skipBlankLines( ls ) :

    i = 0
    for i, l in enumerate( ls ) :
        if( l.strip( ) != '' ) : break
    ls = ls[i:]
    if( ( len( ls ) == 1 ) and ( ls[0].strip( ) == '' ) ) : ls = []
    return( ls )

def getIntegerValue( name, ls ) :

    s = "# %s = " % name
    n = len( s )
    if( ls[0][:n] != s ) : raise Exception( '%s: missing %s info: "%s"' % ( __file__, name, ls[0][:-1] ) )
    value = int( ls[0].split( '=' )[1] )
    return( ls[1:], value )

def getDoubleValue( name, ls ) :

    s = "# %s = " % name
    n = len( s )
    if( ls[0][:n] != s ) : raise Exception( '%s: missing %s info: "%s"' % ( __file__, name, ls[0][:-1] ) )
    value = float( ls[0].split( '=' )[1] )
    return( ls[1:], value )

def compareValues( label, i, v1, v2 ) :

    sv1, sv2 = '%.12e' % v1, '%.12e' % v2
    sv1, sv2 = '%.7e' % float( sv1 ), '%.7e' % float( sv2 )
    if( sv1 != sv2 ) : print '<%s> <%s>' % ( sv1, sv2 )
    if( sv1 != sv2 ) : raise Exception( '%s: values %e and %e diff by %e at %d for label = %s' % ( __file__, v1, v2, v2 - v1, i, label ) )

def getXYData( ls, biSectionMax, accuracy ) :

    ls, length = getIntegerValue( 'length', ls )
    data = [ map( float, ls[i].split( ) ) for i in xrange( length ) ]
    data = pointwiseXY_C.pointwiseXY_C( data, initialSize = len( data ), overflowSize = 10, biSectionMax = biSectionMax, accuracy = accuracy, safeDivide = True )
    ls = ls[length:]
    ls = skipBlankLines( ls )
    return( ls, data )

def checkDivision( count, ls ) :

    ls, biSectionMax = getDoubleValue( 'biSectionMax', ls )
    ls, accuracy = getDoubleValue( 'accuracy', ls )

    ls, addend1 = getXYData( ls, biSectionMax, accuracy )
    ls, addend2 = getXYData( ls, biSectionMax, accuracy )
    ls, add1C = getXYData( ls, biSectionMax, accuracy )
    ls, add2C = getXYData( ls, biSectionMax, accuracy )
    ls, diffC = getXYData( ls, biSectionMax, accuracy )
    add1 = addend1 + addend2
    add2 = addend2 + addend1
    diff = add1 - add1
    if( len( diffC ) != len( diff ) ) : raise Exception( '%s: at %d len( diffC ) = %d != len( diffC ) = %d' % ( __file__, count, len( diffC ), len( diffC ) ) )
    for i, xy in enumerate( diff ) :
        if( xy[1] != 0 ) : raise Exception( "diff at index = %s is %s and not 0." % ( i, xy[1] ) )
    for i, xy in enumerate( diffC ) :
        compareValues( "x diff", count, xy[0], diff[i][0] )
        compareValues( "y diff", count, xy[1], diff[i][1] )
    return( ls )

f = open( os.path.join( CPATH, 'v' ) )
ls = f.readlines( )
f.close( )

count = 0
while( len( ls ) ) :
    count += 1
    ls = checkDivision( count, ls )
