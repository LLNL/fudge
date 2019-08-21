# <<BEGIN-copyright>>
# <<END-copyright>>

"""
This module contains useful fudge math routines that do not fit into any other module.

cmattoon 3/2011: math functions from fudgemisc moved here
"""

from pqu import physicalQuantityWithUncertainty
from fudge.core.utilities import brb
try :
    import numpy
    numpyFloat64 = numpy.float64( 1. )
except :
    numpyFloat64 = 1.

__metaclass__ = type

# P_0(x)  =  1
# P_1(x)  =  x
# P_2(x)  =    ( 3 x^2 - 1 ) / 2
# P_3(x)  =  x ( 5 x^2 - 3 ) / 2
# P_4(x)  =    ( 35 x^4 - 30 x^2 + 3 ) / 8
# P_5(x)  =  x ( 35 x^4 - 30 x^2 + 3 ) / 8
# P_6(x)  =    ( 231 x^6 - 315 x^4 + 105 x^2 - 5 ) / 16
# P_7(x)  =  x ( 429 x^6 - 693 x^4 + 315 x^2 - 35 ) / 16
# P_8(x)  =    ( 6435 x^8 - 12012 x^6 + 6930 x^4 - 1260 x^2 + 35 ) / 128
# P_9(x)  =  x ( 12155 x^8 - 25740 x^6 + 18018 x^4 - 4620 x^2 + 315 ) / 128
# P_10(x) =    ( 46189 x^10 - 109395 x^8 + 90090 x^6 - 30030 x^4 + 3465 x^2 - 63 ) / 256
# 
# P_{n+1} = ( ( 2 n + 1 ) x P_n(x) - n P_{n-1}(x) ) / ( n + 1 )
#

def Legendre( n, x, checkXRange = True ) :
    """Returns the value of the Legendre function of order n at x.  For n <= 10, use analytical form, 
    for n > 10 use the recursive relationship. (This would be way faster in C or using numpy version)"""

    if( n < 0 ) : raise ValueError( "\nError in endlmathmisc.Legendre: n = %d < 0" % n )
    if checkXRange and ( abs( x ) > 1 ) : raise ValueError( "Legendre: |x| > 1; x = %g" % x )
    if ( n ==  0 ) :   return 1.
    elif ( n ==  1 ) : return x
    elif ( n <= 10 ) :
        x2 = x * x
        if   ( n ==  2 ) : return                1.5 * x2 - 0.5
        elif ( n ==  3 ) : return       (        2.5 * x2 -          1.5 ) * x
        elif ( n ==  4 ) : return       (      4.375 * x2 -         3.75 ) * x2 +       0.375
        elif ( n ==  5 ) : return     ( (      7.875 * x2 -         8.75 ) * x2 +       1.875 ) * x
        elif ( n ==  6 ) : return     ( (    14.4375 * x2 -      19.6875 ) * x2 +      6.5625 ) * x2 -      0.3125
        elif ( n ==  7 ) : return   ( ( (    26.8125 * x2 -      43.3125 ) * x2 +     19.6875 ) * x2 -      2.1875 ) * x
        elif ( n ==  8 ) : return   ( ( (   50.2734375 * x2 -     93.84375 ) * x2 +   54.140625 ) * x2 -     9.84375 ) * x2 + 0.2734375 
        elif ( n ==  9 ) : return ( ( ( (   94.9609375 * x2 -    201.09375 ) * x2 +  140.765625 ) * x2 -    36.09375 ) * x2 + 2.4609375 ) * x
        elif ( n == 10 ) : return ( ( ( ( 180.42578125 * x2 - 427.32421875 ) * x2 + 351.9140625 ) * x2 - 117.3046875 ) * x2 + 13.53515625 ) * x2 - 0.24609375
        else: raise Exception( "Better not get here!" )
    Pn = 0.
    Pnp1 = 1.
    n_ = 0
    twoNp1 = 1
    while( n_ < n ) :
        Pnm1 = Pn
        Pn = Pnp1
        n_p1 = n_ + 1
        Pnp1 = ( twoNp1 * x * Pn - n_ * Pnm1 ) / n_p1
        twoNp1 += 2
        n_ = n_p1
    return( Pnp1 )

def LegendreSeries( coefficients ):
    """ return a function to evaluate the Legendre series with given coefficients:
    >>>foo = LegendreSeries([1.0,0,0,0.85])
    >>>foo(0.5) == 1.0*Legendre(0,0.5) + 0.85*Legendre(3,0.5)
    [True]
    Beware: coefficients as stored in ENDF must be multiplied by factor of (2*L+1)/2
    """
    def series(x):
        res = 0 * x
        # define Legendre terms recursively:
        Pn = 0.
        Pnp1 = 1.
        twoNp1 = 1
        for i in range(len(coefficients)):
            Pnm1 = Pn
            Pn = Pnp1
            res += Pn * coefficients[i]
            Pnp1 = (twoNp1 * x * Pn - i * Pnm1) / (i+1)
            twoNp1 += 2
        return res
    return series

def runningZSum( data, xLabel = None, yLabel = None, zLabel = None, normalize = False ) :
    """Returns the running sum of dy * z (normalized to 1 of normalize is True) for each x as an endl3dmath object.
    Data must be list of ( x, list of ( y, z ) )."""

    d3 = []
    for x_yz in data : d3.append( [ x_yz[0], runningYSum( x_yz[1], normalize = normalize ).data ] )
    import endl3dmathClasses
    return endl3dmathClasses.endl3dmath( d3, xLabel = xLabel, yLabel = yLabel, zLabel = zLabel, checkDataType = 0 )

def runningYSum( data, normalize = False ) :
    """Returns the running sum of dx * y (normalized to 1 of normalize is True) as an endl2dmath object. 
    Data must be list of ( x, y )."""

    x1 = None
    runningSum = []
    for xy in data :
        x2 = xy[0]
        y2 = xy[1]
        if ( x1 == None ) :
            Sum = 0.
        else :
            Sum += 0.5 * ( y2 + y1 ) * ( x2 - x1 )
        runningSum.append( [ x2, Sum ] )
        x1 = x2
        y1 = y2
    if( normalize and ( Sum != 0. ) ) :
        for xy in runningSum : xy[1] /= Sum
    import endl2dmathClasses
    return endl2dmathClasses.endl2dmath( runningSum, checkDataType = 0 )

def ZSum( data ) :
    """Returns the area under the curve z(y) for each x as an endl2dmath object. Data must be list of
    ( x, list of ( y, z ) )."""

    d2 = []
    for x_yz in data : d2.append( [ x_yz[0], YSum( x_yz[1] ) ] )
    import endl2dmathClasses
    return endl2dmathClasses.endl2dmath( d2, checkDataType = 0 )

def YSum( data ) :
    "Returns the area under the curve y(x). Data must be list of list( x, y )."

    x1 = None
    for x2, y2 in data :
        if ( x1 == None ) :
            Sum = 0.
        else :
            Sum += ( y2 + y1 ) * ( x2 - x1 )
        x1 = x2
        y1 = y2
    return 0.5 * Sum

def thickenXYList( list, tester, biSectionMax = 6 ) :
    """This functions takes a list of (x,y) points and a function, tester.evaluateAtX, and bi-sectionally adds points to 
    obtain linear-linear tolerance of the returned list and tester.evaluateAtX to tester.relativeTolerance. At most 
    biSectionMax bi-sections are performed between each consecutive pair of inputted points. It is assumed that the 
    inital list of points and the function tester.evaluateAtX agree to tolerance tester.relativeTolerance."""

    def thickenXYList2( xl, yl, xu, yu, newList, tester, level ) :

        if( level == biSectionMax ) : return
        level += 1
        xMid = 0.5  * ( xl + xu )
        yMid = 0.5  * ( yl + yu )
        y = tester.evaluateAtX( xMid )
        dy = abs( y - yMid )
        if( ( dy > abs( y * tester.relativeTolerance ) ) and ( dy > tester.absoluteTolerance ) ) :
            newList.append( [ xMid, y ] )
            thickenXYList2( xl, yl, xMid, y, newList, tester, level )
            thickenXYList2( xMid, y, xu, yu, newList, tester, level )

    if( len( list ) < 2 ) : raise Exception( "len( list ) = %2 < 2" % len( list ) )
    newList = []
    for i, xy in enumerate( list ) :
        x2, y2 = xy
        if( i > 0 ) : thickenXYList2( x1, y1, x2, y2, newList, tester, 0 )
        newList.append( [ x2, y2 ] )
        x1, y1 = x2, y2
    newList.sort( )
    return( newList )

class fastSumOfManyAddends :
    """This class in designed to sum a lot of endl2dmath or fudge2dmath object together efficiently. For example,
    consider the list f2d of 100,000 fudge2dmath objects that are to be summed. One way to do this is as
    s = fudge2dmath( )
    for f in f2d : s = s + f

    In general, this is very inefficient and will take a long time. Using, this class as

    fs = fastSumOfManyAddends( )
    for f in f2d : fs.appendAddend( f )
    s = fs.returnSum( )

    is, in general, much more efficient (i.e., runs a lot faster) and it should never be less efficient.

    While this class was designed for endl2dmath and fudge2dmath objects, it should work for any object
    for which the '+' operation is defined."""

    def __init__( self ) :
        """Constructor for fastSumOfManyAddends."""

        self.clear( )

    def appendAddend( self, addend ) :
        """Adds addend to current sum efficiently."""

        n = len( self.list )
        for i in xrange( n ) :
            if( self.list[i] == None ) :
                self.list[i] = addend
                addend = None
                break
            else :
                addend = addend + self.list[i]
                self.list[i] = None
        if( addend != None ) : self.list.append( addend )

    def clear( self ) :
        """Clears currently summed data."""

        self.list = []

    def returnSum( self ) :
        """Returns the current sum of all addends appended."""

        s = None
        for l in self.list :
            if( l != None ) :
                if( s == None ) :
                    s = l
                else :
                    s = s + l
        return( s )

def checkForNaN( v, str, printErrors = True, indentation = "", messages = None ) :

    if( v != v ) :
        s = "%s: value is nan" % str
        if( type( messages ) != type( None ) ) : messages.append( s )
        if( printErrors ) : printWarning( '\n'.join( checkMessagesToString( s, indentation = indentation ) ) ) 

def isNumber( n ) :

    if( type( 1. ) == type( n ) ) : return( True )
    if( type( 1 ) == type( n ) ) : return( True )
    if( type( numpyFloat64 ) == type( n ) ) : return( True )
    return( False )

def checkNumber( v, str = "", printErrors = True, indentation = "", messages = None, maxAbsFloatValue = None ) :

    checkForNaN( v, str, printErrors = printErrors, indentation = indentation, messages = messages )
    if( not( maxAbsFloatValue is None ) ) :
        if( abs( v ) > maxAbsFloatValue ) :
            s = "%s: abs of number %s exceeds maxAbsFloatValue = %s" % ( str, v, maxAbsFloatValue )
            if( not( messages is None ) ) : messages.append( s )
            if( printErrors ) : printWarning( '\n'.join( checkMessagesToString( s, indentation = indentation ) ) )

def getValue( n ) :

    if( isNumber( n ) ) : return( n )
    if( isinstance( n, physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty ) ) : return( n.getValue( ) )
    raise Exception( 'Invalue number object = %s' % brb.getType( n ) )
