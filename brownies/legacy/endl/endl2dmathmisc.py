# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains miscellaneous routines for supporting the endl2dmath class from module endl2dmathClasses.

Module's global variables::

    endl2d_repr_xFormat : Default format string used for converting an x-element of the data of an endl2dmath object to a string.
    endl2d_repr_yFormat : Default format string used for converting an y-element of the data of an endl2dmath object to a string.
"""

import math

from brownies.legacy.endl import fudgemisc
from brownies.legacy.endl import endlParameters
from fudge.core.math import fudgemath

endl2d_repr_xFormat = "%14.7e"
endl2d_repr_yFormat = "%12.5e"

def get2dmathData( object, callingRoutine, msg ) :
    """If object is a python list, then it is returned. Else, if object is a sub-class 
    of endl2dmath its data is return otherwise a raise is executed."""

    if( type( object ) == type( [] ) ) :
        data = object
    else :
        data = valid2dClassType( object, callingRoutine, msg ).data
    return( data )

def valid2dClassType( object, callingRoutine, msg ) :
    """Returns the first argument, object, if it is a subclass of the endl2dmath class; else, triggers a raise."""

    from . import endl2dmathClasses
    if( isinstance(object, endl2dmathClasses.endl2dmath)) : return(object)
    raise Exception( "\nError in %s: invalid type = %s for %s" % ( callingRoutine, type( object ), msg ) )

def check2dPoint( p ) :
    """check2dPoint( p )\n    Checks that p is the list [ number, number ]. A raise is generated if p is not."""

    if( ( type( p ) != type( [] ) ) or ( len( p ) != 2 ) or not( fudgemath.isNumber( p[0] ) ) or not( fudgemath.isNumber( p[1] ) ) ) :
        raise Exception( "\nError in check2dPoint: data point not list of [ number, number ]:", p )

def check2dData( data, allowNegativeX = True, allowSameX = False, allowZeroX = True, positiveY = True, printWarning = True, printErrors = True, 
    xCloseEps = None, formatCode = 12, maxAbsFloatValue = None ) :
    """Checks that data (or data.data if data is an endl2dmath instance) is a list containing
    valid endl2dmath points. If positiveY is true (default) then a negative y value
    generates a raise.  If positiveY is false then the number of negative y values
    is counted and the count is returned. The x values must be in increasing order
    else a raise is generated. If allowZeroX is False, a zero first x value will generate a 
    raise. Also all x values must be positive unless allowNegativeX
    is true. If printWarning is true than close x value warnings are printed. X values 
    x1 and x2 are considered close if x2 - x1 < xCloseEps * ( x1 + x2 ). If on input xCloseEps
    is None, the it is set to endlParameters.endlEpsx * pow( 0.1, formatCode - 14 ) / 3. 
    If on input formatCode is 12 then it is set to 14.
    This function returns an integer and an array. The integer is the number of negative
    y values and the array is a list of all indices with close x-values. If printWarning
    is 'True' than warnings are printed (warning will not execute a raise). If printErrors 
    is 'True' then errors are printed, and if at least one error exist, a raise will be 
    executed after all data has been checked."""

    points = get2dmathData( data, "check2dData", "" )
    l = len( points )
    w = 0                       # Number of close energy warnings.
    ne = 0                      # Number of points with negative y values.
    badXIndicies = []
    messages = []
    if( xCloseEps is None ) :
        if( formatCode == 12 ) : formatCode = 14
        xCloseEps = endlParameters.endlEpsx / 3. * pow(0.1, (formatCode - 14))
    if( l > 0 ) :
        pl = None
        pn = points[0]
        i = 0
        if( pn[0] == 0. ) and ( not allowZeroX ) :
            s = 'check2dData: zero x value at index = 0'
            messages.append( s )
        elif( ( pn[0] < 0. ) and ( not allowNegativeX ) ) :
            s = 'check2dData: negative x value at index = %d, x = %e' % ( i, pn[0] )
            messages.append( s )
            if( printErrors ) : fudgemisc.printWarning( '\n'.join( fudgemisc.checkMessagesToString( s, indentation = '      ' ) ) )
        while 1 :
            check2dPoint( pn )
            fudgemath.checkNumber( pn[0], "check2dData: x[%d]" % i, messages = messages, indentation = '      ', printErrors = printErrors,
                maxAbsFloatValue = maxAbsFloatValue )
            fudgemath.checkNumber( pn[1], "check2dData: y[%d]" % i, messages = messages, indentation = '      ', printErrors = printErrors,
                maxAbsFloatValue = maxAbsFloatValue )
            if( pl is not None ) :
                checkCloseness = True
                if( pl[0] >= pn[0] ) :
                    checkCloseness = False
                    if( not allowSameX ) :
                        s = 'check2dData: x value x[%d] >= x[%d] (%.8e >= %.8e)' % ( i - 1, i, pl[0], pn[0] )
                        messages.append( s )
                        if( printErrors ) : fudgemisc.printWarning( '\n'.join( fudgemisc.checkMessagesToString( s, indentation = '      ' ) ) )
                if( checkCloseness and ( pn[0] - pl[0] ) < ( ( pn[0] + pl[0] ) * xCloseEps ) ) :
                    sn = endl2d_repr_xFormat % pn[0]
                    sl = endl2d_repr_xFormat % pl[0]
                    if( sn[-4:] == 'e-00' ) : sn = sn[:-4] + 'e+00'
                    if( sl[-4:] == 'e-00' ) : sl = sl[:-4] + 'e+00'
                    if( sn == sl ) :
                        badXIndicies.append( i )
                        s = 'check2dData: x values %.16e and %.16e very close' % ( pl[0], pn[0] )
                        messages.append( s )
                        if( ( w < 2 ) and printWarning ) : fudgemisc.printWarning( '          Warning in %s' % s )
                        w = w + 1
            if( pn[1] < 0 ) :
                if( positiveY ) :
                    s = 'check2dData: negative y value at index = %d (x = %e y = %e)' % ( i, pn[0], pn[1] )
                    messages.append( s )
                    if( printErrors ) : fudgemisc.printWarning( '\n'.join( fudgemisc.checkMessagesToString( s, indentation = '      ' ) ) )
                ne = ne + 1
            i = i + 1
            if( i == l ) : break
            pl = pn
            pn = points[i]
        if( printWarning and ( w > 1 ) ) : fudgemisc.printWarning( "          Warning in check2dData: %d x values very close" % w )
    if( printErrors and ( ( len( messages ) - w ) > 0 ) ) : raise Exception( '          bad 2d data' )
    return ne, badXIndicies, messages

def interpolate2dPoints( interpolation, x, x1, y1, x2, y2 ) :
    """This function interpolates y at x where x1 <= x <= x2 using interpolation."""
    
    y = y1
    if( y1 != y2 ) :
        if( x1 == x2 ) :                                                  # A why to handle a step function.
            y = 0.5 * ( y1 + y2 )
        elif( x == x1 ) :
            pass
        elif( x == x2 ) :
            y = y2
        else : 
            if( interpolation == 0 ) :                                  # linear-linear interpolation.
                y = ( ( x - x1 ) * y2 + ( x2 - x ) * y1 ) / ( x2 - x1 )
            elif ( interpolation == 1 ) :                               # log-linear interpolation.
                y = ( y2 - y1 ) * math.log( x / x1 ) / math.log( x2 / x1 ) + y1
            elif ( interpolation == 2 ) :                               # linear-log interpolation.
                y = y1 * math.exp( math.log( y2 / y1 ) * ( x - x1 ) / ( x2 - x1 ) )
            elif ( interpolation == 3 ) :                               # log-log interpolation.
                y = y1 * math.exp( math.log( y2 / y1 ) * math.log( x / x1 ) / math.log( x2 / x1 ) )
            else :
                raise Exception( "\nError in endl2dmathmisc.interpolate_points: unsupported interpolation = %s" % repr(interpolation) )
    return( y )

def interpolate_X( data, x ) :
    """For internal use only."""

    l = len( data )
    if( l < 1 ) : return( 0. )
    y = 0.
    if( data[0][0] <= x <= data[-1][0] ) :
        i = 0
        while( i < l ) :
            if( x <= data[i][0] ) : break
            i += 1
        if( x == data[i][0] ) :
            y = data[i][1]
        else :
            d1 = data[i-1]
            d2 = data[i]
            x1 = d1[0];
            x2 = d2[0];
            y1 = d1[1];
            y2 = d2[1];
            if( x2 == x1 ) :
                y = y1
            else :
                y = ( y1 * ( x2 - x ) + y2 * ( x - x1 ) ) / ( x2 - x1 )
    return( y )

def dulledges2d( d ) :
    """Calls dullLowerEdge2d( d ) and dullUpperEdge2d( d )."""

    dullLowerEdge2d( d )
    dullUpperEdge2d( d )
    return d

def dullLowerEdge2d( d ) :
    """Dulls the lower edge of d if the y value at that edge is not 0.  d must be an endl2dmath instance."""

    if ( len( d.data ) == 0 ) or ( d.data[0][1] == 0. ) : return d
    x = d.data[0][0]
    if ( x > endlParameters.endlSmallx) :
        if ( len( d.data ) > 1 ) and (d.data[1][0] > (1. + 2. * endlParameters.endlEpsx) * x) :
            d.data[0][1] = d[(1. + endlParameters.endlEpsx) * x]
            d.data[0][0] = (1. + endlParameters.endlEpsx) * x
        d.data.insert(0, [(1 - endlParameters.endlEpsx) * x, 0.])
    return d 

def dullUpperEdge2d( d ) :
    """Dulls the upper edge of d if the y value at that edge is not 0.  d must be an endl2dmath instance."""

    if ( len( d.data ) == 0 ) or ( d.data[-1][1] == 0. ) : return d
    x = d.data[-1][0]
    if ( x < endlParameters.endlLargex) :
        if ( len( d.data ) > 1 ) and (d.data[-2][0] < (1. - 2. * endlParameters.endlEpsx) * x) :
            d.data[-1][1] = d[(1. - endlParameters.endlEpsx) * x]
            d.data[-1][0] = (1. - endlParameters.endlEpsx) * x
        d.data.append([(1 + endlParameters.endlEpsx) * x, 0.])
    return d

def gauss2d( xc, xw, f = 1e-3, xsigmax = 4, xMin = None, xMax = None ) :
    """Returns an endl2dmath instance the has a Gaussian shape centered at xc with
    width xw (i.e., Exp( -( (x - xc) / xw )^2 / 2 ).  The points are generated 
    such that linear interpolation is accurate to f.  The x value's range is limited 
    by xMin, xMax and abs(x - xc) / xw <= xsigmax."""

    from . import endl2dmathClasses
    def addpoints( xmin, xmax, x, a, f, DoPoints ) :

        fs = 1.
        if ( x > 1. ) : fs = -1.
        y = ( 1. + fs * f ) * math.exp( -x * x / 2. )
        xl = -x
        yl = y
        if ( DoPoints[0] ) :
            if ( xl < xmin ) :
                DoPoints[0] = 0
                if ( len( a ) > 0 ) :
                    yl = ( ( xmin - a[0][0] ) * yl +  ( xl - xmin ) * a[0][1] ) / float( xl - a[0][0] )
                    xl = xmin
            if ( xl >= xmin ) : a.insert( 0, [ xl, yl ] )
        if ( DoPoints[1] ) :
            if (  x > xmax ) :
                DoPoints[1] = 0
                if ( len( a ) > 0 ) :
                    y = ( ( xmax - a[-1][0] ) * y +  ( x - xmax ) * a[-1][1] ) / float( x - a[-1][0] )
                    x = xmax
            if (  x <= xmax ) : a.append( [  x, y ] )
    
    def nextpoint( xs_, xe_, f ) :
        xs = xs_
        xe = xe_
        fs = -1
        if ( xs < 1 ) : fs = 1
        x1 = xs
        y1 = math.exp( -x1 * x1 / 2. )
        for i in range( 40 ) :
            x = ( xs + xe ) / 2.
            if ( xe - x < 1e-14 * x ) : break
            s = ( math.exp( ( x1 - x ) * ( x1 + x ) / 2. ) - 1. ) / ( x - x1 )
            z = 0.5 * ( - 1. / s + x1 )
            xr = z - fs * math.sqrt( z * z - 1 )
            r = 1. - ( s * ( xr - x1 ) + 1. ) * y1 * math.exp( xr * xr / 2. )
            if( fs * r < 2. * f ) :
                xs = x
            else :
                xe = x
        return x

    f = max( min( 1e-2, f ), 1e-8 )
    if( xMin is None ) :
        xmin = -xsigmax
    else :
        xmin = ( xMin - xc ) / float( xw )
    if( xMax is None ) :
        xmax = xsigmax
    else :
        xmax = ( xMax - xc ) / float( xw )
    if( xmin < 0 ) and ( xmin < -xsigmax ) : xmin = -xsigmax
    if( xsigmax < xmax ) : xmax = xsigmax
    if ( xmin <= 0 ) and ( xmax >= 0 ) :
        xamin = 0
    else :
        xamin = min( abs( xmin ), abs( xmax ) )
    xamax = max( abs( xmin ), abs( xmax ) )
    DoPoints = [ 1, 1 ]
    a = []
    if ( xamin != 0. ) :                                                    # The curve does not go through x = 0.
        addpoints( xmin, xmax, xamin, a, f, DoPoints )
        x1 = xamin
    else :
        a.append( [0., 1.] )                                                # The center and the next point. 
        s = math.sqrt( 2. * ( math.sqrt( 1. + 6. * f )  - 1. ) / 3. )
        x0 = 2. * s / ( 1. - s * s )
        s = ( math.exp( -x0 * x0 / 2. ) -1. ) / x0
        c = s * ( s * x0 + 1. )
        b = x0 + c
        c = x0 * x0 / 2. + c * x0 + 0.5
        e = ( -b + math.sqrt( b * b + 4. * f * c ) ) / ( 2. * c )
        x1 = x0 + e
        addpoints( xmin, xmax, x1, a, f, DoPoints )

    e = pow( 3. * f / 2., 1. / 3. )                                 # Treat the points around x = 1 as special.
    e = pow( 3. * f / 2. / ( 1. - 9. * e / 8. + 63. * e * e / 40. ), 2. / 3. )
    xs = 0.
    xe = 1.
    for i in range( 40 ) :
        x = ( xs + xe ) / 2.
        if ( xe - x < 1e-14 * x ) : break
        d = 1. - math.exp( ( x * x - 1. ) / 2. ) * ( 1. - ( 1. - e ) * ( x - 1. ) ) + f
        if ( d < 0. ) :
            xs = x
        else :
            xe = x
    x = 1. - .9 * ( 1 - x )
    x1m = x                                                         # Point < 1.
    x1p = 1. + ( 1. - x )                                           # Point > 1.
    x = x1

    if ( x < x1m ) :
        while ( x < 1. ) :
            xs = x
            x = nextpoint( xs, x1m, f )
            if( abs( x1m - x ) < 1e-2 * ( x - xs ) ) : break
            addpoints( xmin, xmax, x, a, f, DoPoints )

    if ( x < x1m ) :
        addpoints( xmin, xmax, x1m, a, f, DoPoints )
        x = x1m
    if ( x < x1p ) :
        addpoints( xmin, xmax, x1p, a, f, DoPoints )
        x = x1p

    while ( x < xamax ) :
        x =  nextpoint( x, xamax + 1, f )
        addpoints( xmin, xmax, x, a, f, DoPoints )
    for e in a : e[0] = e[0] * xw + xc
    return endl2dmathClasses.endl2dmath(a, checkDataType = 0)

def gauss2dDull( xc, xw, f = 1e-3, xsigmax = 4, xMin = None, xMax = None ) :
    """Returns an endl2dmath instance of gauss2d( xc, xw, f, xsigmax, xMin, xMax ) that is dulled using dulledges2d( )."""

    return dulledges2d( gauss2d( xc, xw, f, xsigmax, xMin, xMax ) )

def point2d( x ) :
    """Returns an endl2dmath instance of [ x, 1. ]."""

    from . import endl2dmathClasses
    return endl2dmathClasses.endl2dmath([[x, 1.]], checkDataType = 0)

def point2dDull( x ) :
    """Returns an endl2dmath instance of point2d( x ) that is dulled using dulledges2d( )."""

    return dulledges2d( point2d( x ) )

def flattop2d( xMin, xMax ) :
    """Returns an endl2dmath instance of a flattop from xMin to xMax with height 1.  0 < xMin < xMax."""

    from . import endl2dmathClasses
    if( xMin >= xMax ) : raise Exception( "\nError in flattop2d: ( xMin < 0 ) or xMin >= xMax" )
    return endl2dmathClasses.endl2dmath([[xMin, 1.], [xMax, 1.]], checkDataType = 0)

def flattop2dDull( xMin, xMax ) :
    """Returns the endl2dmath instances returned by flattop2d( xMin, xMax ) and that is dulled using dulledges2d( )."""

    return dulledges2d( flattop2d( xMin, xMax ) )

def ramp2d( xMin, xMax, down = 0 ) :
    """Returns an endl2dmath instance of a ramp from xMin to xMax with height
    0 at xMin and 1 at xMax if down = 0, and 1 at xMin and 0 at xMax otherwise."""

    from . import endl2dmathClasses
    if( xMin >= xMax ) : raise Exception( "\nError in ramp2d: xMin >= xMax" )
    if( down != 0 ) :
        return endl2dmathClasses.endl2dmath([[xMin, 1.], [xMax, 0.]], checkDataType = 0)
    else :
        return endl2dmathClasses.endl2dmath([[xMin, 0.], [xMax, 1.]], checkDataType = 0)

def ramp2dDull( xMin, xMax, down = 0 ) :
    """Returns the endl2dmath instances returned by ramp2d( xMin, xMax, ramp ) and that is dulled
    using dulledges2d( )."""

    return dulledges2d( ramp2d( xMin, xMax, ramp ) )

def triangle2d( x0, x1, x2 ) : 
    """Returns an endl2dmath instance of a triangle that starts at x0 with height 0,
    ramps up to 1 at x1 and ramps down to 0 at x2."""

    from . import endl2dmathClasses
    if( ( x0 > x1 ) or ( x1 > x2 ) ) :
        raise Exception( "\nError in triangle2d: triangle points must have x0 <= x1 <= x2" )
    return endl2dmathClasses.endl2dmath([[x0, 0.], [x1, 1.], [x2, 0.]], checkDataType = 0)

def convertFunctionToLinLin( f, xMin, xMax, accuracy ):
    '''
    Returns a new endl2dmath object that represents the function f(x) as a lin-lin table to 
    the specified absolute accuracy on the interval xMin to xMax.
    '''

    from . import endl2dmathClasses
    def convertFunctionToLinLin2( f, xMin, yMin, xMax, yMax, accuracy, data, level ) :

        if( level > 16 ) : return
        xMid, yMid = 0.5 * ( xMin + xMax ), 0.5 * ( yMin + yMax )
        y = f( xMid )
        if( abs( yMid - y ) > abs( y * accuracy ) ) :
            convertFunctionToLinLin2( f, xMin, yMin, xMid, y, accuracy, data, level + 1 )
            data.append( [ xMid, y ] )
            convertFunctionToLinLin2( f, xMid, y, xMax, yMax, accuracy, data, level + 1 )

    data = [ [ xMin, f( xMin ) ] ]
    convertFunctionToLinLin2( f, xMin, data[0][1], xMax, f( xMax ), accuracy, data, 0 )
    data.append( [ xMax, f( xMax ) ] )
    return(endl2dmathClasses.endl2dmath(data = data))

def convertLogLogToLinLin( self, accuracy, diSectionMax = 3 ) :
    """Returns a new endl2dmath object of self converted from log-log to lin-lin to specified accuracy."""

    if( accuracy < 1e-4 ) : accuracy = 1e-4
    if( diSectionMax > 8 ) : diSectionMax = 8

    def convertLogLogToLinLin2( self, x1, y1, x2, y2, level ) :
        """For internal use."""

        level_ = level + 1
        if( level_ > diSectionMax ) : return
        u2 = x2 / x1
        v2 = y2 / y1
        logXs = math.log( u2 )
        logYs = math.log( v2 )
        logYXs = logYs / logXs
        u = logYXs * ( u2 - v2 ) / ( ( 1 - logYXs ) * ( v2 - 1 ) )
        vLog = math.pow( u, logYXs )
        vLin = ( u2 - u + v2 * ( u - 1 ) ) / ( u2 - 1 )
        if( abs( vLin / vLog  - 1. ) > accuracy ) :
            x = x1 * u
            y = y1 * vLog
            self.setValue( x, y )
            convertLogLogToLinLin2( self, x1, y1, x, y, level_ )
            convertLogLogToLinLin2( self, x, y, x2, y2, level_ )

    new = self.copyData( )
    x1 = None
    for x2, y2 in self.data :
        if( x1 is not None ) :
            if( ( x1 <= 0 ) or ( y1 <= 0 ) or ( x2 <= 0 ) or ( y2 <= 0 ) ) : 
                raise Exception( "One of x1 = %e y1 = %e x2 = %d or y2 = %e is <= 0." % ( x1, y1, x2, y2 ) )
            if( ( x1 != x2 ) and ( y1 != y2 ) ) : convertLogLogToLinLin2( new, x1, y1, x2, y2, 0 )
        x1, y1 = x2, y2
    return( new )
