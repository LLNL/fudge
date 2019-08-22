# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

"""
This module contains the endl2dmath class.  Also see the supporting module endl2dmathmisc.py.

It also contains the variables doSafeDivide and doSafeDivideAlways (see __div__ for details).
"""
#
#  Things to fix;
#   1) What to do if self.data is an empty list (= []).
#   2) Handle case when self.data has two identical E data (i.e., self.data[i][0] = self.data[i+1][0]).
#   3) endl2dmath.__init__ should do a deepcopy of data
#

import os
import math
from math import *  # Needed for scalex
import copy

from fudge.core import fudgemisc
from fudge.core.utilities import fudgeFileMisc, fudgeExceptions
from fudge.vis.gnuplot import plotbase
from fudge.core.utilities import subprocessing
from xData import XYs as XYsModule

__metaclass__ = type

numpyImported = False
try :
    import numpy
    numpyImported = True
except :
    fudgemisc.printWarning( "Warning from endl2dmathClasses.py: numpy not imported" )

import endlParameters
import endl1dmathClasses
import endl2dmathmisc

doSafeDivide = False
doSafeDivideAlways = False
safeDivideMaxValue = 1e20

class endl2dmath :
    """
    The class endl2dmath is the base class for all 2d data (i.e., data consisting of
    (x, y) pairs or 2 column data). The data is stored as a python list of [ x, y ] values.
    For example, the 2 column data::

      2.5  4.5
      3.1  3.6
      3.4  4.4
    is stored as self.data = [ [ 2.5, 4.5 ], [ 3.1, 3.6 ], [ 3.4, 4.4 ] ].

    Members::

        data               A python list of [ x, y ] pairs.
        columns            Always 2 (i.e., 2 column data).
        interpolation      Type of interpolation used (see below).
        xLabel             X label used for plotting.
        yLabel             Y label used for plotting.

    Interpolation values and meaning::

        xylog   plot-type
        ----------------------
          0     linear-linear
          1     log-linear
          2     linear-log
          3     log-log
    """

    def __init__( self, data = None, checkDataType = False, xLabel = None, yLabel = None, label = "unknown", interpolation = 0, 
        template = None, allowSameX = False ) :
        """Returns an endl2dmath instance. Data must be of type list[ number, number ]."""

        if( data is None ) : data = []
        if( ( interpolation < 0 ) or ( interpolation > 3 ) ) :
            raise Exception( "\nError in endl2dmath.__init__: Bad interpolation value = %s" % `interpolation` )
        self.data = endl2dmathmisc.get2dmathData( data, "endl2dmath.__init__", "data" )
        if ( checkDataType ) : dummy = endl2dmathmisc.check2dData( self.data, positiveY = 0, allowSameX = allowSameX )
        if( not( template is None ) ) : endl2dmathmisc.valid2dClassType( template, "endl2dmath.__init__", "template" )
        self.columns = 2
        self.label = label
        self.xLabel = xLabel
        if( xLabel is None ) :
            try :
                self.xLabel = template.xLabel
            except :
                pass
        self.yLabel = yLabel
        if( yLabel is None ) :
            try :
                self.yLabel = template.yLabel
            except :
                pass
        self.interpolation = interpolation
        if(   self.interpolation == 0 ) :
            self.xInterpolation = 'linear,linear'
        elif( self.interpolation == 1 ) :
            self.xInterpolation = 'log,linear'
        elif( self.interpolation == 2 ) :
            self.xInterpolation = 'linear,log'
        elif( self.interpolation == 3 ) :
            self.xInterpolation = 'log,log'
        if ( not( template is None ) and ( interpolation is None ) ) : self.interpolation = template.interpolation

    def __getitem__( self, i ) :
        """Returns self's y value evaluated at i.  If i is outside the range
    of self's i domain then 0. is returned. This will be changed in version 4 to return
    self.data[i] (i.e., the (i+1)^th element of self.data). Use getValue or getValueInDomain
    instead of this routine."""

        return self.getValue( i )

    def __setitem__( self, x, y ) :
        """Sets self's y value evaluated at x to y.  If there is no point at x 
        then one is added. x and y must be numbers. Like __getitem__ this will 
        also be changed in version 4. Use setValue instead of this routine."""

        self.setValue( x, y )

    def __len__( self ) :
        "Returns the number of (x,y) pairs in the self's data member."

        return len( self.data )

    def __repr__( self ) :
        """Returns a printable string of the data in self. This method uses endl2dmathmisc.endl2d_repr_xFormat
        and endl2dmathmisc.endl2d_repr_yFormat to convert each point to a string."""

        xy = "%s %s" % ( endl2dmathmisc.endl2d_repr_xFormat, endl2dmathmisc.endl2d_repr_yFormat )
        s = [ xy % ( x, y ) for x, y in self.data ]
        s = '\n'.join( s ) + '\n'
        return s

    def __neg__( self ) :
        "Returns an endl2dmath instance that has all y values of self negated."

        d = self.copyData( )
        for xy in d.data : xy[1] = -xy[1]
        return d

    def __abs__( self ) :
        "Returns an endl2dmath instance whose y values are the absolute value of self's y values."

        d = self.copyData( )
        for xy in d.data : xy[1] = math.fabs( xy[1] )
        return d

    def __add__( self, v ) :
        """Returns an endl2dmath instance that is first set to self.union( v ).
    Then the y values of self are set to self's y value plus v's y value evaluated at
    each union point (v must be a number or another endl2dmath instance)."""

        if ( type( v ) == type( 1 ) ) or ( type( v ) == type( 1. ) ) :
            d = self.copyData( )
            for p in d.data : p[1] = p[1] + v
        else :
            data = endl2dmathmisc.valid2dClassType( v, "endl2dmath.__add__", "addition" )
            if ( len( self ) < 1 ) :
                d = data.copyData( )
            elif ( len( data ) < 1 ) :
                d = self.copyData( )
            else :
                d = self.union( data )
                j = 0
                dArray = data.data
                x1, y1 = dArray[0]
                x2, y2 = dArray[1]
                xEnd  = dArray[-1][0]
                j = 1
                n = len( dArray )
                for p in d.data :
                    x = p[0]
                    if( x < x1 ) :
                        vp = 0.
                    elif( x > xEnd ) :
                        vp = 0.
                    else :
                        if( x2 < x ) : 
                            j += 1
                            x1, y1 = x2, y2
                            x2, y2 = dArray[j]
                        vp = endl2dmathmisc.interpolate2dPoints( self.interpolation, x, x1, y1, x2, y2 )
                    p[1] = p[1] + vp
                if  ( self.data[0][0]  == data.data[-1][0] ) :
                    x = self.data[0][0]
                    d.setValue( x, 0.5 * ( self.data[0][1] + data.data[-1][1] ) )
                    j = 0
                    dx = endlParameters.endlEpsx * x
                    xp = x - 2. * dx
                    for xy in d.data :
                        if( xy[0] >= xp ) : break
                        j += 1
                    xp += dx
                    if( d.data[j][0] == x ) :
                        d.setValue( xp, data.getValue( xp ) )
                        j += 1
                    xp = x + 2. * dx
                    j += 1
                    if( d.data[j][0] > xp ) :
                        xp -= dx
                        d.setValue( xp, self.getValue( xp ) )
                elif( self.data[-1][0] == data.data[0][0]  ) :
                    x = data.data[0][0]
                    d.setValue( x, 0.5 * ( self.data[-1][1] + data.data[0][1] ) )
                    j = 0
                    dx = endlParameters.endlEpsx * x
                    xp = x - 2. * dx
                    for xy in d.data :
                        if( xy[0] >= xp ) : break
                        j += 1
                    xp += dx
                    if( d.data[j][0] == x ) :
                        d.setValue( xp, self.getValue( xp ) )
                        j += 1
                    xp = x + 2. * dx
                    j += 1
                    if( d.data[j][0] > xp ) :
                        xp -= dx
                        d.setValue( xp, data.getValue( xp ) )
            d.interpolation = min( self.interpolation, data.interpolation )
        return d

    def __radd__( self, v ) :
        """Same as __add__( self, v )."""

        return self + v

    def __sub__( self, v ) :
        """Returns an endl2dmath instance that is first set to self.union( v ).
    Then the y values of self are set to self's y value minus v's y value evaluated at 
    each union point (v must be a number or another endl2dmath instance)."""

        if ( type( v ) == type( 1 ) ) or ( type( v ) == type( 1. ) ) :
            d = self.copyData( )
            for p in d.data : p[1] = p[1] - v
        else :
            data = endl2dmathmisc.valid2dClassType( v, "endl2dmath.__sub__", "subtraction" )
            if ( len( self ) < 1 ) :
                d = -data
            elif ( len( data ) < 1 ) :
                d = self.copyData( )
            else :
                d = self.union( data )
                j = 0
                for p in d.data :
                    vp, j = data.GetValueNIndex( p[0], j )
                    p[1] = p[1] - vp
            d.interpolation = min( self.interpolation, data.interpolation )
        return d

    def __rsub__( self, v ) :
        """Returns an endl2dmath instance that is the same as self except every
    y value is subtracted from v (v must be a number)."""

        if ( type( v ) == type( 1 ) ) or ( type( v ) == type( 1. ) ) :
            d = self.copyData( )
            for e in d.data : e[1] = v - e[1]
        else :
            raise Exception( "\nError in endl2dmath.__rsub__: invalid type for subtraction" )
        return d

    def __mul__( self, v ) :
        """Returns an endl2dmath instance that is the same as self except every
    y value is multiplied by v (v must be a number)."""

        if ( type( v ) == type( 1 ) ) or ( type( v ) == type( 1. ) ) :
            v = float( v )
            d = self.copyData( )
            for e in d.data : e[1] = v * e[1]
        else :
            data = endl2dmathmisc.valid2dClassType( v, "endl2dmath.__mul__", "multiplication" )
            if ( len( self ) < 1 ) :        # Nothing to multiply
                d = self.copyData( )
            elif ( len( data ) < 1 ) :
                d = data.copyData( )
            if ( len( self ) == 1 ) :       # Weird situation, but will treat (kind of) as a delta function.
                d = self.copyData( )
                d.data[0][1] = d.data[0][1] * data.getValue( d.data[0][0] )
            elif ( len( data ) == 1 ) :
                d = data.copyData( )
                d.data[0][1] = d.data[0][1] * self.getValue( d.data[0][0] )
            else :
                d = endl2dmath.union( self, data, xDomainUnionOnly = True )
                dArray = data.data
                if( len( d ) < 2 ) :
                    for xy in d.data : xy[1] = xy[1] * data.getValue( xy[0] )
                else :
                    j = 0
                    x0 = d.data[0][0]
                    for x, y in dArray :                # This is needed as the slice above may make x0 > dArray[0][0].
                        if( x >= x0 ) : break           # Should break on x == x0, as d contains all x-values in self and v.
                        j += 1
                    x1, y1 = dArray[j]
                    if( ( x1 > x0 ) and ( j > 0 ) ) :
                        j -= 1
                        x1, y1 = dArray[j]
                    j += 1
                    x2, y2 = dArray[j]
                    xEnd  = dArray[-1][0]
                    for p in d.data :
                        x = p[0]
                        if( x > xEnd ) :
                            vp = 0.
                        else :
                            if( x2 < x ) : 
                                j += 1
                                x1, y1 = x2, y2
                                x2, y2 = dArray[j]
                            vp = endl2dmathmisc.interpolate2dPoints( self.interpolation, x, x1, y1, x2, y2 )
                        p[1] = p[1] * vp
                    trimOnMultiply = True
                    if( hasattr( self, 'trimOnMultiply' ) ) : trimOnMultiply = self.trimOnMultiply
                    if( hasattr( v, 'trimOnMultiply' ) ) : trimOnMultiply = trimOnMultiply and v.trimOnMultiply
                    if( trimOnMultiply ) : d.trim( )
            d.interpolation = min( self.interpolation, data.interpolation )
        return d

    def __rmul__( self, v ) :
        """Same as __mul__( self, v )."""

        return self.__mul__( v )

    def __div__( self, v ) :
        """Returns an endl2dmath instance that is the same as self except every
    y value is divide by v (v must be a number or endl2dmath instance). If either the 
    variable endl2dmathClasses.doSafeDivide or endl2dmathClasses.doSafeDivideAlways 
    is True, then the y-value where a divide-by-zero occurs is replaced by 
    float( 'NaN' ). Otherwise, a raise is triggered. Note, endl2dmathClasses.doSafeDivide
    is reset to False everytime __div__ is called, while endl2dmathClasses.doSafeDivideAlways
    is not."""

        global doSafeDivide, doSafeDivideAlways
        _doSafeDivide = doSafeDivide or doSafeDivideAlways
        doSafeDivide = False
        if ( type( v ) == type( 1 ) ) or ( type( v ) == type( 1. ) ) :
            v = float( v )
            d = self.copyData( )
            for e in d.data : e[1] = e[1] / v
        else :
            data = endl2dmathmisc.valid2dClassType( v, "endl2dmath.__div__", "division" )
            if ( len( self ) < 1 ) :        # Nothing to divide
                d = self.copyData( )
            elif ( len( data ) < 1 ) :
                raise Exception( "\nError in endl2dmath.__div__: denominator is an empty object" )
            else :
                d = data.slicex( xMin = self.xMin( ), xMax = self.xMax( ) )
                d = self.union( d )
                j = 0
                for p in d.data :
                    vp, j = data.GetValueNIndex( p[0], j )
                    try :
                        p[1] = p[1] / vp
                    except ZeroDivisionError :
                        if( _doSafeDivide ) :
                            p[1] = float( 'NaN' )
                        else :
                            raise Exception( "\nError in endl2dmath.__div__: Divide by zero at x = %.12e" % p[0] )
                    except :
                        raise
            d.interpolation = min( self.interpolation, data.interpolation )
        return d

    def safeDivide( self, other, fMaxValue = 10. ) :
        """Like __div__ except it attemps to handle divide by removing the offending point and 
        adding near by points on either side. This routine also detects zeros in the denominator
        and adds near by points on either side. See __div__ for more information."""

        global doSafeDivide
        def getMaxValue( y ) :
            "For internal use only. Returns safeDivideMaxValue if argument is positive and -safeDivideMaxValue otherwise."

            if( y < 0. ) : return( -safeDivideMaxValue )
            return( safeDivideMaxValue )

        def addPoint( self, safeData, flag, x1, x2, newOther, yMax ) :
            "For internal use only. Adds a point near the zero of the denominator."

            eps = 1e-5
            if( abs( x2 - x1 ) < eps * ( abs( x1 ) + abs( x2 ) ) ) : return
            if( flag == 1 ) :   # Add point just after x1.
                x = x1 + eps * ( x2 - x1 )
            else :              # Add point just before x2.
                x = x2 - eps * ( x2 - x1 )
            yo = newOther.getValue( x )
            ys = self.getValue( x )
            if( yo == 0. ) :
                y = yMax
                if( ys < 0. ) : y = -yMax
            else :
                y = ys / yo
            safeData.append( [ x, y ] )

        if ( type( other ) == type( 1 ) ) or ( type( other ) == type( 1. ) ) :
            other = float( other )
            results = self.copyData( )
            if( doSafeDivide and ( other == 0. ) ) :
                for e in results.data :   e[1] = getMaxValue( e[1] )
            else :
                for e in results.data : e[1] = e[1] / other
        else :
            if( len( other ) == 0 ) :
                newOther = other
            else :
                other = other.copyData( )
                other.trim( )
                newOther = []
                y1 = other.data[0][1]
                for x2, y2 in other.data :
                    if( y2 * y1 < 0 ) :                                     # Add the point between x1 and x2 where other's y-value is zero.
                        xZero = x1 - y1 * ( x2 - x1 ) / float( y2 - y1 )
                        newOther.append( [ xZero, 0. ] )
                    newOther.append( [ x2, y2 ] )
                    x1 = x2
                    y1 = y2
            newOther = endl2dmath( newOther, checkDataType = 0, template = self )
            doSafeDivide = True                                             # Note, this is automatically reset to false in __div__
            results = self / newOther
            yMin = None
            for xy in results.data :                                        # Find yMin, yMax and set all invalid (i.e., 'nan') y-values to None.
                x = xy[0]
                y = xy[1]
                if( not ( ( y > -1. ) or ( y < 1. ) ) ) :                   # Here I am relying on a bug? in python.
                    xy[1] = None
                    y = None
                if( not( y is None ) ) :
                    if( yMin is None ) :
                        yMin = y
                        yMax = y
                    else :
                        yMin = min( yMin, y )
                        yMax = max( yMax, y )
            yMax = max( abs( yMin ), abs( yMax ) ) * fMaxValue
            safeData = []
            if( len( results ) ) :
                x1, y1 = results.data[0]
                mode = 0
                if( y1 is None ) : mode = 1
                for x2, y2 in results.data :
                    if( y2 is None ) :
                        if( mode == 0 ) : addPoint( self, safeData, -1, x1, x2, newOther, yMax )  # Add point just before x2.
                        mode = 1
                    else :
                        if( mode == 1 ) : addPoint( self, safeData,  1, x1, x2, newOther, yMax )  # Add point just after x1.
                        mode = 0
                        safeData.append( [ x2, y2 ] )
                    x1 = x2
            results = endl2dmath( safeData, checkDataType = 0, template = results )
        return( results )

    def __rdiv__( self, v ) :
        """Returns an endl2dmath instance that is the same as self except every
    y value is v divided by y value (v must be a number)."""


        if ( type( v ) == type( 1 ) ) or ( type( v ) == type( 1. ) ) :
            v = float( v )
            d = self.copyData( )
            for e in d.data : e[1] = v / e[1]
        else :
            raise Exception( "\nError in endl2dmath.__rdiv__: numerator must be a number." )

    def __pow__( self, p ) :
        "Returns an endl2dmath instance whose y values are self's y values raised to the power p."

        d = self.copyData( )
        for xy in d.data : xy[1] = math.pow( xy[1], p )
        return d

    def __rpow__( self, p ) :
        "Returns an endl2dmath instance whose y values are p raised to the power of self's y values."

        d = self.copyData( )
        for xy in d.data : xy[1] = math.pow( p, xy[1] )
        return d

    def GetValueNIndex( self, x, i = 0 ) :
        """Returns the tuple ( y, index ) where y is the value of self evaluated at x and
    index is the starting index to be used in the next call when ascending the list.
    I is the starting index."""

        try :
            x = float( x )
        except :
            raise Exception( "\nError in endl2dmath.GetValueNIndex: Index x not a number" )
        if ( type( i ) == type( 1 ) ) :
            y = 0.
            if ( len( self ) > 0 ) :
                if ( self.data[0][0] <= x <= self.data[-1][0] ) :
                    l = len( self )
                    while ( i < l ) :
                        if ( x <= self.data[i][0] ) : break
                        i = i + 1
                    if ( x == self.data[i][0] ) :
                        y = self.data[i][1]
                    else :
                        xy1 = self.data[i-1]
                        x1 = xy1[0]
                        y1 = xy1[1]
                        xy2 = self.data[i]
                        x2 = xy2[0]
                        y2 = xy2[1]
                        y = endl2dmathmisc.interpolate2dPoints( self.interpolation, x, x1, y1, x2, y2 )
            return ( y, i )
        else :
            raise Exception( "\nError in endl2dmath.GetValueNIndex: i not an integer" )

    def toString( self, format = None ) :
        "Returns the string returned by the endl2dmath's __repr__ function. This can be useful when endl2dmath is used as a base class."

        endl2d_repr_xFormat = endl2dmathmisc.endl2d_repr_xFormat
        endl2d_repr_yFormat = endl2dmathmisc.endl2d_repr_yFormat
        if( not( format is None ) ) :
            endl2dmathmisc.endl2d_repr_xFormat = format
            endl2dmathmisc.endl2d_repr_yFormat = format
        s = endl2dmath.__repr__( self )
        endl2dmathmisc.endl2d_repr_xFormat = endl2d_repr_xFormat
        endl2dmathmisc.endl2d_repr_yFormat = endl2d_repr_yFormat
        return( s )

    def toStringWithPrefixSuffix( self, Prefix = "", Suffix = "" ) :
        "Returns a printable string of the data in self with Prefix and Suffix append to each line."

        xy = "%s %s" % ( endl2dmathmisc.endl2d_repr_xFormat, endl2dmathmisc.endl2d_repr_yFormat )
        s = [ xy % ( x, y ) for x, y in self.data ]
        return fudgemisc.stringWithPrefixSuffix( s, Prefix, Suffix )

    def clip( self, yMin, yMax ) :
        """Returns an endl2dmath object whose y values are self's y values clipped to be between yMin and yMax."""

        if( yMin >= yMax ) : raise Exception( "\nError in endl2dmath.clip: yMin must be less than yMax." )
        d = []
        Mode = 0
        for x, y in self.data :
            if(   y < yMin ) :
                yp = yMin
                if( ( Mode == 2 ) or ( Mode == 3 ) ) :
                    if( x != x1 ) :
                        xp = ( yp - y1 ) * ( x - x1 ) / ( y - y1 ) + x1
                        d.append( [ xp, yp ] )
                Mode = 1
            elif( y <= yMax ) :
                yp = y
                if  ( Mode == 1 ) :
                    if( x != x1 ) :
                        xp = ( yMin - y1 ) * ( x - x1 ) / ( y - y1 ) + x1
                        d.append( [ xp, yMin ] )
                elif( Mode == 3 ) :
                    if( x != x1 ) :
                        xp = ( yMax - y1 ) * ( x - x1 ) / ( y - y1 ) + x1
                        d.append( [ xp, yMax ] )
                Mode = 2
            else :
                yp = yMax
                if( ( Mode == 1 ) or ( Mode == 2 ) ) :
                    if( x != x1 ) :
                        xp = ( yp - y1 ) * ( x - x1 ) / ( y - y1 ) + x1
                        d.append( [ xp, yp ] )
                Mode = 3
            d.append( [ x, yp ] )
            x1 = x
            y1 = y
        return( endl2dmath( d, checkDataType = 0, template = self ) )

    def cmp( self, other, f = 5.e-5 ) :
        "Compares self to other and prints out data where they differ by more than f."

        l1 = len( self )
        l2 = len( other )
        i1, i2 = 0, 0
        while ( i1 < l1 ) and ( i2 < l2 ) :
            if ( self.data[i1][0] == other.data[i2][0] ) :
                d = self.relativediff( other, i1, i2 )
                if ( abs( d ) > f ) : print "i1 = %5d i2 = %5d   %.5e %.5e %.5e %+.5e %+.5e" % \
                (i1, i2, self.data[i1][0], self.data[i1][1], other.data[i2][1], self.data[i1][1] - other.data[i2][1], d)
                i1 = i1 + 1
                i2 = i2 + 1
            elif ( self.data[i1][0] < other.data[i2][0] ) :
                print "Extra E in  self: i1 = %5d %.5e %.5e" % ( i1, self.data[i1][0], self.data[i1][1] )
                i1 = i1 + 1
            else :
                print "Extra E in other: i2 = %5d %.5e %.5e" % ( i2, other.data[i2][0], other.data[i2][1] )
                i2 = i2 + 1

    def convolution( self, other, dxMax = None, convolutionMethod = None, interpolationAccuracy = 1e-2 ) :
        """Returns an endl2dmath instance that is a convolution of self and other. convolutionMethod
        can be None or 'pythonCodingV1'. If None and module fudgeConvolutions exist then function
        fudgeConvolutions.convolution will be used; otherwise, the slow, original python code will be 
        used. If interpolationAccuracy is not None then the endl2dmath.thin method is called with
        interpolationAccuracy = interpolationAccuracy."""

        def convoluteGH( start, end, x, g, h, hxMax, lastx, convolute ) :
            """Performs the convolution of g * h starting at x in g. The x-range of g must be greater
            than or equal to that of h."""

            dx = h[-1][0] - h[0][0]
            xEnd = x + dx
            if( xEnd > g[end][0] ) : xEnd = g[end][0]           # In case of round-off problem from previous line.
            while( g[end][0] < xEnd ) : end += 1
            value = 0.
            gx1, gy1 = g[start]
            start += 1
            gx2, gy2 = g[start]
            if( gx1 < x ) :
                gyp1 = ( gy2 * ( x - gx1 ) + gy1 * ( gx2 - x ) ) / float( gx2 - gx1 )
            else :
                gyp1 = gy1
            hx1 = None
            for hx2, hy2 in h :
                x2 = x + hx2
                if( x2 > xEnd ) : x2 = xEnd                     # In case of round-off problem from previous line.
                if( not( hx1 is None ) ) :
                    doMore = True
                    while( doMore ) :
                        if( gx2 < x2 ) :
                            xp2, gyp2 = gx2, gy2
                            hyp2 = ( hy2 * ( ( gx2 - x ) - hx1 ) + hy1 * ( hx2 - ( gx2 - x ) ) ) / float( hx2 - hx1 )
                            gx1, gy1 = gx2, gy2
                            if( gx2 < xEnd ) : start += 1
                            gx2, gy2 = g[start]
                        else :
                            doMore = False
                            xp2 = x2
                            if( gx2 == x2 ) :
                                gyp2 = gy2
                                gx1, gy1 = gx2, gy2
                                if( gx2 < xEnd ) : start += 1
                                gx2, gy2 = g[start]
                            else :
                                gyp2 = ( gy2 * ( x2 - gx1 ) + gy1 * ( gx2 - x2 ) ) / float( gx2 - gx1 )
                            hyp2 = hy2
                        value += ( xp2 - xp1 ) * ( ( gyp2 + gyp1 ) * ( hyp2 + hyp1 ) + gyp1 * hyp1 + gyp2 * hyp2 )
                        xp1 = xp2
                        gyp1 = gyp2
                        hyp1 = hyp2
                xp1 = x1 = x2
                hx1 = hx2
                hyp1 = hy1 = hy2
            value /= 6.
            convolute.append( [ x + hxMax, value ] )

        if( convolutionMethod not in [ None, 'pythonCodingV1' ] ) :
            fudgemisc.printWarning( 'Bad convolutionMethod argument to method endl2dmath.convolution in module endl2dmathClasses. ' + \
                'Using convolutionMethod = None' )
            convolutionMethod = None
        if( convolutionMethod is None ) :
            try :
                import fudgeConvolutions
            except :
                convolutionMethod = 'pythonCodingV1'
            else :
                if( dxMax is None ) : dxMax = 0.
                convolute = fudgeConvolutions.convolution( self.data, other.data, dxMax = dxMax )
        if( convolutionMethod == 'pythonCodingV1' ) :
            if( numpyImported == False ) : raise fudgeExceptions.ENDL_numpyException( 
                "numpy was not in python search path, cannot use endl2dmath.convolution." )
            convolute = []
            if( ( len( self ) > 1 ) and ( len( other ) > 1 ) ) :
                dx1 = self.xMax( ) - self.xMin( )
                dx2 = other.xMax( ) - other.xMin( )
                if( dx1 >= dx2 ) :
                    dx = dx2
                    g = self
                    h = other
                else :
                    dx = dx1
                    g = other
                    h = self
                if( dx > 0. ) :
                    h = numpy.array( h.data )
                    g = numpy.array( g.data )
                    hxMin = h[0][0]
                    hxMax = h[-1][0]
                    gxMin = g[0][0]
                    gxMax = g[-1][0]
                    up = 0
                    down = len( h ) - 1
                    while( up < down ) :                            # Reverse x-data as we need h(y-x).
                        x, y = h[down]
                        h[down] = h[up]
                        h[up] = x, y
                        up += 1
                        down -= 1
                    for xy in h : xy[0] = hxMax - xy[0]
                    if( dxMax is None ) : dxMax = dx / 50.
                    start = 0
                    n = len( g ) - 1
                    x = g[0][0]
                    end = 1
                    while( ( g[end][0] < x + dx ) and ( end < n ) ) : end += 1
                    convoluteGH( start, end, x, g, h, hxMax, x, convolute )
                    while( start < n ) :
                        lastx = x
                        x = min( x + dxMax, g[n][0] - dx )
                        if( abs( x - lastx ) <= 1e-10 * ( abs( x ) + abs( lastx ) ) ) : break
                        if( g[start+1][0] <= x ) :
                            start += 1
                            x = g[start][0]
                        if( lastx == x ) : continue
                        while( ( g[end][0] < x + dx ) and ( end < n ) ) : end += 1
                        convoluteGH( start, end, x, g, h, hxMax, lastx, convolute )
        convolute = endl2dmath( convolute, checkDataType = 0 )
        if( not( interpolationAccuracy is None ) ) :
            convolute = convolute.thin( interpolationAccuracy = interpolationAccuracy )
        return( convolute )

    def copyData( self ) :
        """Returns an endl2dmath instance that is a copy, and not a reference, of self."""

        d = []
        for xy in self.data : d.append( [ xy[0], xy[1] ] )
        return endl2dmath( d, checkDataType = 0, template = self )

    def exp( self ) :
        "Returns an endl2dmath instance that has all y values of self exponentiated."

        d = self.copyData( )
        for xy in d.data : xy[1] = math.exp( xy[1] )
        return d

    def getDimensions( self ) :
        """Returns the dimensions (2 for endl2dmath) for this type of data."""

        return( 2 )

    @property
    def dimension( self ) :

        return( self.getDimensions( ) )

    def domain( self ) :
        "Returns the domain (minumum and maximum x values) for self."

        return( self.xMin( ), self.xMax( ) )

    def getValue( self, x ) :
        """Returns self's y value evaluated at x.  If x is outside the range
        of self's x domain then 0. is returned. Also see getValueInDomain."""

        return self.GetValueNIndex( x )[0]

    def getValueInDomain( self, x ) :
        """Returns self's y value evaluated at x.  If x is outside the range
        of self's x domain then None is returned. Also see getValue."""

        if( ( len( self ) == 0 ) or ( x < self.xMin( ) ) or ( x > self.xMax( ) ) ) : return( None )
        return self.GetValueNIndex( x )[0]

    def integrateX( self, xArray = None, weight = None, normalize = 0 ) :
        """Returns an endl1dmath object of length len( xArray ) - 1 whose value at index i
    is the integral of self from xArray[i] to xArray[i+1].  If xArray is none then 
    the x values of self are used for xArray.  If Normailzie is non-zero the value 
    at index i is divided by the x interval xArray[i+1] - xArray[i]. If weight is not 
    equal to None, weigth must be an endl2dmath object, and the integral is of self * weight."""

        if( ( xArray is None ) and not( weight is None ) ) : xArray = self.xArray( )
        if( xArray is None ) :
            d = []
            ix = 0
            for xy in self.data :
                if ( ix != 0 ) :
                    if ( normalize ) :
                        d.append( 0.5 * ( xy[1] + y1 ) )
                    else :
                        d.append( 0.5 * ( xy[0] - x1 ) * ( xy[1] + y1 ) )
                x1 = xy[0]
                y1 = xy[1]
                ix += 1
        else :
            nx = len( xArray )
            if ( nx < 2 ) :
                d = []
            elif ( ( len( self.data ) < 2 ) or ( self.data[0][0] >= xArray[-1] ) or ( self.data[-1][0] <= xArray[0] ) ) :
                d = ( nx - 1 ) * [0.]
            else :
                if( weight is None ) :
                    ixy = []
                    for x in xArray : ixy.append( [ x, 0. ] )
                    ixy = endl2dmath( ixy, checkDataType = False )
                    s = self.copyData( )
                    if ( ( self.data[ 0][1] != 0. ) and ( self.data[ 0][0] > xArray[ 0] ) ) : s.data.insert( 0, [ self.data[ 0][0], 0. ] )
                    if ( ( self.data[-1][1] != 0. ) and ( self.data[-1][0] < xArray[-1] ) ) : s.data.append(    [ self.data[-1][0], 0. ] )
                    ixy = s.union( ixy )
                    iixy = 0
                    for xy in ixy.data :                    # Find index to start integration in ixy.
                        if ( xy[0] == xArray[0] ) : break
                        iixy += 1
                    xy1 = ixy.data[iixy]
                    iixy += 1
                    ix = 1
                    d = []
                    while ( ix < nx ) :                     # Do integration for each energy group.
                        x = xArray[ix]
                        sum = 0
                        while ( xy1[0] < x ) :
                            xy2 = ixy.data[iixy]
                            sum += ( xy2[0] - xy1[0] ) * ( xy2[1] + xy1[1] )
                            xy1 = xy2
                            iixy += 1
                        if ( normalize ) : sum /= ( xArray[ix] - xArray[ix-1] )
                        d.append( 0.5 * sum )
                        ix += 1
                else :
                    endl2dmathmisc.valid2dClassType( weight, "endl2dmath.integrateX", "weight" )
                    xMin = max( self.xMin( ), weight.xMin( ) )
                    xMax = min( self.xMax( ), weight.xMax( ) )
                    a = self.slicex( xMin = xMin, xMax = xMax )
                    c = weight.slicex( xMin = xMin, xMax = xMax )
                    a = a.union( c )
                    b = a.copyData( )
                    b.map( c )                              # self and weight are now a and b with the same x values.
                    ixy = []
                    for x in xArray : ixy.append( [ x, 0. ] )
                    ixy = endl2dmath( ixy, checkDataType = False )
                    s = self.copyData( )
                    if ( ( self.data[ 0][1] != 0. ) and ( self.data[ 0][0] > xArray[ 0] ) ) : s.data.insert( 0, [ self.data[ 0][0], 0. ] )
                    if ( ( self.data[-1][1] != 0. ) and ( self.data[-1][0] < xArray[-1] ) ) : s.data.append(    [ self.data[-1][0], 0. ] )
                    ixy = s.union( ixy )
                    iixy = 0
                    for xy in ixy.data :                    # Find index to start integration in ixy.
                        if ( xy[0] == xArray[0] ) : break
                        iixy += 1
                    xy1 = ixy.data[iixy]
                    iixy += 1
                    ix = 1
                    d = []
                    while ( ix < nx ) :                     # Do integration for each energy group.
                        x = xArray[ix]
                        sum = 0
                        while ( xy1[0] < x ) :
                            xy2 = ixy.data[iixy]
                            sum += ( xy2[0] - xy1[0] ) * ( xy2[1] + xy1[1] )
                            xy1 = xy2
                            iixy += 1
                        if ( normalize ) : sum /= ( xArray[ix] - xArray[ix-1] )
                        d.append( 0.5 * sum )
                        ix += 1
        return endl1dmathClasses.endl1dmath( d )

    def integrateOneFunction( self, xMin = None, xMax = None ) :
        """Returns the integral of self from xMin to xMax."""

        if( xMin is None ) : xMin = self.xMin( )
        if( xMax is None ) : xMax = self.xMax( )
        xMin = max( xMin, self.xMin( ) )
        xMax = min( xMax, self.xMax( ) )
        self_ = XYsModule.XYs( self.data, accuracy = 1e-5, axes = XYsModule.XYs.defaultAxes( ) )
        return( float( self_.integrate( xMin, xMax ) ) )

    def integrateTwoFunctions( self, other, xMin = None, xMax = None ) :
        """Returns the integral of self and other from xMin to xMax."""

        endl2dmathmisc.valid2dClassType( other, "endl2dmath.integrateTwoFunctions", "other" )
        if( xMin is None ) : xMin = self.xMin( )
        if( xMax is None ) : xMax = self.xMax( )
        xMin = max( xMin, self.xMin( ), other.xMin( ) )
        xMax = min( xMax, self.xMax( ), other.xMax( ) )
        self_ = XYsModule.XYs( self.data, accuracy = 1e-5, axes = XYsModule.XYs.defaultAxes( ) )
        other_ = XYsModule.XYs( other.data, accuracy = 1e-5, axes = XYsModule.XYs.defaultAxes( ) )
        return( self_.integrateTwoFunctions( other_, xMin, xMax ) )

    def integrateThreeFunctions( self, other1, other2, xMin = None, xMax = None ) :
        """Returns the integral of self, other1 and other2 from xMin to xMax."""

        endl2dmathmisc.valid2dClassType( other1, "endl2dmath.integrateThreeFunctions", "other1" )
        endl2dmathmisc.valid2dClassType( other2, "endl2dmath.integrateThreeFunctions", "other2" )
        if( xMin is None ) : xMin = self.xMin( )
        if( xMax is None ) : xMax = self.xMax( )
        xMin = max( xMin, self.xMin( ), other1.xMin( ), other1.xMin( ) )
        xMax = min( xMax, self.xMax( ), other1.xMax( ), other2.xMax( ) )
        self_ = XYsModule.XYs( self.data, accuracy = 1e-5, axes = XYsModule.XYs.defaultAxes( ) )
        other1_ = XYsModule.XYs( other1.data, accuracy = 1e-5, axes = XYsModule.XYs.defaultAxes( ) )
        other2_ = XYsModule.XYs( other2.data, accuracy = 1e-5, axes = XYsModule.XYs.defaultAxes( ) )
        return( self_.integrateThreeFunctions( other1_, other2_, xMin, xMax ) )

    def union( self, other, xDomainUnionOnly = False ) :
        """Returns an endl2dmath instance with its x values being a union of self's
    and other's x values.  The y values are self's y values mapped onto
    the union's x values. If xDomainUnionOnly is True then only the x domain
    common to self and other is returned."""

        d1 = self.data
        n1 = len( d1 )
        d2 = endl2dmathmisc.valid2dClassType( other, "endl2dmath.union", "intersecting" ).data
        n2 = len( d2 )
        d = []
        if( ( n1 != 0 ) or ( n2 != 0 ) ) :
            i1 = 0
            i2 = 0
            if( n1 != 0 ) :
                x1 = d1[i1][0]
            else :
                x1 = d2[i2][0]
            if( n2 != 0 ) :
                x2 = d2[i2][0]
            else :
                x2 = d1[i1][0]
            if( xDomainUnionOnly ) :
                if( x1 < x2 ) :
                    while( i1 < n1 ) :
                        x1 = d1[i1][0]
                        if( x2 <= x1 ) : break
                        i1 += 1
                elif( x1 > x2 ) :
                    while( i2 < n2 ) :
                        x2 = d2[i2][0]
                        if( x1 <= x2 ) : break
                        i2 += 1
            j = 0
            while( ( i1 < n1 ) and ( i2 < n2 ) ) :
                if( x1 <= x2 ) :
                    d.append( [ x1, d1[i1][1] ] )
                    if( x1 == x2 ) :
                        i2 += 1
                        if ( i2 < n2 ) : x2 = d2[i2][0]
                    i1 += 1
                    if( i1 < n1 ) : x1 = d1[i1][0]
                else :
                    ( y, j ) = self.GetValueNIndex( x2, j )
                    d.append( [ x2, y ] )
                    i2 += 1
                    if ( i2 < n2 ) : x2 = d2[i2][0]
            if( not xDomainUnionOnly ) :
                while( i1 < n1 ) :
                    d.append( [ d1[i1][0], d1[i1][1] ] )
                    i1 += 1
                while( i2 < n2 ) :
                    x2 = d2[i2][0]
                    ( y, j ) = self.GetValueNIndex( x2, j )
                    d.append( [ x2, y ] )
                    i2 += 1
        d = endl2dmath( d, checkDataType = 0, template = self )
        return d

    def log( self ) :
        "Returns an endl2dmath instance that has all y values of self log_e-ed."

        d = self.copyData( )
        for xy in d.data : xy[1] = math.log( xy[1] )
        return d

    def log10( self ) :
        "Returns an endl2dmath instance that has all y values of self log10_e-ed."

        d = self.copyData( )
        for xy in d.data : xy[1] = math.log10( xy[1] )
        return d

    def max( self ) :
        "Returns the maximum y value of self of None if self is empty."

        if( len( self ) == 0 ) :
            m = None
        else :
            m = self.data[0][1]
            for x, y in self.data :
                if( y > m ) : m = y
        return m

    def min( self ) :
        "Returns the minimum y value of self or None if self is empty."

        if( len( self ) == 0 ) :
            m = None
        else :
            m = self.data[0][1]
            for x, y in self.data :
                if( y < m ) : m = y
        return m

    def map( self, other ) :
        """Sets self's y values to other evaluated at self's x values."""

        j = 0
        for d in self.data : ( d[1], j ) = other.GetValueNIndex( d[0], j )

    def label( self ) :
        """Returns self's label."""

        return( self.label )

    def normalize( self ) :
        """Returns an endl2dmath instance that is the same as self normalized to 1. This, is
        the integral of dx y(x) is 1."""

        d = self.copyData( )
        sum = 0.
        x1 = None
        for x2, y2 in d.data :
            if( not( x1 is None ) ) : sum += ( x2 - x1 ) * ( y2 + y1 )
            x1 = x2
            y1 = y2
        sum /= 2.
        if( sum != 0 ) :
            for xy in d.data : xy[1] /= sum
        return( d )


    def overstrike( self, other, leftFill = 1, rightFill = 1 ) :
        """Returns an endl2dmath instance that is the same as self, except that the data from
    other.xMin to other.xMax is replaced by other's data.  Other must be an endl2dmath object.
    If leftFill is true then a point at other.xMin * ( 1. - endlParameters.endlEpsx ) is inserted 
    with y value of self at that x value.  If rightFill is true then a point at other.xMax * 
    ( 1. + endlParameters.endlEpsx ) is inserted  with y value of self at that x value."""

        d = []
        endl2dmathmisc.valid2dClassType( other, "endl2dmath.overstrike", "" )
        if ( len( other ) == 0 ) :
            for xy in self.data : d.append( [ xy[0], xy[1] ] )
        else :
            i = 0
            n = len( self.data )
            xMin = other.data[0][0] * ( 1. - endlParameters.endlEpsx )
            while ( i < n ) :                                   # Add self's points that are less than other's xMin.
                xy = self.data[i]
                if ( xy[0] >= xMin ) : break
                d.append( [ xy[0], xy[1] ] )
                i += 1
            if ( leftFill and ( i > 0 ) and ( i < n ) ) : d.append( [ xMin, self.getValue( xMin ) ] )
            for xy in other.data : d.append( [ xy[0], xy[1] ] ) # Add in other's points
            xMax = other.data[-1][0] * ( 1. + endlParameters.endlEpsx )
            while ( i < n ) :
                xy = self.data[i]
                if ( xy[0] > xMax ) : break
                i += 1
            if ( rightFill and ( i > 0 ) and ( i < n ) ) : d.append( [ xMax, self.getValue( xMax ) ] )
            while ( i < n ) :                                   # Add self's points that are greater than other's xMax.
                xy = self.data[i]
                d.append( [ xy[0], xy[1] ] )
                i += 1
        return endl2dmath( d, checkDataType = 0, template = self )

    def plot( self, xMin = None, xMax = None, yMin = None, yMax = None, xylog = 0, xLabel = None, yLabel = None, 
        title = None, style = "lines" ) :
        """
        This routine is like qplot (quick plot) except it spawns an interactive plot.
        qplot is faster while plot is more flexible.

        xylog interpolation values and meaning::

            xylog   plot-type
           -----------------------
              0     linear-linear
              1     log-linear
              2     linear-log
              3     log-log
        """

        if ( xLabel is None ) and not( self.xLabel is None ) : xLabel = self.xLabel
        if ( yLabel is None ) and not( self.yLabel is None ) : yLabel = self.yLabel
        dt = plotbase.parsePlotOptions( xMin, xMax, yMin, yMax, xLabel, yLabel, title )
        f = fudgeFileMisc.fudgeTempFile( )
        for p in self.data : f.write( "%15.7e %14.6e\n" % ( p[0], p[1] ) )
        f.close( )
        p = os.path.join( __file__.split( 'fudge/legacy/' )[0], "fudge", "vis", "gnuplot", "endl2dplot.py" )
        s = [ 'python', p, 'xylog', str( xylog ) ] + dt + [ f.getName( ) ]
        subprocessing.spawn( s )

    def qplot( self, xMin = None, xMax = None, yMin = None, yMax = None, xylog = 0, xLabel = None, yLabel = None,
        title = None, style = "lines" ) :
        """Also see plot( ).

        xylog interpolation values and meaning::

            xylog   plot-type
           -----------------------
              0     linear-linear
              1     log-linear
              2     linear-log
              3     log-log
        """

        import Gnuplot
        xylog = int( xylog )            # Allow argument to be a string
        if( not( xMin is None ) ) : xMin = float( xMin )
        if( not( xMax is None ) ) : xMax = float( xMax )
        if( not( yMin is None ) ) : yMin = float( yMin )
        if( not( yMax is None ) ) : yMax = float( yMax )

        self.g = Gnuplot.Gnuplot( )
        self.g( 'set style data %s' % style )
        xylog = xylog % 4
        if ( len( self.data ) > 0 ) :
            mx = self.data[0][0]
            my = self.data[0][1]
            for d in self.data :
                mx = min( mx, d[0] )
                my = min( my, d[1] )
            if ( mx <= 0. ) and ( xylog % 2 == 1 ) : xylog = xylog - 1
            if ( my <= 0. ) and ( xylog > 1 ) : xylog = xylog - 2
        if   ( xylog == 1 ) : self.g( 'set logscale x' )
        elif ( xylog == 2 ) : self.g( 'set logscale y' )
        elif ( xylog == 3 ) : self.g( 'set logscale xy' )
        if( not( xMin is None ) or not( xMax is None ) ) :
            xMin = `xMin`
            if ( xMin == "None" ) : xMin = "*"
            xMax = `xMax`
            if ( xMax == "None" ) : xMax = "*"
            self.g( "set xrange [ %s : %s ]" % ( xMin, xMax ) )
        if( not( yMin is None ) or not( yMax is None ) ) :
            yMin = `yMin`
            if ( yMin == "None" ) : yMin = "*"
            yMax = `yMax`
            if ( yMax == "None" ) : yMax = "*"
            self.g( "set yrange [ %s : %s ]" % ( yMin, yMax ) )
        if( ( xLabel is None ) and not( self.xLabel is None ) ) : xLabel = `self.xLabel`
        if( not( xLabel is None ) ) : self.g( "set xlabel %s" % `xLabel` )
        if( ( yLabel is None ) and not( self.yLabel is None ) ) : yLabel = `self.yLabel`
        if( not( yLabel is None ) ) : self.g( "set ylabel %s" % `yLabel` )
        if( not( title  is None ) ) : self.g( "set title %s" % `title` )
        self.g.plot( self.data )

    def relativediff( self, other, i1, i2 ) :
        """Compares self's y value at i1 to other's y value at i2.  Mainly for internal use."""

        d = ( self.data[i1][1] - other.data[i2][1] )
        if ( abs( d ) > 0 ) : d = 2 * d / ( abs( self.data[i1][1] ) + abs( other.data[i2][1] ) )
        return d

    def removeClosePoints( self, verbose = 0 ) :
        """Calls endl2dmathmisc.check2dData and removes all points that check2dData claims are too close."""

        ne, badXIndicies, messages = endl2dmathmisc.check2dData( self,  printErrors = False, printWarning = False )
        badXIndicies.reverse( )
        for i in badXIndicies :
            if( verbose ) : fudgemisc.printWarning( "    removing point at index %d with x value = %e" % ( i, self.data[i][0] ) )
            del self.data[i]

    def reversex( self ) :
        """Returns an endl2dmath instance whose list of [x,y] has been reversed. For example, the data
        [ [1, 2], [3, 3], [4, 5], [6, 3] ] becomes [ [6, 3], [4, 5], [3, 3], [1, 2] ]."""

        d = self.copyData( )
        d.data.reverse( )
        return( d )


    def scalex( self, f, xLabel = None ) :
        """Returns an endl2dmath instance whose x data is set to f(x) where the x in f(x)
        is self's x data and f(x) is a string representing a function. The y data is not altered.
        Two examples of f(x) are 'exp(x)' and 'sqrt( pow( x, 2 ) - 1 ) / x'. 
        Sets the returned objects xLabel to xLabel if is it not None."""

        d = self.copyData( )
        for xy in d.data :
            x = xy[0]
            xy[0] = eval( f )
        if( len( d ) ) :
            inverted = True
            x1 = d.data[0][0]
            for x2, y in d.data :
                if( x1 < x2 ) :
                    inverted = False
                    break;
            if( inverted ) : d = d.reversex( )
        if( not( xLabel is None ) ) : d.xLabel = xLabel
        return( d )

    def set( self, other ) :
        "Sets self's data to other's data."

        self.data = other.data

    def setValue( self, x, y ) :
        """Sets self's y value evaluated at x to y.  If there is no point at x 
        then one is added. x and y must be numbers."""

        np = [ x, y ]
        endl2dmathmisc.check2dPoint( np )
        if ( len( self.data ) == 0 ) :
            self.data = [ [ x, y ] ]
        else :
            if ( x <= self.data[-1][0] ) :
                i = 0
                for p in self.data :
                    if ( x <= p[0] ) :
                        if ( x == p[0] ) :
                            p[1] = y
                        else :
                            self.data.insert( i, np )
                        return
                    i = i + 1
            else :
                self.data.append( np )

    def slice( self, xMin = None, xMax = None ) :
        """Returns an endl2dmath instance that includes all points of self that lie between
    xMin and xMax inclusively.  If xMin = None (xMax = None) then the lowest (highest)
    x value of self is used for xMin (xMax)."""

        d = []
        if ( len( self.data ) > 0 ) :
            if( xMin is None ) : xMin = self.data[0][0]
            if( xMax is None ) : xMax = self.data[-1][0]
            i = 0
            imax = len( self )
            while ( i < imax ) and ( self.data[i][0] < xMin ) : i = i + 1
            while ( i < imax ) and ( self.data[i][0] <= xMax ) : d.append( self.data[i] ); i = i + 1
        return endl2dmath( d, checkDataType = 0, template = self )

    def sliceDull( self, xMin = None, xMax = None ) :
        """Same as slice except that the end points of the returned endl2dmath instance are 
    dulled if they are not the same as the end points of self."""

        s = self.slice( xMin, xMax )
        if ( len( s ) > 0 ) :
            if( s.data[ 0][0] > self.data[ 0][0] ) : endl2dmathmisc.dullLowerEdge2d( s )
            if( s.data[-1][0] < self.data[-1][0] ) : endl2dmathmisc.dullUpperEdge2d( s )
        return s
        

    def slicex( self, xMin = None, xMax = None, addZeroPoint = False ) :
        """Same as slice except if self does not contain a point at xMin (xMax) then it is added with the y 
    value being self evalauted at xMin (xMax). If self evalauted at xMin (xMax) is zero then no point is added,
    unless addZeroPoint is True."""

        if ( len( self ) > 0 ) :
            if( xMin is None ) : xMin = self.data[0][0]
            if( xMax is None ) : xMax = self.data[-1][0]
        d = self.slice( xMin, xMax )
        if ( len( self ) > 0 ) :
            if( ( xMin < self.xMax( ) ) and ( xMax > self.xMin( ) ) ) :
                if( ( len( d ) == 0 ) or ( d.data[0][0] > xMin ) ) :
                    v = self.getValue( xMin )
                    if( ( v != 0. ) or addZeroPoint ) : d.data.insert( 0, [ xMin, v ] )
                if( ( len( d ) == 0 ) or ( d.data[-1][0] < xMax ) ) :
                    v = self.getValue( xMax )
                    if( ( v != 0. ) or addZeroPoint ) : d.data.append( [ xMax, v ] )
        return d

    def slicexDull( self, xMin = None, xMax = None ) :
        """Same as slicex except the end points of the returned endl2dmath instance are 
    dulled if they are not the same as the end points of self."""

        s = self.slicex( xMin, xMax )
        if ( len( s ) > 0 ) :
            if( s.data[ 0][0] > self.data[ 0][0] ) : endl2dmathmisc.dullLowerEdge2d( s )
            if( s.data[-1][0] < self.data[-1][0] ) : endl2dmathmisc.dullUpperEdge2d( s )
        return s

    def sqrt( self ) :
        "Returns an endl2dmath instance that has all y values of self sqrt-ed."

        d = self.copyData( )
        for xy in d.data : xy[1] = math.sqrt( xy[1] )
        return d

    def thicken( self, dx, Log = 0 ) :
        """Returns an endl2dmath instance that is the same as self except that more
    points are added, if needed, to insure that the x spacing between points is less
    than dx if Log = 0 or the ratio between x values of points is less than dx if 
    Log != 0. The actual x spacing may be smaller than dx but not larger."""

        if( dx <= 0. ) : raise Exception( "\nError in endl2dmath.thicked: dx = %e <= 0." % dx )
        if( ( Log != 0 ) and ( dx <= 1. ) ) : raise Exception( "\nError in endl2dmath.thicked: dx = %e <= 1. for Log != 0" % dx )
        d = endl2dmath( [], checkDataType = 0, template = self )
        i = 0
        j = 0
        n = len( self )
        if( n > 0 ) :
            x2 = float( self.data[i][0] )
            ( y, j ) = self.GetValueNIndex( x2, j )
            d.data.append( [ x2, y ] )
            i += 1
        while( i < n ) :
            x1 = x2
            x2 = float( self.data[i][0] )
            if( Log == 0 ) :
                rm = ( x2 - x1 ) / dx
                m = int( rm )
                if( abs( rm - m ) > 1e-6 ) : m += 1
                if( m > 0 ) : s = ( x2 - x1 ) / m
            else :
                m = int( math.log( x2 / x1 ) / math.log( dx ) )
                if( m > 0 ) : s = math.exp( math.log( x2 / x1 ) / ( m + 1 ) )
            k = 1
            x = x1
            while( k <= m ) :
                if( Log == 0 ) :
                    x += s
                else :
                    x *= s
                if( abs( x - x2 ) <= 0.01 * dx ) : break
                ( y, j ) = self.GetValueNIndex( x, j )
                d.data.append( [ x, y ] )
                k += 1
            i += 1
            ( y, j ) = self.GetValueNIndex( x2, j )
            d.data.append( [ x2, y ] )
        return d

    def thin( self, interpolationAccuracy = 0.01, thinningMethod = None ) :
        """Returns an endl2dmath instance that is the same as self except the points are
    thinned while keeping linear interpolation accurate to interpolationAccuracy. interpolationAccuracy
    is forced into the range 1e-1 to 1e-8.  thinningMethod can be None, 'c_CodingV1' or 'pythonCodingV1'.
    If thinningMethod is None or 'c_CodingV1' then thin tries to use the thinning function in fudge2dThin.
    If this fails or thinningMethod is 'pythonCodingV1', a slower python code is used."""

        interpolationAccuracy = max( 1e-8, min( 1e-1, interpolationAccuracy ) )

        if( thinningMethod not in [ None, 'c_CodingV1', 'pythonCodingV1' ] ) :
            fudgemisc.printWarning( 'Bad thinningMethod argument to method endl2dmath.thin in module endl2dmathClasses. Using thinningMethod = None' )
            thinningMethod = None

        if( thinningMethod in [ None, 'c_CodingV1' ] ) :
            try :
                import fudge2dThin
            except :
                thinningMethod = 'pythonCodingV1'
            else :
                datap = fudge2dThin.fudge2dThin( self.data, interpolationAccuracy = interpolationAccuracy )

        if( thinningMethod == 'pythonCodingV1' ) :
            n = len( self )
            thinPoint = n * [ True ]

            i = -1
            x1 = None
            x2 = None
            multiPoints = 0
            for x3, y3 in self.data :       # Remove middle point if surrounding points have same x-value.
                if( not( x1 is None ) ) :
                    if( x1 == x2 == x3 ) :
                        thinPoint[i] == False
                        multiPoints += 1
                x1 = x2
                x2 = x3
                i += 1
            if( multiPoints == 0 ) :
                data = self.data
            else :
                data = []
                i = 0
                for xy in self.data :
                    if( thinPoint[i] ) : data.append( xy )
                    i += 1
                n = len( data )
                thinPoint = n * [ True ]

            if( n > 0 ) :                   # Do not thin end points.
                thinPoint[0] = False
                thinPoint[-1] = False

            i = -1
            y1 = None
            y2 = None
            for x3, y3 in data :            # Do not thin local minima and maxima.
                if( not( y1 is None ) and thinPoint[i] ) :
                    if( ( y2 - y1 ) * ( y3 - y2 ) < 0. ) : thinPoint[i] = False
                y1 = y2
                y2 = y3
                i += 1

            i = 0
            x1 = None
            for x2, y2 in data :            # Do not thin step function.
                if( x1 == x2 ) :
                    thinPoint[i-1] = False
                    thinPoint[i] = False
                x2 = x1
                i += 1

            iPrior = 0
            nm2 = n - 1
            while( iPrior < nm2 ) :           # Now thin the other points.
                x1, y1 = data[iPrior]
                iNext = iPrior + 1
                iCurrent = iNext
                while( thinPoint[iCurrent] ) :
                    iNext += 1
                    if( iNext > n ) : break
                    x3, y3 = data[iNext]
                    iCurrent = iPrior + 1
                    while( iCurrent != iNext ) :
                        x2, y2 = data[iCurrent]
                        if( abs( ( y3 - y1 ) * ( x2 - x1 ) - ( y2 - y1 ) * ( x3 - x1 ) ) >= abs( interpolationAccuracy * y2 * ( x3 - x1 ) ) ) :
                            i = iCurrent
                            f = 0.
                            while( i < iNext ) :
                                x2, y2 = data[i]
                                y2p = ( y3 - y1 ) * ( x2 - x1 ) / ( x3 - x1 ) + y1
                                if( y2 == y2p == 0. ) :
                                    iCurrent = i
                                else :
                                    fp = abs( y2p - y2 ) / ( abs( y2 ) + abs( y2p ) )
                                    if( fp > f ) :
                                        f = fp
                                        iCurrent = i
                                i += 1
                            thinPoint[iCurrent] = False
                        if( not thinPoint[iCurrent] ) : break
                        iCurrent += 1
                iPrior = iCurrent

            datap = []
            i = 0
            for x, y in data :
                if( not thinPoint[i] ) : datap.append( [ x, y ] )
                i += 1
        return( endl2dmath( datap, checkDataType = 0, template = self ) )

    def toInterpolation( self, interpolation, accuracy, diSectionMax = 3 ) :
        """Returns a new endl2dmath object with data converted to the specified interpolation and accuracy."""

        if( ( interpolation < 0 ) or ( interpolation > 3 ) ) : raise Exception( "Invalid interpolation = %d" % interpolation )
        if( self.interpolation != interpolation ) :
            if( self.interpolation == 3 ) :
                if( interpolation == 0 ) :
                    d = endl2dmathmisc.convertLogLogToLinLin( self, accuracy, diSectionMax = diSectionMax )
                else :
                    raise Exception( "Conversion of data from interpolation = %d to interpolation = %d is currently not supported" % \
                        ( self.interpolation, interpolation ) )
            else :
                raise Exception( "Conversion of data from interpolation = %d to interpolation = %d is currently not supported" % \
                    ( self.interpolation, interpolation ) )
        else :
            d = self.copyData( )
        return( d )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0, cls = None ) :
        """
        This routine is designed to make an endl2dmath instance work with plotting packages like plot2d.py. This routine is not 
        compatible with legacy fudge as the returned instance is of type XYs and not of type endl2dmath.
        """

        if( accuracy == None ) : accuracy = 1e-3
        linear = self.toInterpolation( 0, accuracy )
        linearXYs = XYsModule.XYs( linear.data, accuracy = accuracy, axes = XYsModule.XYs.defaultAxes( ) )
        if( cls is not None ) : linearXYs = cls.returnAsClass( linearXYs )
        return( linearXYs )

    def copyDataToXYs( self ) :

        return( self.data )

    def copyDataToXsAndYs( self ) :

        xs, ys = [], []
        for x, y in self.data :
            xs.append( x )
            ys.append( y )
        return( xs, ys )

    def trim( self ) :
        """Removes extra leading and trailing data whose y values are zero."""

        i = -1                                  # First, remove zeros from beginning.
        for e in self.data :
            if ( e[1] != 0. ) : break
            i = i + 1
        if ( i > 0 ) : self.data = self.data[i:]
        i = len( self.data ) - 1                # Now, remove zeros from end.
        if ( i >= 0 ) :
            j = i - 1
            while ( i >= 0 ) and ( self.data[i][1] == 0. ) : i = i - 1
            if ( i == -1 ) : self.data = []
            elif ( i < j ) : self.data = self.data[:i + 2]

    def vectorScale( self, vector ) :
        """Scales the y value of self at x by vector evaluated at x.  Only self's data within the x range of vector are modified."""

        endl2dmathmisc.valid2dClassType( vector, "endl2dmath.vectorScaleAdd", "vector" )
        d = self.copyData( )
        i = 0
        for xy in d.data :
            ( y, i ) = vector.GetValueNIndex( xy[0], i )
            xy[1] = y * xy[1]
        return d

    def vectorScaleAdd( self, vector ) :
        """Scales the y value of self at x by vector evaluated at x.
        Only self's data within the x range of vector are modified."""

        endl2dmathmisc.valid2dClassType( vector, "endl2dmath.vectorScaleAdd", "vector" )
        d = self.copyData( )
        i = 0
        xMax = vector.xMax( )
        if( not( xMax is None ) ) :
            for xy in d.data :
                if ( xy[0] > xMax ) : break
                ( y, i ) = vector.GetValueNIndex( xy[0], i )
                xy[1] += y * xy[1]
        return d

    def xArray( self ) :
        "Returns an endl1dmath instances of the x values in self."

        xa = []
        for xy in self.data : xa.append( xy[0] )
        return endl1dmathClasses.endl1dmath( xa )

    def xMax( self ) :
        "Returns the maximum x value of self if self contains data, otherwise it returns None."

        if ( len( self.data ) > 0 ) : return self.data[-1][0]
        return None

    def xMin( self ) :
        "Returns the minimum x value of self if self contains data, otherwise it returns None."

        if ( len( self.data ) > 0 ) : return self.data[0][0]
        return None

    def yMax( self ) :
        "Returns the maximum y value of self if self contains data, otherwise it returns None."

        return( self.yMinMax( )[1] )

    def yMin( self ) :
        "Returns the minimum y value of self if self contains data, otherwise it returns None."

        return( self.yMinMax( )[0] )

    def yMinMax( self ) :
        "Returns the tuple of the minimum and maximum y values of self if self contains data, otherwise it returns ( None, None )."

        the_min = None
        if ( len( self.data ) > 0 ) : the_min = self.data[0][1]
        the_max = the_min
        for x, y in self.data:
            the_min = min( the_min, y )
            the_max = max( the_max, y )
        return( the_min, the_max )
