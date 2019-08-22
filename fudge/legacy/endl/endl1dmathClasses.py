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

"""
This module contains the endl1dmath class. Also see the supporting module endl1dmathmisc.py.
"""

import os, math, copy, types
from fudge.core import fudgemisc
import endl1dmathmisc
from fudge.core.utilities import fudgeFileMisc, subprocessing
from fudge.vis.gnuplot import plotbase

__metaclass__ = type

class endl1dmath :
    """
    This class is designed to support operations on a list of numbers (e.g., [ number, number, ..., number ]) that may be useful to scientists.
    A number is an integer or a float.

    Members ::

        data        | A python list of numbers.
        columns     | This is always 1 as this is 1d data.
        yLabel      | Labelled used as the y-Label when plotting.

    Some examples are,

    >>> my1d = endl1dmath( [ 1, 2, 3 ] )
    >>> print my1d
     1.0000000e+00
     2.0000000e+00
     3.0000000e+00
    >>> new1d = 2. * my1d + 3.14
    >>> print new1d
     5.1400000e+00
     7.1400000e+00
     9.1400000e+00
    >>> f1d = new1d.mapFunction( "math.exp(x)" )
    >>> print f1d
     1.7071577e+02
     1.2614284e+03
     9.3207651e+03
    """

    def __init__( self, data = None, checkDataType = False, yLabel = None, label = "unknown", toFloat = 0 ) :
        """Returns an endl1dmath object. Data must be a python list of numbers (e.g., [ number, number, ..., number ])."""

        if( data == None ) : data = []
        self.data = endl1dmathmisc.get1dmathData( data, "endl1dmath.__init__", "data" )
        if ( checkDataType ) : endl1dmathmisc.check1dData( data )
        self.columns = 1
        self.yLabel = yLabel
        self.label = label
        if( toFloat ) : self.toFloat( )
        self.xInterpolation = 'flat'

    def __getitem__( self, i ) :
        """Returns self's value at index i."""

        return self.data[i]

    def __setitem__( self, i, v ) :
        """Sets self's value at index i to v."""

        if ( type( v ) != type( 1 ) and type( v ) != type( 1. ) ) : raise Exception( "\nError in endl1dmath.__setitem__: invalid data type." )
        self.data[i] = v

    def __len__( self ) :
        "Returns the number of points in self's data (i.e., return len( self.data ))."

        return len( self.data )

    def __repr__( self ) :
        "Returns a printable string of self's data. Uses endl1dmathmisc.endl1d_repr_xFormat to convert each point to a string."

        s = ""
        for v in self.data : s = s + endl1dmathmisc.endl1d_repr_xFormat % v + "\n"
        return s

    def __cmp__( self, other ) :

        data = endl1dmathmisc.valid1dClassType( other, "endl1dmath.__cmp__", "comparison" )
        if( self.data < other.data ) : return( -1 )
        if( self.data == other.data ) : return( 0 )
        return( 1 )

    def __neg__( self ) :
        "Returns an endl1dmath object whose points are the negation of self's points."

        d = []
        for i in range( len( self.data ) ) : d.append( -self.data[i] )
        return endl1dmath( d, checkDataType = 0, yLabel = self.yLabel )

    def __add__( self, v ) :
        """Returns an endl1dmath object whose points are self's points plus v.  v can be a number or
        an endl1dmath object that is the same length as self."""

        d = []
        if ( type( v ) == type( 1 ) ) or ( type( v ) == type( 1. ) ) :
            for p in self.data : d.append( p + v )
        else :
            data = endl1dmathmisc.valid1dClassType( v, "endl1dmath.__add__", "addition" )
            if ( len( self.data ) != len( data ) ) : raise Exception( "\nError in endl1dmath.__add__: data lengths differ." )
            i = 0
            for p in self.data :
                d.append( p + data[i] )
                i += 1
        return endl1dmath( d, checkDataType = 0, yLabel = self.yLabel )

    def __radd__( self, v ) :
        """Same as __add__( self, v )."""

        return self + v

    def __sub__( self, v ) :
        """Returns an endl1dmath object that is self's points minus v.  v can be a number or
        an endl1dmath object that must be the same length as self."""

        d = []
        if ( type( v ) == type( 1 ) ) or ( type( v ) == type( 1. ) ) :
            for p in self.data : d.append( p - v )
        else :
            data = endl1dmathmisc.valid1dClassType( v, "endl1dmath.__sub__", "subtraction" )
            if ( len( self.data ) != len( data ) ) : raise Exception( "\nError in endl1dmath.__add__: data lengths differ." )
            i = 0
            for p in self.data :
                d.append( p - data[i] )
                i += 1
        return endl1dmath( d, checkDataType = 0, yLabel = self.yLabel )

    def __rsub__( self, v ) :
        """Returns an endl1dmath object that is v minus self's points.  v must be a number."""

        d = []
        if ( type( v ) == type( 1 ) ) or ( type( v ) == type( 1. ) ) :
            for p in self.data : d.append( v - p )
        else :
            raise Exception( "\nError in endl1dmath.__rsub__: invalid type for subtraction" )
        return endl1dmath( d, checkDataType = 0, yLabel = self.yLabel )

    def __mul__( self, v ) :
        """Returns an endl1dmath object that is self's points multiplied by v. v can be a number or
        an endl1dmath object that must be the same length as self."""

        d = []
        if ( type( v ) == type( 1 ) ) or ( type( v ) == type( 1. ) ) :
            v = float( v )
            for p in self.data : d.append( p * v )
        else :
            data = endl1dmathmisc.valid1dClassType( v, "endl1dmath.__mul__", "multiplication" )
            if ( len( self.data ) != len( data ) ) : raise Exception( "\nError in endl1dmath.__mul__: data lengths differ." )
            i = 0
            for p in self.data :
                d.append( p * data[i] )
                i += 1
        return endl1dmath( d, checkDataType = 0, yLabel = self.yLabel )

    def __rmul__( self, v ) :
        """Same as __mul__( self, v )."""

        return self.__mul__( v )

    def __div__( self, v ) :
        """Returns an endl1dmath object that is self's points divided by v.  v can be a number or
        an endl1dmath object that must be the same length as self."""

        d = []
        if ( type( v ) == type( 1 ) ) or ( type( v ) == type( 1. ) ) :
            v = float( v )
            for p in self.data : d.append( p / v )
        else :
            data = endl1dmathmisc.valid1dClassType( v, "endl1dmath.__div__", "division" )
            if ( len( self.data ) != len( data ) ) : raise Exception( "\nError in endl1dmath.__div__: data lengths differ." )
            i = 0
            for p in self.data :
                d.append( p / data[i] )
                i += 1
        return endl1dmath( d, checkDataType = 0, yLabel = self.yLabel )

    def check( self ) :
        "Check that self's data is a list of numbers."

        endl1dmathmisc.check1dData( self.data )

    def copyData( self ) :
        """Returns an endl1dmath object that is a copy, and not a reference, of self."""

        d = []
        for x in self.data : d.append( x )
        return endl1dmath( d, checkDataType = 0, yLabel = self.yLabel )

    def getDimensions( self ) :
        """Returns the dimensions (1 for endl1dmath) for this type of data."""

        return( 1 )

    def mapFunction( self, function ) :
        """
        Maps 'function' onto every point of self's data. 'Function' must be a string of the form 'f(x)' where
        f(x) evaluates to a number (e.g. 'math.sin(x)'). That is, self.data[i] = f( self.data[i] ).  For example, to
        exponentiate the points of instance my1d enter

        >>> my1d.mapFunction( 'exp(x)' )
        """

        d = []
        for x in self.data :
            exec( 'd.append( ' + function + ' )' )
        return endl1dmath( d, checkDataType = 1, yLabel = self.yLabel )

    def plot( self, xMin = None, xMax = None, yMin = None, yMax = None, xylog = 0, xLabel = None, yLabel = None, title = None ) :
        """
        This routine is like qplot (quick plot) except it spawns an interactive plot.  qplot is faster while plot is more flexible.

        ===== =============
        xylog   plot-type
        ===== =============
        0     linear-linear
        1     log-linear
        2     linear-log
        3     log-log
        ===== =============
        """

        if( xLabel == None ) : xLabel = 'indice'
        if( yLabel == None ) : yLabel = self.yLabel
        dt = plotbase.parsePlotOptions( xMin, xMax, yMin, yMax, xLabel, yLabel, title )
        f = fudgeFileMisc.fudgeTempFile( )
        i = 0
        for v in self.data :
            f.write( "%d %14.8e\n" % ( i, v ) )
            i += 1
        f.close( )
        p = os.path.join( __file__.split( '/fudge/core/' )[0], "fudge", "vis", "gnuplot", "endl2dplot.py" )
        s = [ 'python', p, 'xylog', str( xylog ) ] + dt + [ f.getName( ) ]
        subprocessing.spawn( s )

    def qplot( self, xMin = None, xMax = None, yMin = None, yMax = None, xylog = 0, xLabel = None, yLabel = None, title = None ) :
        """
        Also see plot( ).

        ===== =============
        xylog   plot-type
        ===== =============
        0     linear-linear
        1     log-linear
        2     linear-log
        3     log-log
        ===== =============
        """
        import Gnuplot
        xylog = int( xylog )            # Allow argument to be a string
        if( xMin != None ) : xMin = float( xMin )
        if( xMax != None ) : xMax = float( xMax )
        if( yMin != None ) : yMin = float( yMin )
        if( yMax != None ) : yMax = float( yMax )

        self.g = Gnuplot.Gnuplot( )
        self.g( 'set style data linespoints' )
        xylog = xylog % 4
        if ( len( self.data ) > 0 ) :
            my = self.data[0]
            for d in self.data : my = min( my, d )
            if ( my <= 0. ) and ( xylog > 1 ) : xylog = xylog - 2
        if   ( xylog == 1 ) : self.g( 'set logscale x' )
        elif ( xylog == 2 ) : self.g( 'set logscale y' )
        elif ( xylog == 3 ) : self.g( 'set logscale xy' )
        if ( xMin != None ) or ( xMax != None ) :
            xMin = `xMin`
            if ( xMin == "None" ) : xMin = "*"
            xMax = `xMax`
            if ( xMax == "None" ) : xMax = "*"
            self.g( "set xrange [ %s : %s ]" % ( xMin, xMax ) )
        if ( yMin != None ) or ( yMax != None ) :
            yMin = `yMin`
            if ( yMin == "None" ) : yMin = "*"
            yMax = `yMax`
            if ( yMax == "None" ) : yMax = "*"
            self.g( "set yrange [ %s : %s ]" % ( yMin, yMax ) )
        if ( xLabel != None ) : self.g( "set xlabel %s" % xLabel )
        if ( yLabel != None ) : self.g( "set ylabel %s" % yLabel )
        if ( title  != None ) : self.g( "set title %s" % `title` )
        d = []
        i = 0
        for v in self.data :
            d.append( [ i, v ] )
            i += 1
        self.g.plot( d )

    def set( self, other, duplicate = 0, checkDataType = 1 ) :
        """Sets self.data = other.data. If duplicate is true (i.e., not 0) than a copy of other's data 
        is made so that other's and self's data can be modified independently."""

        if( type( other ) == type( [] ) ) :
            data = other
        else :
            data = endl1dmathmisc.valid1dClassType( other, "endl1dmath.set", "setting" )
        if ( checkDataType ) : endl1dmathmisc.check1dData( data )
        if( duplicate ) :
            self.data = []
            for d in data : self.data.append( d )
        else :
            self.data = data

    def toFloat( self ) :
        "Converts every point of self to a float."

        map( float, self.data )

    def toString( self, format = None ) :
        "Returns the string returned by the endl1dmath's __repr__ function."

        endl1d_repr_xFormat = endl1dmathmisc.endl1d_repr_xFormat
        if( format != None ) : endl1dmathmisc.endl1d_repr_xFormat = format
        s = endl1dmath.__repr__( self )
        endl1dmathmisc.endl1d_repr_xFormat = endl1d_repr_xFormat
        return( s )

    def toStringWithPrefixSuffix( self, Prefix = "", Suffix = "" ) :
        "Returns a printable string of the points in self. Each point is a separate line surrounded by Prefix and Suffix."

        s = [ endl1dmathmisc.endl1d_repr_xFormat % x for x in self.data ]
        return fudgemisc.stringWithPrefixSuffix( s, Prefix, Suffix )
