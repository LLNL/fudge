# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the endl3dmath class.  Also see the supporting module endl3dmathmisc.py.
"""

import os
import sys

from fudge.core.utilities import brb, fudgeFileMisc
from LUPY import subprocessing
from fudge.vis.gnuplot import plotbase
from fudge.core.math import fudgemath
from brownies.legacy.endl import fudgemisc
from brownies.legacy.endl import endl2dmathmisc, endl3dmathmisc, endl1dmathClasses

__metaclass__ = type

class endl3dmath :
    """
    This class is designed to support operations on a list of 3-column number data that may be useful to scientists.
    (A number is an integer or a float.) In this documenation the three columns are referred to as x, y and z.
    The data is to be stored as a python list of [ x, yz-data-for-x ], where yz-data-for-x is a list of
    [ y, z ] data. For example, consider the 3 column data::

     1 1 1
     1 2 2
     1 3 1
     2 1 1
     2 2 0
     2 4 0
     3 2 0
     3 3 1
     3 4 2

    This data is stored using python lists as ::


     [ [ 1, [ [ 1, 1 ], [ 2, 2 ], [ 3, 1 ] ] ],
       [ 2, [ [ 1, 1 ], [ 2, 0 ], [ 4, 0 ] ] ]
       [ 3, [ [ 2, 0 ], [ 3, 1 ], [ 4, 2 ] ] ] ]

    Members::

            data        | A python list of numbers.
            columns     | This is always 3 as this is 1d data.
            xLabel      | Labelled used as the x-Label when plotting.
            yLabel      | Labelled used as the y-Label when plotting.
            zLabel      | Labelled used as the z-Label when plotting.
    """

    def __init__( self, data = None, checkDataType = False, xLabel = None, yLabel = None, zLabel = None, label = "unknown", interpolation = 0, template = None ) :
        """Creates an endl3dmath instance.  Data must be of type list[ number, list[ number, number ] ].
        The default plotting labels can be set with xLabel, yLabel and zLabel."""

        if( data is None ) : data = []
        if( interpolation != 0 ) : raise Exception( "\nError in endl3dmath.__init__: Bad interpolation value = %s" % repr(interpolation) )
        self.interpolation = interpolation
        self.data = endl3dmathmisc.get3dmathData(data, 'endl3dmath.__init__', 'data')
        if ( template is not None ) : endl3dmathmisc.valid3dClassType(template, "endl3dmath.__init__", "template")
        self.columns = 3
        self.xLabel = xLabel
        if ( ( template is not None ) and ( xLabel is None ) ) : self.xLabel = template.xLabel
        self.yLabel = yLabel
        if ( ( template is not None ) and ( yLabel is None ) ) : self.yLabel = template.yLabel
        self.zLabel = zLabel
        if ( ( template is not None ) and ( zLabel is None ) ) : self.zLabel = template.zLabel
        self.label = label
        if(   self.interpolation == 0 ) :
            self.xInterpolation = 'linear,linear,linear'
        else :
            raise Exception( 'Unsupported interpolation = %d for 3d data' % self.interpolation )

    def __getitem__( self, xy ) :
        """See getValue for current functioning. This method will be changed in version 4 to return
        self.data[i] (i.e., the (i+1)^th element of self.data)."""

        return( self.getValue( xy[0], xy[1] ) )

    def __len__( self ) :
        """Returns the number of items in the top list of self (i.e. the number
    of differce x values)."""

        return len( self.data )

    def __repr__( self ) :
        """Returns a printable string of the data in self.  This method uses endl3dmathmisc.endl3d_repr_xFormat, 
        endl3dmathmisc.endl3d_repr_yFormat and endl3dmathmisc.endl3d_repr_zFormat to convert each point to a string."""

        syz = "%%s %s %s" % (endl3dmathmisc.endl3d_repr_yFormat, endl3dmathmisc.endl3d_repr_zFormat)
        sxyz = []
        for x_etal in self.data :
            x = endl3dmathmisc.endl3d_repr_xFormat % x_etal[0]
            for etal in x_etal[1] : sxyz.append( syz % ( x, etal[0], etal[1] ) )
        s = '\n'.join( sxyz ) + '\n'
        return s

    def __add__( self, other ) :
        """Returns an endl3dmath instance that is the addition of self and other. If other is a number, then other is added
        to all z-values. If other is an endl3dmath instance, then 
            1) a union of the x-values from self and other is formed,
            2) for each x-value, the yz-values from self and other are obtained and added as two endl2dmath objects
                which forms the addition at x.
        Note that a renormalization of the returned instance may be needed for normalized distributions."""

        if( fudgemath.isNumber( other ) ) :
            other = float( other )
            d = self.copyData( )
            for x, yzs in d.data : 
                for yz in yzs : yz[1] = yz[1] + other
        else :
            data = endl3dmathmisc.valid3dClassType(other, "endl3dmath.__add__", "addition")
            if( len( self ) < 1 ) :
                d = data.copyData( )
            elif( len( data ) < 1 ) :
                d = self.copyData( )
            else :
                try :
                    unitBaseSelf = self.unitbase
                except :
                    unitBaseSelf = False
                try :
                    unitBaseOther = other.unitbase
                except :
                    unitBaseOther = False

                xs = [ x for x, yz in self.data ]
                for x, yz in other.data :
                    if( x not in xs ) : xs.append( x )
                xs.sort( )
                d = []
                for x in xs :
                    yz = self.getAtX( x, unitBase = unitBaseSelf, endl2dmathObject = True ) + other.getAtX( x, unitBase = unitBaseOther, 
                        endl2dmathObject = True )
                    d.append( [ x, yz.data ] )
                d = endl3dmath( d, interpolation = min( self.interpolation, data.interpolation ) )
        return( d )

    def __sub__( self, other ) :
        """Returns an endl3dmath instance that is the subtraction of other from self. If other is a number, then other is subracted
        from all z-values. If other is an endl3dmath instance, then 
            1) a union of the x-values from self and other is formed,
            2) for each x-value, the yz-values from self and other are obtained and subtracted as two endl2dmath objects.
        Note that a renormalization of the returned instance may be needed for normalized distributions."""

        if( fudgemath.isNumber( other ) ) :
            other = float( other )
            d = self.copyData( )
            for x, yzs in d.data : 
                for yz in yzs : yz[1] = yz[1] - other
        else :
            data = endl3dmathmisc.valid3dClassType(other, "endl3dmath.__sub__", "substraction")
            if( len( self ) < 1 ) :
                d = data.copyData( )
            elif( len( data ) < 1 ) :
                d = self.copyData( )
            else :
                try :
                    unitBaseSelf = self.unitbase
                except :
                    unitBaseSelf = False
                try :
                    unitBaseOther = other.unitbase
                except :
                    unitBaseOther = False

                xs = [ x for x, yz in self.data ]
                for x, yz in other.data :
                    if( x not in xs ) : xs.append( x )
                xs.sort( )
                d = []
                for x in xs :
                    yz = self.getAtX( x, unitBase = unitBaseSelf, endl2dmathObject = True ) - other.getAtX( x, unitBase = unitBaseOther, 
                        endl2dmathObject = True )
                    d.append( [ x, yz.data ] )
                d = endl3dmath( d, interpolation = min( self.interpolation, data.interpolation ) )
        return( d )

    def __mul__( self, other ) :
        """Returns an endl3dmath instance that is the multiplication of self by other. Currently, other must be a number."""

        if( fudgemath.isNumber( other ) ) :
            d = self.copyData( )
            for x, yzs in d.data : 
                for yz in yzs : yz[1] = yz[1] * other
        else :
            raise Exception( 'multiplier of an endl3dmath instance can only be a number, %s is not allowed' % brb.getType( other ) )

    def __div__( self, other ) :
        """Returns an endl3dmath instance that is the division of self by other. Currently, other must be a number."""

        if( fudgemath.isNumber( other ) ) :
            d = self.copyData( )
            for x, yzs in d.data : 
                for yz in yzs : yz[1] = yz[1] / other
        else :
            raise Exception( 'divisor of an endl3dmath instance can only be a number, %s is not allowed' % brb.getType( other ) )

    __truediv__ = __div__   # Python 3

    def copyData( self ) :
        """Returns an endl3dmath instance that is a copy, and not a reference, of self."""

        d3 = []
        for x_yz in self.data :
            d2 = []
            for yz in x_yz[1] : d2.append( [ yz[0], yz[1] ] )
            d3.append( [ x_yz[0], d2 ] )
        return endl3dmath( d3, checkDataType = 0, template = self )

    def copyDataToW_XYs( self, wUnit = None, xUnit = None, yUnit = None ) :
        """A function designed to work with fudgeMultiPlots.py. Not for general use."""

        import copy

        return( copy.copy( self.data ) )

    def getAtX(self, x, unitBase = False, extrapolation = endl3dmathmisc.noExtrapolation, endl2dmathObject = False) :
        """Calls endl3dmathmisc.interpolate3d to get self's yz data at x."""

        return(endl3dmathmisc.interpolate3d(x, self.data, unitBase = unitBase, extrapolation = extrapolation, endl2dmathObject = endl2dmathObject))

    def setAtX( self, x, data ) :
        """Sets self's yz data at x to data. Data can be a endl2dmath object."""

        if( ( type( x ) == type( 1 ) ) or ( type( x ) == type( 1. ) ) ) :
            points = endl2dmathmisc.get2dmathData(data, 'endl3dmathClasses.setAtX', 'points')
            i = 0
            for xyz in self.data :
                x_ = xyz[0]
                if( x == x_ ) :
                    xyz[1] = points
                    return
                elif( x < x_ ) :
                    break
                i += 1
            if( len( self.data ) == 0 ) :
                self.data = [ [ x, points ] ]
            else :
                self.data.insert( i, [ x, points ] )
        else :
            raise Exception( '\nError in endl3dmath.setAtX: Bad x value. Must be a number.' )

    def getDimensions( self ) :
        """Returns the dimensions (3 for endl3dmath) for this type of data."""

        return( 3 )

    @property
    def dimension( self ) :

        return( self.getDimensions( ) )

    def getValue( self, x, y = None ) :
        """
        Returns an endl2dmath instance of self evaluated at x if y is None. Otherwise, returns, as a float, the value of self at x, y.
        If x (y) is outside the range of self's x (y) domain then 0. is returned.
        """

        endl2dmath_atX = self.getAtX( x, endl2dmathObject = True )
        if( y is None ) : return( endl2dmath_atX )
        return( endl2dmath_atX.getValue( y ) )

    def toString( self, format = None ) :
        """Returns the string returned by the endl3dmath's __repr__ function. This can be useful when endl3dmath is used as a base class."""

        endl3d_repr_xFormat = endl3dmathmisc.endl3d_repr_xFormat
        endl3d_repr_yFormat = endl3dmathmisc.endl3d_repr_yFormat
        endl3d_repr_zFormat = endl3dmathmisc.endl3d_repr_zFormat
        if( format is not None ) :
            endl3dmathmisc.endl3d_repr_xFormat = format 
            endl3dmathmisc.endl3d_repr_yFormat = format
            endl3dmathmisc.endl3d_repr_zFormat = format
        s = endl3dmath.__repr__( self )
        endl3dmathmisc.endl3d_repr_xFormat = endl3d_repr_xFormat
        endl3dmathmisc.endl3d_repr_yFormat = endl3d_repr_yFormat
        endl3dmathmisc.endl3d_repr_zFormat = endl3d_repr_zFormat
        return( s )

    def toStringWithPrefixSuffix( self, Prefix = "", Suffix = "" ) :
        """Returns a printable string of the data in self with Prefix and Suffix append to each line.
        See __repr__ to change the output format."""

        syz = "%%s %s %s" % (endl3dmathmisc.endl3d_repr_yFormat, endl3dmathmisc.endl3d_repr_zFormat)
        sxyz = []
        for x_etal in self.data :
            x = endl3dmathmisc.endl3d_repr_xFormat % x_etal[0]
            for etal in x_etal[1] : sxyz.append( syz % ( x, etal[0], etal[1] ) )
        return fudgemisc.stringWithPrefixSuffix( sxyz, Prefix, Suffix )

    def plot( self, xMin = None, xMax = None, yMin = None, yMax = None, zMin = None, zMax = None, xyzlog = 0, 
        xLabel = None, yLabel = None, zLabel = None, title = None, xrot = None, zrot = None ) :
        """
        Plots the data.

        xyzlog values and meaning::

            xyzlog   plot-type for x-y-z axis
           -----------------------------------
              0     linear-linear-linear
              1     log-linear-linear
              2     linear-log-linear
              3     log-log-linear
              4     linear-linear-log
              5     log-linear-log
              6     linear-log-log
              7     log-log-log
        """

        if ( xLabel is None ) and ( self.xLabel is not None ) : xLabel = self.xLabel
        if ( yLabel is None ) and ( self.yLabel is not None ) : yLabel = self.yLabel
        if ( zLabel is None ) and ( self.zLabel is not None ) : zLabel = self.zLabel
        dt = plotbase.parsePlotOptions( xMin, xMax, yMin, yMax, xLabel, yLabel, title, zMin = zMin, zMax = zMax, zLabel = zLabel, \
            xrot = xrot, zrot = zrot )
        f = fudgeFileMisc.FudgeTempFile( )
        f.write( endl3dmath.__repr__( self ) )
        f.close( )
        p = os.path.join( __file__.split( 'fudge/legacy/' )[0], "fudge", "vis", "gnuplot", "endl3dplot.py" )
        python = sys.executable
        s = [ python, p, 'xyzlog', str( xyzlog ) ] + dt + [ f.getName( ) ]
        subprocessing.spawn( s )

    def xArray( self ) :
        """Returns an endl1dmath instances of the x values in self."""

        xa = []
        for xy in self.data : xa.append( xy[0] )
        return endl1dmathClasses.endl1dmath(xa)

    def xMax( self ) :
        """Returns the maximum x value of self if self contains data, otherwise it returns None."""

        if ( len( self.data ) > 0 ) : return self.data[-1][0]
        return None

    def xMin( self ) :
        """Returns the minimum x value of self if self contains data, otherwise it returns None."""

        if ( len( self.data ) > 0 ) : return self.data[0][0]
        return None

    def copyDataToGridWsAndXsAndYs( self ) :
        """A method to be like the W_XYs' method of the same name. Used in the vis.matplotlib.plot2d module."""

        xs = [ x for x, yzs in self.data ]
        ys = set( )
        for x, yzs in self.data :
            for y, z in yzs : ys.add( y )
        ys = sorted( ys )
        zs = []
        for y in ys :
            zs.append( [ self.getValue( x, y ) for x in xs ] )
        return( xs, ys, zs )
