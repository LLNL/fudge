# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the endl4dmath class. Also see the supporting module endl4dmathmisc.py.
"""

import os
import sys

from fudge.core.utilities import fudgeFileMisc
from LUPY import subprocessing
from fudge.vis.gnuplot import plotbase
from brownies.legacy.endl import fudgemisc
from brownies.legacy.endl import endl4dmathmisc, endl3dmathmisc, endl1dmathClasses

__metaclass__ = type

class endl4dmath :
    """
This class is designed to support operations on a list of 4-column number data that may be useful to scientists.
(A number is an integer or a float.) In this documenation the four columns are referred to as t, x, y and z.
The data is to be stored as a python list of [ t, xyz-data-for-t ], where xyz-data-for-t is a list of
[ x, yz-data-for-x ] and yz-data-for-x is a list of [ y, z ] data. For example, consider the 4 column data::

 1 1 1 1
 1 1 2 2
 1 1 3 1
 1 2 1 1
 1 2 2 0
 1 2 4 0
 3 1 1 1
 3 1 2 3
 4 2 2 0
 4 2 3 1
 4 2 4 2

For this data there are 3 different t-values.  The length of this data is 3, for the 3 different t-values.
The first t-value has the data::

 1 1 1
 1 2 2
 1 3 1
 2 1 1
 2 2 0
 2 4 0

for which there are 2 x values, and the first having the data::
 1 1
 2 2 
 3 1 

All this data is stored using python lists as ::

 [ [ 1, [ [ 1, [ [ 1, 1 ],
                 [ 2, 2 ],
                 [ 3, 1 ]
               ]
          ],
          [ 2, [ [ 1, 1 ],
                 [ 2, 0 ],
                 [ 4, 0 ]
               ]
          ]
        ]
   ]
   [ 3, [ [ 1, [ [ 1, 1 ],
                 [ 2, 3 ]
               ]
          ]
        ]
   ]
   [ 4, [ [ 2, [ [ 2, 0 ],
                 [ 3, 1 ],
                 [ 4, 2 ]
               ]
          ]
        ]
   ]
 ]

Note, the xyz data is store like the data in an endl3dmath class. That is, endl3dmath data looks like a list of [ t, 3d-data ].
However, the 3d part is not an endl3dmath class but is structured like the data in an endl3dmath class. Thus, to plot
the xyz data for the (i+1)^th t-value of an endl4dmath instance fourD, create an endl3dmath object of the data as endl3dmath( fourD.data[i][1] )
and use the endl3dmath plot method.

Members ::
        data        | A python list of numbers.
        columns     | This is always 4 as this is 1d data.
        tLabel      | Labelled used as the t-Label when plotting.
        xLabel      | Labelled used as the x-Label when plotting.
        yLabel      | Labelled used as the y-Label when plotting.
        zLabel      | Labelled used as the z-Label when plotting.
"""

    def __init__(self, data, checkDataType = False, tLabel = None, xLabel = None, yLabel = None, 
        zLabel = None, label = "unknown", interpolation = 0, template = None ) :
        """Creates an endl4dmath instance. Data must be of type list[ number, list[ number, list[ number, number ] ] ].
        The default plotting labels can be set with tLabel, xLabel, yLabel and zLabel."""

        if( interpolation != 0 ) : raise Exception( "\nError in endl4dmath.__init__: Bad interpolation value = %s" % repr(interpolation) )
        self.interpolation = interpolation
        self.data = endl4dmathmisc.get4dmathData(data, "endl4dmath.__init__", "data")
        if ( template is not None ) : endl4dmathmisc.valid4dClassType(template, "endl4dmath.__init__", "template")
        self.columns = 4
        self.tLabel = tLabel
        if ( ( template is not None ) and ( tLabel is None ) ) : self.tLabel = template.tLabel
        self.xLabel = xLabel
        if ( ( template is not None ) and ( xLabel is None ) ) : self.xLabel = template.xLabel
        self.yLabel = yLabel
        if ( ( template is not None ) and ( yLabel is None ) ) : self.yLabel = template.yLabel
        self.zLabel = zLabel
        if ( ( template is not None ) and ( zLabel is None ) ) : self.zLabel = template.zLabel
        self.label = label
        if(   self.interpolation == 0 ) :
            self.xInterpolation = 'linear,linear,linear,linear'
        else :
            raise Exception( 'Unsupported interpolation = %d for 4d data' % self.interpolation )

    def __getitem__( self, txy ) :
        """Returns self's z value evaluated at txy = (t, x, y). If t (x or y) is outside the range of self's t (x or y) domain then 0. is returned."""

        t = txy[0]
        x = txy[1]
        y = txy[2]
        l = len( self.data )
        if( l < 1 ) : return( 0. )
        z = 0.
        if( self.data[0][0] <= t <= self.data[-1][0] ) :
            i = 0
            while( i < l ) :
                if( t <= self.data[i][0] ) : break
                i += 1
            if( t == self.data[i][0] ) :
                z = endl3dmathmisc.interpolate_XY(self.data[i][1], x, y)
            else :
                d1 = self.data[i-1]
                d2 = self.data[i]
                z1 = endl3dmathmisc.interpolate_XY(d1[1], x, y)
                z2 = endl3dmathmisc.interpolate_XY(d2[1], x, y)
                t1 = d1[0]
                t2 = d2[0]
                z = ( z1 * ( t2 - t ) + z2 * ( t - t1 ) ) / ( t2 - t1 )
        return( z )

    def __len__( self ) :                                   # class endl4dmath
        """Returns the number of items in the top list of self (i.e., the number of t-values)."""

        return len( self.data )

    def __repr__( self ) :                                  # class endl4dmath
        """Returns a printable string of the data in self. This method uses endl4dmathmisc.endl4d_repr_tFormat, 
        endl4dmathmisc.endl4d_repr_xFormat, endl4dmathmisc.endl4d_repr_yFormat and endl4dmathmisc.endl4d_repr_zFormat
        to convert each point to a string."""

        strData = endl4dmathmisc.string4dData(self.data)
        return( '\n'.join( strData ) + '\n' )

    def copyData( self ) :
        """Returns an endl4dmath instance that is a copy, and not a reference, of self."""

        d4 = []
        for t_xyz in self.data :
            d3 = []
            for x_yz in t_xyz[1] :
                d2 = []
                for yz in x_yz[1] : d2.append( [ yz[0], yz[1] ] )
                d3.append( [ x_yz[0], d2 ] )
            d4.append( [ t_xyz[0], d3 ] )
        return endl4dmath( d4, checkDataType = 0, template = self )

    def getDimensions( self ) :
        """Returns the dimensions (4 for endl4dmath) for this type of data."""

        return( 4 )

    @property
    def dimension( self ) :

        return( self.getDimensions( ) )

    def toString( self, format = None ) :
        """Returns the string returned by the endl4dmath's __repr__ function. This can be useful when endl4dmath is used as a base class."""

        endl4d_repr_tFormat = endl4dmathmisc.endl4d_repr_tFormat
        endl4d_repr_xFormat = endl4dmathmisc.endl4d_repr_xFormat
        endl4d_repr_yFormat = endl4dmathmisc.endl4d_repr_yFormat
        endl4d_repr_zFormat = endl4dmathmisc.endl4d_repr_zFormat
        if( format is not None ) :
            endl4dmathmisc.endl4d_repr_tFormat = format
            endl4dmathmisc.endl4d_repr_xFormat = format
            endl4dmathmisc.endl4d_repr_yFormat = format
            endl4dmathmisc.endl4d_repr_zFormat = format
        s = endl4dmath.__repr__( self )
        endl4dmathmisc.endl4d_repr_tFormat = endl4d_repr_tFormat
        endl4dmathmisc.endl4d_repr_xFormat = endl4d_repr_xFormat
        endl4dmathmisc.endl4d_repr_yFormat = endl4d_repr_yFormat
        endl4dmathmisc.endl4d_repr_zFormat = endl4d_repr_zFormat
        return( s )

    def toStringWithPrefixSuffix( self, Prefix = "", Suffix = "" ) :
        """Returns a printable string of the data in self with Prefix and Suffix append to each line.
        See __repr__ to change the output format."""

        ssyz = "%%s %%s %s %s" % (endl4dmathmisc.endl4d_repr_yFormat, endl4dmathmisc.endl4d_repr_zFormat)
        stxyz = []
        for t_etal in self.data :
            t = endl4dmathmisc.endl4d_repr_tFormat % t_etal[0]
            for x_etal in t_etal[1] :
                x = endl4dmathmisc.endl4d_repr_xFormat % x_etal[0]
                for etal in x_etal[1] : stxyz.append( ssyz % ( t, x, etal[0], etal[1] ) )
        return fudgemisc.stringWithPrefixSuffix( stxyz, Prefix, Suffix )

    def ls( self ) :                                        # class endl4dmath
        """Prints the number of items in the top list of self and then the number of items in each item of the top list."""

        print( "list of %d lists of len (" % len(self.data) )
        c = ""
        s = ""
        for i in self.data :
            s += "%s%d" % ( c, len( i[1] ) )
            c = ", "
        print( s, ")" )

    def ll( self ) :                                        # class endl4dmath
        """Same as __repr__( )."""

        print( endl4dmath.__repr__( self ) )

    def plot( self, xMin = None, xMax = None, yMin = None, yMax = None, zMin = None, zMax = None, xyzlog = 0, \
        tLabel = None, xLabel = None, yLabel = None, zLabel = None, title = None, tScaleLabel = None, \
        xrot = None, zrot = None, style = None ) :                        # class endl4dmath
        """
        Plots the data.

        xyzlog values and meaning are ::

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

        if ( tLabel is None ) and ( self.tLabel is not None ) : tLabel = self.tLabel
        if ( xLabel is None ) and ( self.xLabel is not None ) : xLabel = self.xLabel
        if ( yLabel is None ) and ( self.yLabel is not None ) : yLabel = self.yLabel
        if ( zLabel is None ) and ( self.zLabel is not None ) : zLabel = self.zLabel

        dt = plotbase.parsePlotOptions( xMin, xMax, yMin, yMax, xLabel, yLabel, title, zMin = zMin, zMax = zMax, \
            zLabel = zLabel, tLabel = tLabel, tScaleLabel = tScaleLabel, xrot = xrot, zrot = zrot, style = style )
        f = fudgeFileMisc.FudgeTempFile( )
        format = fudgemisc.getFormat( self )
        f.write( endl4dmath.toString( self, format = format ) )
        f.close( )
        p = os.path.join( __file__.split( '/fudge/legacy/' )[0], "fudge", "vis", "gnuplot", "endl4dplot.py" )
        python = sys.executable
        s = [ python, p, 'xyzlog', str( xyzlog ) ] + dt + [ f.getName( ) ]
        subprocessing.spawn( s )

    def tArray( self ) :
        """Returns an endl1dmath instances of the t values in self."""

        ta = []
        for t, xyz in self.data : ta.append( t )
        return endl1dmathClasses.endl1dmath(ta)

    def tMax( self ) :
        """Returns the maximum t value of self if self contains data, otherwise it returns None."""

        if ( len( self.data ) > 0 ) : return self.data[-1][0]
        return None

    def tMin( self ) :
        """Returns the minimum t value of self if self contains data, otherwise it returns None."""

        if ( len( self.data ) > 0 ) : return self.data[0][0]
        return None
