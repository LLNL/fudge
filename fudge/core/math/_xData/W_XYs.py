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
This module contains the class W_XYs. This class treats a list of (w_i, xys) pairs as if they were the function xys(w)
where xys is an instance of class XYs.

    axes            Description of the w, x and y data (e.g., label, interpolation, units).
"""

import axes
import XYs
from xData.ancestry import ancestry
from pqu import PQU
from fudge.core.math import fudgemath
from fudge.core.utilities import brb, subprocessing

__metaclass__ = type

noExtrapolationToken = 'noExtrapolation'
flatExtrapolationToken = 'flatExtrapolation'
monikerW_XYs = 'W_XYs'

class W_XYs( ancestry ) :

    type = monikerW_XYs
    xData = monikerW_XYs
    moniker = monikerW_XYs

    def __init__( self, axes_, data = None, index = None, value = None, parent = None, isPrimaryXData = None ) :
        """
        If the data argument is not None, it must be the argument for either the method setFromList or setFromW_XYs.
        """

        attribute = None
        if( parent is not None ) :
            attribute = 'index'
            if( index is None ) : raise Exception( 'index argument cannot be None when parent is not None' )
        ancestry.__init__( self )

        axes.isValidAxes( axes_ )
        axes_.checkDimension( 3 )

        self.axes = axes_.copy( )           # Even secondary must have axes
        self.index = index
        if( value is not None ) : value = fudgemath.toFloat( value )
        self.value = value
        if( isPrimaryXData is not None ) : self.isPrimaryXData = isPrimaryXData

        self.xys = []
        if( isinstance( data, ( list, tuple ) ) ) :
            self.setFromList( data )
        elif( isinstance( data, W_XYs ) ) :
            self.setFromW_XYs( data )

    def __len__( self ) :

        return( len( self.xys ) )

    def __getitem__( self, index ) :

        return( self.xys[index] )

    def __setitem__( self, index, xys ) :

        index_, xys_ = self._set_insertCommon( index, xys, xys.value )
        if( index_ is not None ) :
            if( index_ > 0 ) :
                if( xys_.value <= self.xys[index_-1].value ) : 
                    raise Exception( 'xys.value = %s is <= prior xys.value = %s' % ( xys_.value, self.xys[index_-1].value ) )
            if( index_ < ( len( self ) - 1 ) ) :
                if( xys_.value >= self.xys[index_+1].value ) :
                    raise Exception( 'xys.value = %s is >= next xys.value = %s' % ( xys_.value, self.xys[index_+1].value ) )
            self.xys[index_] = xys_

    def append( self, xys_ ) :

        self[len( self )] = xys_

    def insert( self, index, xys, w = None ) :
        """
        Inserts xys at index. If w is None, w is take from the value of xys.
        """

        if( w is None ) : w = xys.value
        index_, xys_ = self._set_insertCommon( index, xys, w )
        if( index_ is not None ) :
            if( index_ > 0 ) :
                if( w <= self.xys[index_-1].value ) : raise Exception( 'w = %s is <= prior xys.value = %s' % ( w, self.xys[index_-1].value ) )
            if( w >= self.xys[index_].value ) : raise Exception( 'w = %s is >= next xys.value = %s. index = %d' % ( w, self.xys[index_].value, index_ ) )
            self.xys.insert( index_, xys_ )
            for i1 in xrange( index_ + 1, len( self ) ) : self.xys[i1].index = i1

    def insertAtW( self, w, xys ) :
        """
        Inserts xys at the appropriate index for w. The inserted XYs instance will have value w, even if xys as a value.
        """

        index = -1               # Set in case self is empty and next line does not set index or xys.
        for index, xys_ in enumerate( self ) :
            if( xys_.value >= w ) : break
        if( index == -1 ) :
            index = 0
        else :
            if( xys_.value == w ) : del self.xys[index]
            if( xys_.value < w ) : index += 1               # Happends when w is greater than last items value.
        self.insert( index, xys, w = w )

    def _set_insertCommon( self, index, xys, w ) :
        """For internal use only."""

        if( not( isinstance( xys, XYs.XYs ) ) ) : raise Exception( 'right-hand-side must be instance of XYs; it is %s' % brb.getType( xys_ ) )
        n = len( self )
        if( n < index ) : raise IndexError( 'index = %s while length is %s' % ( index, n ) )
        index_ = index
        if( index_ < 0 ) : index_ += n
        if( index_ < 0 ) : raise IndexError( 'index = %s' % index )
        xys_ = xys.copy( parent = self, index = index_, value = w, axes_ = axes.referenceAxes( self ), isPrimaryXData = False )
        if( n == 0 ) :
            self.xys.append( xys_ )
            return( None, None )
        elif( n == index_ ) :
            if( w <= self.xys[-1].value ) : raise Exception( 'w = %s is <= prior xys.value = %s' % ( w, self.xys[-1].value ) )
            self.xys.append( xys_ )
            return( None, None )
        return( ( index_, xys_ ) )

    def copy( self, index = None, value = None, parent = None, axes = None, isPrimaryXData = None ) :

        if( index is None ) : index = self.index
        if( value is None ) : value = self.value
        if( axes is None ) : axes = self.axes
        n = W_XYs( axes, index = index, value = value, parent = parent, isPrimaryXData = isPrimaryXData )
        for i, xys in enumerate( self ) : n[i] = xys
        return( n )

    def copyDataToW_XYs( self, wUnit = None, xUnit = None, yUnit = None ) :

        wScale = 1.0
        if( wUnit is not None ) : wScale = PQU.PQU( '1 ' + self.axes[0].getUnit( ) ).getValueAs( wUnit )
        w_xys = []
        for xys in self : w_xys.append( [ wScale * xys.value, xys.copyDataToXYs( xUnit = xUnit, yUnit = yUnit ) ] )
        return( w_xys )

    def copyDataToGridWsAndXsAndYs( self, wUnit = None, xUnit = None, yUnit = None ) :

        from fudge.core.utilities.brb import uniquify
        wScale = 1.0
        if( wUnit is not None ) : wScale = PQU.PQU( '1 ' + self.axes[0].getUnit( ) ).getValueAs( wUnit )
        w_xys = self.copyDataToW_XYs( wUnit, xUnit, yUnit )
        ws = []
        xs = []
        ys = []
        for w_xy in w_xys:
            ws.append( w_xy[0] )
            xs += [ xy[0] for xy in w_xy[1] ]
        xs = sorted( uniquify( xs ) )
        for x in xs:
            ys.append( [ self.getValue( w, x ) for w in ws ] )
        return [ ws, xs, ys ]

    def getDimensions( self ) :
        """Returns the dimensions (3 for W_XYs) for this type of data."""

        return( 3 )

    def getValue( self, w, x = None ) :
        """
        Returns an XYs instance of self evaluated at w if x is None. Otherwise, returns, as a float, the value of self at w, x.
        """

        XYs_atW = self.interpolateAtW( w )
        if( x is None ) : return( XYs_atW )
        return( XYs_atW.getValue( x ) )

    def integrate( self, wMin = None, wMax = None, xMin = None, xMax = None ) :

        if( len( self ) == 0. ) : return( 0. )
        wis = [ [ xys.value, xys.integrate( xMin = xMin, xMax = xMax ) ] for xys in self ]
        axes_ = axes.defaultAxes( dimension = 2, independentInterpolation = axes.linearToken, dependentInterpolation = axes.linearToken )
        return( XYs.XYs( axes_, wis, self[0].getAccuracy( ) ).integrate( wMin, wMax ) )

    def interpolateAtW( self, w, unitBase = False, extrapolation = noExtrapolationToken ) :
        """Returns the interpolation of two XYs, xys1 and xys2, at w as a list of [x,y] pairs where wl = xys1.value <= w <= xys2.value = wu.
        If w is outside the domain of wl to wu, are raise is executed. Linear interpolation is performed on the xy data.  If unitBase is
        True, then unit base interpolation is performed."""

        if( len( self ) == 0 ) : raise Exception( "No data to interpolate" )
        if( w < self[0].value ) :
            if( extrapolation == flatExtrapolationToken ) :
                return( self[0].copy( value = w ) )
            else :
                raise Exception( "Interpolation point = %s less than %s" % ( w, self[0].value ) )
        if( w > self[-1].value ) :
            if( extrapolation == flatExtrapolationToken ) :
                return( self[-1].copy( value = w ) )
            else :
                raise Exception( "Interpolation point = %s greater than %s" % ( w, self[-1].value ) )
        for index, xy2 in enumerate( self ) :
            if( xy2.value >= w ) : break
        if( w == xy2.value ) : return( xy2.copy( value = w ) )
        xy1 = self[index-1]
        if( unitBase ) :
            xy = XYs.pointwiseXY_C.unitbaseInterpolate( w, xy1.value, xy1, xy2.value, xy2 )
        else :
            f = ( xy2.value - w ) / ( xy2.value - xy1.value )
            xy = f * xy1 + ( 1. - f ) * xy2
        xyp = xy1.copy( value = w )
        xyp.setData( xy )
        return( xyp )

    def normalize( self, insitu = True ) :

        n = self
        if( not( insitu ) ) : n = self.copy( )
        for xys in n.xys : xys.normalize( insitu = True )
        return( n )

    def setFromList( self, w_xys, accuracy = 1e-3 ) :
        """
        Replaces self's data with the data in w_xys. w_xys must be a list of [ w_i, xys_i ] pairs with w_i < w_{i+1}. xys_i can be
        either a list of [ x_j, y_j ] pairs with x_j < x_{j+1} or an XYs instance. For example:
            w_xys = [ [ 1e-6, [ [ -1, 0.5 ], [ 1, 0.5 ] ] ], [ 20, [ [ -1, 0. ], [ 0, 0. ], [ 1, 2. ] ] ] ].
        """

        self.xys = []
        for w1, xys in w_xys :
            if( isinstance( xys, XYs.XYs ) ) :
                xys_ = xys.copy( value = w1 )
            else :
                xys_ = XYs.XYs( axes.referenceAxes( self ), xys, accuracy, value = w1 )
            self.append( xys_ )

    def setFromW_XYs( self, w_xys ) :
        """
        Replaces self's data with the data in w_xys. w_xys must be a W_XYs instance.
        """

        if( not( isinstance( w_xys, W_XYs ) ) ) : raise TypeError( 'instance of of type W_XYs. Is type "%d"' % brb.getType( w_xys ) )
        w_xys_ = [ ( xys.value, xys ) for xys in w_xys ]
        self.setFromList( w_xys_ )

    def thin( self, accuracy ) :
        """Remove extra points in both the W and X-directions, while maintaining specified accuracy.
        Returns a thinned copy (original is unchanged)."""
        def interpolate_wxys( interp, x, x1, y1, x2, y2 ):
            if interp=='lin-lin':
                return ( y1 * ( x2 - x ) + y2 * ( x - x1 ) ) / ( x2 - x1 )

        def thin2( thinned, kThin, accuracy, ilow, ihigh ):
            """Recursively bisect w_xys, keep track of which xys can be safely removed
            using binary array 'kThin'."""
            dYMax, dYRMax, iMax = 0, 0, 0
            y1 = thinned[ilow]; y2 = thinned[ihigh]
            x1 = y1.value; x2 = y2.value
            if ilow+1 >= ihigh: return
            for i in range(ilow+1,ihigh):
                y = interpolate_wxys( 'lin-lin', thinned[i].value, x1, y1, x2, y2 )
                s = 0.5 * ( abs(y) + abs(thinned[i]) )
                dY = abs( y - thinned[i] )
                dYR = 0
                dYR = dY / s    # need safeDivide for this step
                # dY and dYR are XYs objects. How does the maximum compare to accuracy?
                dY = dY.yMax(); dYR = dYR.yMax()
                if ( dYR > dYRMax ) or ( ( dYR >= 0.9999 * dYRMax ) and ( dY > dYMax ) ):
                    iMax = i;
                    if dY > dYMax: dYMax = dY
                    if dYR > dYRMax: dYRMax = dYR
            if dYRMax < accuracy:
                for i in range(ilow+1,ihigh): kThin[i] = True
            else:
                thin2( thinned, kThin, accuracy, ilow, iMax )
                thin2( thinned, kThin, accuracy, iMax, ihigh )

        try:
            thinned = self.copy()
            for i in range(len(thinned)):   # first thin each XYs obj
                xys = thinned[i].thin( accuracy )
                xys.value = thinned[i].value
                thinned[i] = xys

            # may wish to remove middle xys object if surrounding objects are all equal. Leave for later, though
            if thinned.axes[0].interpolation.independent==axes.linearToken:
                length = len(thinned)
                kThin = [False] * len(thinned)  # boolean array
                thin2( thinned, kThin, accuracy, 0, length-1 )  # find which W-values can be removed

                # need one more copy
                rThinned = self.__class__( self.axes )
                for i,thin in enumerate(kThin):
                    if not thin:
                        rThinned.append( thinned[i] )
                return rThinned #, kThin
        except Exception as e:
            print "Exception encountered during thinning:",e

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQU.PQU( '1 ' + self.domainUnit( ) ).getValueAs( unitTo ) )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( PQU.valueOrPQ( self.xys[0].value, unitFrom = self.axes[0].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( PQU.valueOrPQ( self.xys[-1].value, unitFrom = self.axes[0].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainGrid( self, unitTo = None ) :

        scale = self.domainUnitConversionFactor( unitTo )
        return( [ scale * xys.value for xys in self ] )

    def domainUnit( self ) :

        return( self.axes[0].getUnit( ) )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0, cls = None ) :

        from fudge.core.math import miscellaneous

        independent, dependent, qualifier = self.axes[0].interpolation.getInterpolationTokens( )
        if( independent not in [ axes.linearToken ] ) : raise Exception( 'independent interpolation = %s not supported' % dependent )
        if( dependent not in [ axes.linearToken, axes.flatToken ] ) : raise Exception( 'dependent interpolation = %s not supported' % dependent )
        axes_ = self.axes.copy( )
        axes_[0].setInterpolation( axes.interpolationXY( axes.linearToken, axes.linearToken, qualifier ) )
        axes_[1].setInterpolation( axes.interpolationXY( axes.linearToken, axes.linearToken ) )
        if( cls is None ) : cls = W_XYs
        if( hasattr( self, 'isPrimaryXData' ) and self.isPrimaryXData ) :
            n = cls( axes_, self.productFrame )
        else:
            n = cls( axes_ )
        w_prior = None
        independentXY, dependentXY, qualifier = self.axes[1].interpolation.getInterpolationTokens( )
        for w_xys in self :
            if( ( independentXY != axes.linearToken ) or ( dependentXY != axes.linearToken ) ) :
                value = w_xys.value
                w_xys = w_xys.toPointwise_withLinearXYs( accuracy, lowerEps = lowerEps, upperEps = upperEps )
                w_xys.value = value
            if( w_prior is not None ) :
                if( independent == axes.flatToken ) :
                    value = miscellaneous.shiftFloatDownABit( w_xys.value, lowerEps )
                    n.append( w_prior.copy( value = value ) )
            n.append( w_xys )
            w_prior = w_xys
        n.value = self.value
        return( n )

    def toXML( self, tag = 'xData', indent = '', incrementalIndent = '  ', pairsPerLine = 100, xyFormatter = None, xySeparater = ' ' ) :

        return( '\n'.join( self.toXMLList( tag = tag, indent = indent, incrementalIndent = incrementalIndent, pairsPerLine = pairsPerLine, xyFormatter = xyFormatter, 
            xySeparater = xySeparater ) ) )

    def toXMLList( self, tag = 'xData', indent = '', incrementalIndent = '  ', pairsPerLine = 100, xyFormatter = None, xySeparater = ' ' ) :

        if( hasattr( self, 'tag' ) ) : tag = self.tag
        if( xyFormatter is None ) : xyFormatter = XYFormatter
        indent2, indexValueStr  = indent + incrementalIndent, ''
        if( self.value is not None ) : indexValueStr = ' value="%s"' % self.value
        if( self.index is not None ) : indexValueStr += ' index="%s"' % self.index
        xDataString = ''
        if( hasattr( self, 'isPrimaryXData' ) ) :
            if( self.isPrimaryXData ) : xDataString = ' xData="%s"' % self.xData
        else :
            xDataString = ' type="%s"' % self.type
        extraXMLAttributeString = ''
        if( hasattr( self, 'extraXMLAttributeString' ) ) : extraXMLAttributeString = ' ' + self.extraXMLAttributeString( )
        xmlString = [ '%s<%s%s%s%s>' % ( indent, tag, indexValueStr, xDataString, extraXMLAttributeString ) ] 
        xmlString += self.axes.toXMLList( indent = indent2 )
        for xys in self.xys : xmlString += xys.toXMLList( tag = self.axes[0].getLabel( ), indent = indent2, incrementalIndent = incrementalIndent, \
            pairsPerLine = pairsPerLine, xyFormatter = xyFormatter, xySeparater = xySeparater, oneLine = True )
        xmlString[-1] += '</%s>' % tag
        return( xmlString )

    def toString( self, prefix = '' ) :

        lines = []
        for xys in self :
            for xyStr in xys.toString( ).split( '\n' ) :
                wStr = "%s%16.8e %s" % ( prefix, xys.value, xyStr )
                lines.append( wStr )
            del lines[-1]
        return( '\n'.join( lines ) )

    def plot( self, xMin = None, xMax = None, yMin = None , yMax = None, zMin = None , zMax = None, xyzlog = 0, title = '' ) :

        from fudge.core import fudgemisc
        from fudge.core.utilities import fudgeFileMisc
        from fudge.vis.gnuplot import plotbase
        import os

        def getUnitlessNumber( value, unit, default ) :

            if( value is None ) : value = default
            if( fudgemath.isNumber( value ) ) : return( float( value ) )
            if( type( '' ) == type( value ) ) :
                try :
                    return( float( value ) )
                except :
                    value = PQU.PQU( value )
            if( isinstance( value, PQU.PQU ) ) : return( value.getValueAs( unit ) )
            raise Exception( 'Cannot convert %s to a unitless number' % str( value ) )

        xLabel = self.axes[0].plotLabel( )
        yLabel = self.axes[1].plotLabel( )
        zLabel = self.axes[2].plotLabel( )

        xMin = getUnitlessNumber( xMin, self.axes[0].getUnit( ), self.domainMin( ) )
        xMax = getUnitlessNumber( xMax, self.axes[0].getUnit( ), self.domainMax( ) )
        for xys in self :
            yMin2 = getUnitlessNumber( yMin, self.axes[1].getUnit( ), xys.xMin( ) )
            yMin = min( yMin2, xys.xMin( ) )
            yMax2 = getUnitlessNumber( yMax, self.axes[1].getUnit( ), xys.xMax( ) )
            yMax = max( yMax2, xys.xMax( ) )
            zMin2 = getUnitlessNumber( zMin, self.axes[1].getUnit( ), xys.yMin( ) )
            zMin = min( zMin2, xys.yMin( ) )
            zMax2 = getUnitlessNumber( zMax, self.axes[1].getUnit( ), xys.yMax( ) )
            zMax = max( zMax2, xys.yMax( ) )

        dt = plotbase.parsePlotOptions( xMin, xMax, yMin, yMax, xLabel, yLabel, title, zMin = zMin, zMax = zMax, zLabel = zLabel )
        f = fudgeFileMisc.fudgeTempFile( )
        f.write( self.toString( ) )
        f.close( )
        p = os.path.join( __file__.split( '/fudge/core/' )[0], "fudge", "vis", "gnuplot", "endl3dplot.py" )
        s = [ 'python', p, 'xyzlog', str( xyzlog ) ] + dt + [ f.getName( ) ]
        subprocessing.spawn( s )

    @classmethod
    def parseXMLNode( cls, xdataElement, xPath=[], linkData={} ):
        """ translate W_XYs and its subclasses from xml """
        xPath.append( xdataElement.tag )
        attrs = dict( xdataElement.items() )
        if "xData" in attrs:
            assert attrs.pop("xData") == "W_XYs"
            attrs['isPrimaryXData'] = True
        else: attrs['isPrimaryXData'] = False
        axes_ = axes.parseXMLNode( xdataElement[0], xPath )
        wxys_ = cls( axes_, **attrs )
        wxys_.tag = xdataElement.tag
        for xys in xdataElement[1:]:
            # doesn't work yet:
            #wxys_.append( XYs.XYs.parseXMLNode( xys, xPath ) )
            xPath.append( '%s[@value="%s"]' % (xys.tag, xys.get('value') ) )
            axes_ = axes.referenceAxes( parent=wxys_ )
            data = map(float, xys.text.split())
            data = zip( data[::2], data[1::2] )
            wxys_.append( XYs.XYs( axes_, data, float(xys.get("accuracy")), index=int(xys.get("index")),
                value=float(xys.get("value")), parent=wxys_ ) )
            xPath.pop()
        xPath.pop()
        return wxys_

    @staticmethod
    def defaultAxes( labelsUnits = {} ) :

        return( axes.defaultAxes( dimension = 3, labelsUnits = labelsUnits ) )

def XYFormatter( x, y ) :

    return( '%s %s' % ( PQU.toShortestString( x ), PQU.toShortestString( y ) ) )

if( __name__ == '__main__' ) :

    vl_w_xy = axes.axes( dimension = 3 )
    vl_w_xy[0] = axes.axis( 'energy_in', 0, 'eV', interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
    vl_w_xy[1] = axes.axis( 'mu', 0, '', interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
    vl_w_xy[2] = axes.axis( 'P(energy_in|mu)', 0, '' )
    vl_w_xy.checkDimension( 3 )

    vl_xy = axes.axes( )
    vl_xy[0] = vl_w_xy[1]
    vl_xy[1] = vl_w_xy[2]

    w_xys1 = W_XYs( vl_w_xy )
    w_xys1[0] = XYs.XYs( vl_xy, [ [ 1, 0 ], [ 3, 2 ], [ 4, 1 ] ], 1e-3, biSectionMax = 7, value = 1.3 )
    w_xys1[1] = XYs.XYs( vl_xy, [ [ 1, 0 ], [ 4, 2 ] ], 1e-3, biSectionMax = 7, value = 1.4 )
    w_xys1[2] = XYs.XYs( vl_xy, [ [ 1, 1 ], [ 2, 2 ], [ 4, 3 ], [ 7, 1 ] ], 1e-3, biSectionMax = 7, value = 2.5 )
    w_xys1[3] = XYs.XYs( vl_xy, [ [ 1, 0 ], [ 2.5, 2 ], [ 3, 3 ], [ 7, 0 ] ], 1e-3, biSectionMax = 7, value = 3.2 )

    print
    print w_xys1.domainMin( ), w_xys1.domainMin( 'eV' ), w_xys1.domainMin( 'MeV' ), w_xys1.domainMin( asPQU = True ), w_xys1.domainMin( 'MeV', True )
    print w_xys1.domainMax( ), w_xys1.domainMax( 'eV' ), w_xys1.domainMax( 'MeV' ), w_xys1.domainMax( asPQU = True ), w_xys1.domainMax( 'MeV', True )
    print w_xys1.toXML( indent = '    ', pairsPerLine = 5 )

    w_xys2 = w_xys1.copy( )
    w_xys2.append( XYs.XYs( vl_xy, [ [ 1, 0 ], [ 2, 3 ], [ 2.5, 2 ], [ 6, 1 ] ], 1e-3, biSectionMax = 7, value = 3.21 ) )

    print
    print w_xys1.toXML( indent = '    ', pairsPerLine = 5 )
    print
    print w_xys2.toXML( indent = '    ', pairsPerLine = 5 )

    print
    w_xys2.normalize( )
    print w_xys2.toXML( indent = '    ', pairsPerLine = 5 )
