# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.core.utilities import brb
import XYs, axes, LegendreSeries
from fudge.core import ancestry

__metaclass__ = type

monikerRegions = 'regions'
monikerRegion = 'region'

class regions( ancestry.ancestry ) :

    def __init__( self, regionClass, axes_, parent = None, isPrimaryXData = True ) :

        ancestry.ancestry.__init__( self, monikerRegions, parent )
        self.regionClass = regionClass
        if( not( axes.isValidAxes( axes_ ) ) ) : raise Exception( 'axes_ argument is not instance of axes class: it is %s' % brb.getType( axes_ ) )
        self.axes = axes_.copy( parent = self )
        self.isPrimaryXData = isPrimaryXData
        self.regions = []

    def __len__( self ) :
        """Returns the number of regions in self."""

        return( len( self.regions ) )

    def __getitem__( self, i ) :
        """Returns the (i-1)^th region of self."""

        return( self.regions[i] )

    def __setitem__( self, i, region_ ) :
        """Set region (the right-hand-side) to the :math:`{i-1}^{th}` region of self. Region must be an instance of self's regionClass.
        If :math:`i > 0`, the following must also be met:
        
            - all the prior regions must already exists,
            - region's minimum domain value must be equal to the prior region's maximum domain value.
        """

        if( not( isinstance( region_, self.regionClass ) ) ) : raise Exception( 'Region must be an instance of %s. It is %s' % \
            ( self.regionClass.__class__, brb.getType( region_ ) ) )
        n = len( self )
        if( not( 0 <= i <= n ) ) : raise IndexError( 'Index = %s not in range 0 <= index <= %d' % ( i, n ) )
        XYs.raiseNotSameUnit( self.axes[0], region_.axes[0] )
        XYs.raiseNotSameUnit( self.axes[1], region_.axes[1] )
        if( len( region_ ) < 2 ) : raise Exception( 'Region must contain at least 2 points; it has %s.' % len( region_ ) )
        if( isinstance( region_, ( XYs.XYs, LegendreSeries.W_XYs_LegendreSeries ) ) ) :
            region = region_.copy( parent = self, moniker = monikerRegion, index = i, isPrimaryXData = False )
        else :
            region = region_.copy( parent = self, moniker = monikerRegion, index = i )
        if( len( self ) == 0 ) :
            self.regions.append( region )
        else :
            if( i > 0 ) :
                if( self.regions[i-1].domainMax( ) != region.domainMin( ) ) : raise Exception( "Prior region's domainMax %s != new region's domainMin = %s" \
                    % ( self.regions[i-1].domainMax( ), region.domainMin( ) ) )
            if( i < ( n - 1 ) ) :
                if( self.regions[i+1].domainMin( ) != region.domainMax( ) ) : raise Exception(  "Next region's domainMin %s != new region's domainMax = %s" \
                    % ( self.regions[i-1].domainMin( ), region.domainMax( ) ) )
            if( i == n ) :
                self.regions.append( region )
            else :
                self.regions[i] = region

    def append( self, region ) :

            self[len( self )] = region

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.regions[0].domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.regions[-1].domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomainUnit( self ) :

        return( self.axes[0].getUnit( ) )

    def toXMLList( self, indent = "" ) :

        tag = self.moniker
        if( hasattr( self, 'tag' ) ) : tag = self.tag       # See note 1
        xDataString = ''
        if( self.isPrimaryXData and ( hasattr( self, 'xData' ) ) ) : xDataString = ' xData="%s"' % self.xData   # See note 1
        indent2 = indent + '  '
        xmlString = [ '%s<%s%s>' % ( indent, tag, xDataString ) ]
        xmlString += self.axes.toXMLList( indent2 )
        for region in self.regions : xmlString += region.toXMLList( indent = indent2, pairsPerLine = 100 )
        xmlString[-1] += '</%s>' % tag
        return( xmlString )

class regionsXYs( regions ) :

    xData = 'regionsXYs'

    def __init__( self, axes_, parent = None, isPrimaryXData = True ) :

        regions.__init__( self, XYs.XYs, axes_, parent = parent, isPrimaryXData = isPrimaryXData )

    def toPointwiseLinear( self, accuracy, lowerEps, upperEps, removeOverAdjustedPoints = False, axes_ = None ) :
        """
        Converts the regions of self into a single ``XYs.XYs`` instance that has 'linear,linear' interpolation. At the
        boundary between two abutting regions, the x-values are the same, which is not allowed for an ``XYs.XYs`` instance.
        The arguments ``lowerEps``, ``upperEps`` are used to smear the x-values at a boundary as follows. Let :math:`(x_l, y_l)` and
        :math:`(x_u, y_u)` be the abutting points for two abutting regions. If :math:`y_l = y_u` then the point :math:`(x_u, y_u)` is removed.
        Otherwise, if( lowerEps > 0 ) the point :math:`(x_l, y_l)` is moved to :math:`x = x_l * ( 1 - lowerEps )` (or :math:`x = x_l * ( 1 + lowerEps )`
        if :math:`x_l < 0`) and the :math:`y` value is interpolated at :math:`x`. If :math:`x` is less than the x-value of the point below :math:`(x_l, y_l)`
        and ``removeOverAdjustedPoints`` is True then the point :math:`(x_l, y_l)` is removed; otherwise, a raise is executed. Similarly 
        for upperEps and the point :math:`(x_u, y_u)`.
        """

        def getAdjustedX( x, eps ) :

            if( x == 0. ) :
                x_ = eps
            elif( x < 0 ) :
                x_ = x * ( 1. - eps )
            else :
                x_ = x * ( 1. + eps )
            return( x_ )

        if( len( self.regions ) == 0 ) :
            if( axes_ is None ) : axes_ = axes.defaultAxes( )
            return( XYs( axes_, [], accuracy = accuracy ) )

        xUnit, yUnit = self.regions[0].axes[0].getUnit( ), self.regions[0].axes[1].getUnit( )
        if( lowerEps < 0. ) : raise Exception( 'lowerEps = %s must >= 0.' % lowerEps )
        if( upperEps < 0. ) : raise Exception( 'upperEps = %s must >= 0.' % upperEps )
        if( ( lowerEps == 0. ) and ( upperEps == 0. ) ) : raise Exception( 'lowerEps and upperEps cannot both be 0.' )
        xys = []
        for region in self.regions :
            if( region.getAccuracy( ) > accuracy ) : accuracy = region.getAccuracy( )
        interpolationStr = '%s,%s' % ( axes.linearToken, axes.linearToken )
        for iRegion, region in enumerate( self.regions ) :
            region_ = region.changeInterpolation( axes.linearToken, axes.linearToken, accuracy = accuracy, lowerEps = 2 * lowerEps, upperEps = 2 * upperEps )
            region_ = region_.copyDataToXYs( )
            if( iRegion > 0 ) : 
                x12, y12 = xys[-1]
                x21, y21 = region_[0]
                if( y12 == y21 ) :              # Remove first point of region as it is the same as the last point.
                    region_ = region_[1:]
                else :
                    if( lowerEps != 0. ) :
                        x11, y11 = xys[-2]
                        x = getAdjustedX( x12, -lowerEps )
                        if( x <= x11 ) :
                            if( removeOverAdjustedPoints ) :
                                del xys[-1]
                            else :
                                raise Exception( 'Adjustment at %s makes new x = %s >= prior x = %s; eps = %s' % ( x12, x, x11, lowerEps ) )
                        else :
                            xys[-1] = [ x, XYs.pointwiseXY_C.interpolatePoint( interpolationStr, x, x11, y11, x12, y12 ) ]
                    if( upperEps != 0. ) :
                        x22, y22 = region_[1]
                        x = getAdjustedX( x21, upperEps )
                        if( x >= x22 ) :
                            if( removeOverAdjustedPoints ) :
                                del region_[0]
                            else :
                                raise Exception( 'Adjustment at %s makes new x = %s >= next x = %s; eps = %s' % ( x21, x, x22, upperEps ) )
                        else :
                            region_[0] = [ x, XYs.pointwiseXY_C.interpolatePoint( interpolationStr, x, x21, y21, x22, y22 ) ]
            xys += region_
        newAxes = self.axes.copy( standAlone = True )
        newAxes[0].setInterpolation( axes.interpolationXY( axes.linearToken, axes.linearToken ) )
        pointwise = XYs.XYs( newAxes, xys, accuracy = accuracy )
        return( pointwise )

    @classmethod
    def parseXMLNode( cls, regionsElement, linkData={} ):
        attrs = dict( regionsElement.items() )
        axes_ = axes.parseXMLNode( regionsElement[0] )
        regions_ = cls( axes_ )
        for region in regionsElement[1:]:
            xys = XYs.XYs.parseXMLNode( region, linkData )
            xys.parent = regions_; xys.tag = region.tag; xys.isPrimaryXData=False
            regions_.append( xys )
        return regions_
