# <<BEGIN-copyright>>
# <<END-copyright>>

import math
from fudge.core.utilities import brb
from fudge.core.ancestry import ancestry
from fudge.core.math import fudgemath
import axes, XYs, W_XYs
from pqu import physicalQuantityWithUncertainty

__metaclass__ = type

monikerXYs_LegendreSeries = 'XYs_LegendreSeries'
monikerW_XYs_LegendreSeries = 'W_XYs_LegendreSeries'
monikerV_W_XYs_LegendreSeries = 'V_W_XYs_LegendreSeries'

noExtrapolationToken = 'noExtrapolation'
flatExtrapolationToken = 'flatExtrapolation'

class XYs_LegendreSeries( ancestry ) :
    """This class stores and manipulates an angular pdf (i.e. pdf(mu) where mu is the cos of the angle) 
represented as Legendre coefficients. The pdf and Legendre coefficients are related by

    pdf(mu) = Sum_over_l_of ( l + 0.5 ) * C_l * P_l(mu)

where the sum is from l = 0 to lMax, lMax is the highest Legendre coefficients in the instance, C_l 
is the Legendre coefficient for Legendre order l and P_l(mu) is the Legendre polynomial of order l."""

    moniker = monikerXYs_LegendreSeries

    def __init__( self, unit, coefficients, index = None, value = None, parent = None ) :

        ancestry.__init__( self, self.moniker, parent )
        self.unit = unit
        self.index = index
        self.value = value
        self.coefficients = []
        for c_l in coefficients : self.coefficients.append( float( c_l ) )

    def __len__( self ) :
        """Returns the number of Legendre coefficients in the instance (i.e., lMax)."""

        return( len( self.coefficients ) )

    def __getitem__( self, l ) :
        """Returns the (l+1)^th Legendre coefficient."""

        return( self.coefficients[l] )

    def __setitem__( self, l, c_l ) :
        """Sets the (l+1)^th Legendre coefficient to c_l. l must be less than or equal to lMax."""

        if( l == len( self ) ) :
            self.coefficients.append( float( c_l ) )
        else :
            self.coefficients[l] = float( c_l ) 

    def __str__( self ) :
        """Returns a string representation of the Legendre coefficients of self."""

        return( ' '.join( [ "%g" % c_l for c_l in self.coefficients ] ) )

    def copy( self, index = None, value = None, parent = None ) :
        """Creates a new XYs_LegendreSeries that is a copy of self except for its parent. The new 
        instance's index and value members are changes if index or value arguments are not None 
        respectively. The new instance's parent determined by the parent argument."""

        if( index is None ) : index = self.index
        if( value is None ) : value = self.value
        n = XYs_LegendreSeries( self.unit, self.coefficients, index = index, value = value, parent = parent )
        return( n )

    def getValue( self, mu ) :
        """Using the Legendre coefficients, this method calculates pdf(mu) and returns it."""

        P = 0.
        for l, c_l in enumerate( self.coefficients ) : P += ( l + 0.5 ) * c_l * fudgemath.Legendre( l, mu, checkXRange = False ) 
        return( P )

    def invert( self ) :
        """This method returns a XYs_LegendreSeries instance that is the mirror in of self about mu = 0.
        That is, returns a Legendre series which represent self's pdf(-mu)."""

        c = self.copy( )
        sign = 1
        for l in xrange( len( c ) ) :
            c.coefficients[l] *= sign
            sign *= -1
        return( c )

    def isIsotropic( self ) :
        """Returns True if self is isotropic."""

        for coefficient in self[1:] :
            if( coefficient != 0. ) : return( False )
        return( True )

    def toPointwiseLinear( self, accuracy ) :
        """The method constructs the pdf(mu) versus mu and returns it as a XYs instance. The accuracy of the 
        reconstruction (hence the number of points in the returned XYs) is determined by the accuracy argument."""

        if( accuracy < 1e-6 ) : accuracy = 1e-6
        if( accuracy > 0.1 ) : accuracy = 1.0
        P, n = [], 100
        for i in xrange( n ) :
            mu = -1. + ( 2. * i ) / n
            P.append( [ mu, self.getValue( mu ) ] )
        P.append( [ 1., self.getValue( 1. ) ] )
        axes_ = axes.axes( )
        axes_[0] = axes.axis( 'mu', 0, '', interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )  # What about frame??????
        axes_[1] = axes.axis( 'P(mu)', 1, self.unit )
        P = XYs.XYs( axes_, P, accuracy )
        return( P.thin( accuracy = accuracy ) )

    def toXML( self, tag = 'xData', indent = '', incrementalIndent = '  ' ) :
        """This method returns the XML string representation of self."""

        return( '\n'.join( self.toXMLList( tag = tag, indent = indent, incrementalIndent = incrementalIndent ) ) )

    def toXMLList( self, tag = 'xData', indent = '', incrementalIndent = '  ', includeType = True ) :
        """This method returns self as a list of string which converts to the XML string representation of self via
        the python statement '\n'.join( xmlString )."""

        from pqu.physicalQuantityWithUncertainty import toShortestString
        indent2 = indent + incrementalIndent
        attributes = ''
        if( self.value is not None ) : attributes += ' value="%s"' % self.value
        if( self.index is not None ) : attributes += ' index="%s"' % self.index
        if( includeType ) : attributes += ' type="%s"' % self.moniker
        xmlString = [ '%s<%s%s length="%s">' % ( indent, tag, attributes, len( self ) ) ]
        xmlString += [ toShortestString(c_l) for c_l in self.coefficients ]
        xmlString[-1] += '</%s>' % tag
        return( [ ' '.join( xmlString ) ] )

class W_XYs_LegendreSeries( ancestry ) :

    moniker = monikerW_XYs_LegendreSeries
    xData = monikerW_XYs_LegendreSeries

    def __init__( self, axes_, index = None, value = None, parent = None, isPrimaryXData = True ) :

        ancestry.__init__( self, self.moniker, parent )
        self.index = index
        self.value = value
        self.axes = axes_.copy( parent = self )
        self.LegendreSeries_s = []
        self.isPrimaryXData = isPrimaryXData

    def __len__( self ) :

        return( len( self.LegendreSeries_s ) )

    def __getitem__( self, index ) :

        return( self.LegendreSeries_s[index] )

    def __setitem__( self, index, LegendreSeries_ ) :

        if( not( isinstance( LegendreSeries_, XYs_LegendreSeries ) ) ) : 
            raise Exception( 'right-hand-side must be instance of XYs_LegendreSeries; it is %s' % brb.getType( LegendreSeries_ ) )
        if( type( LegendreSeries_.value ) != type( 1. ) ) : 
            raise Exception( 'LegendreSeries value must be a float; it is a %s' % brb.getType( LegendreSeries_.value ) )
        n = len( self )
        if( n < index ) : raise IndexError( 'index = %s while length is %s' % ( index, n ) )
        if( index < 0 ) : raise IndexError( 'index = %s' % index )
        LegendreSeries_ = LegendreSeries_.copy( index = index, value = LegendreSeries_.value, parent = self )
        if( n == 0 ) :
            self.LegendreSeries_s.append( LegendreSeries_ )
        elif( n == index ) :
            if( LegendreSeries_.value <= self.LegendreSeries_s[-1].value ) : 
                raise Exception( 'LegendreSeries.value = %s is <= prior value = %s' % ( LegendreSeries_.value, self.LegendreSeries_s[-1].value ) )
            self.LegendreSeries_s.append( LegendreSeries_ )
        else :
            if( index > 0 ) :
                if( LegendreSeries_.value <= self.LegendreSeries_s[index - 1].value ) :
                    raise Exception( 'LegendreSeries.value = %s is <= prior value = %s' % ( LegendreSeries_.value, self.LegendreSeries_s[index - 1].value ) )
            if( index < ( n - 1 ) ) :
                if( LegendreSeries_.value >= self.LegendreSeries_s[index + 1].value ) :
                    raise Exception( 'LegendreSeries.value = %s is >= next value = %s' % ( LegendreSeries_.value, self.LegendreSeries_s[index + 1].value ) )
            self.LegendreSeries_s[index] = LegendreSeries_

    def append( self, LegendreSeries_ ) :

        self[len( self )] = LegendreSeries_

    def copy( self, parent = None, index = None, value = None, moniker = None, isPrimaryXData = True ) :  # Name is not used here but needed for now for regions.

        if( index is None ) : index = self.index
        if( value is None ) : value = self.value
        n = W_XYs_LegendreSeries( self.axes, index = index, value = value, parent = parent, isPrimaryXData = isPrimaryXData )
        for i, LegendreSeries_ in enumerate( self ) : n[i] = LegendreSeries_
        return( n )

    def getLegendreSeriesAtW( self, w, extrapolation = noExtrapolationToken ) :

        if( self[0].value >= w ) :
            return( self[0] )
        if( self[-1].value <= w ) :
            return( self[-1] )
        for LegendreSeries2 in self :
            if( LegendreSeries2.value >= w ) : break
            LegendreSeries1 = LegendreSeries2
        if( LegendreSeries2.value == w ) : return( LegendreSeries2 )
        independent, dependent, qualifier = self.axes[0].interpolation.getInterpolationTokens( )
        coefficients1 = LegendreSeries1.coefficients[:]
        coefficients2 = LegendreSeries2.coefficients[:]
        n = len( coefficients2 ) - len( coefficients1 )
        if( n > 0 ) :
            coefficients1 += n * [ 0. ]
        if( n < 0 ) :
            coefficients2 += -n * [ 0. ]
        if( independent == axes.linearToken ) :
            s = ( w - LegendreSeries1.value ) / ( LegendreSeries2.value - LegendreSeries1.value )
        elif( independent == axes.logToken ) :
            s = math.log( w / LegendreSeries1.value ) / math.log( LegendreSeries2.value / LegendreSeries1.value )
        else :
            raise Exception( 'Unsupported interpolation = %s' % independent  )
        if( dependent != axes.linearToken ) : raise Exception( 'Unsupported interpolation = %s' % dependent  )
        coefficients = []
        for i, Cl_1 in enumerate( coefficients1 ) :
            coefficients.append( s * ( coefficients2[i] - Cl_1 ) + Cl_1 )
        return( XYs_LegendreSeries( self.axes[-1].getUnit( ), coefficients, value = w ) )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( physicalQuantityWithUncertainty.valueOrPQ( self.LegendreSeries_s[0].value, unitFrom = self.axes[0].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( physicalQuantityWithUncertainty.valueOrPQ( self.LegendreSeries_s[-1].value, unitFrom = self.axes[0].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomainUnit( self ) :

        return( self.axes[0].getUnit( ) )

    def invert( self ) :
        """This method returns a W_XYs_LegendreSeries instance whose Legendre series are a mirror in of themselfs about mu = 0.
        That is, returns Legendre series which represent self's pdf(w,-mu)."""

        n = W_XYs_LegendreSeries( self.axes )
        for i, LegendreSeries_ in enumerate( self ) : n[i] = LegendreSeries_.invert( )
        return( n )

    def isIsotropic( self ) :

        for LegendreSeries_s in self :
            if( not( LegendreSeries_s.isIsotropic( ) ) ) : return( False )
        return( True )

    def maxLegendreOrder( self ) :

        lMax = 1
        for LegendreSeries_s in self : lMax = max( lMax, len( LegendreSeries_s ) )
        return( lMax - 1 )

    def toLegendreLinear( self, accuracy ) :

        independent, dependent, qualifier = self.axes[0].interpolation.getInterpolationTokens( )
        if( dependent not in [ axes.linearToken, axes.logToken ] ) : raise Exception( 'dependent interpolation = %s not supported' % dependent )
        if( independent not in [ axes.linearToken, axes.flatToken ] ) : raise Exception( 'independent interpolation = %s not supported' % independent )

        axes__ = self.axes.copy( standAlone = True )
        axes__[0].setInterpolation( axes.interpolationXY( axes.linearToken, axes.linearToken ) )

        pwl = self.toPointwiseLinear( accuracy )
        n = W_XYs_LegendreSeries( axes__ )
        for e_in in pwl : n.append( self.getLegendreSeriesAtW( e_in.value ) )
        return( n )

    def toPointwiseLinear( self, accuracy, axes_ = None ) :

        from fudge.core.utilities import mathMisc
        from fudge.vis.gnuplot import fudgeMultiPlots

        def logFill( n, LS1, w_xys1, LS2, w_xys2, level = 0 ) :

            def getLSMid( ls1, ls2 ) :

                lsMid = XYs_LegendreSeries( ls1.unit, ls1.coefficients )
                for l in xrange( len( ls1 ), len( ls2 ) ) : lsMid[l] = 0
                return( lsMid )

            if( level == 4 ) : return           # 4 is arbitary.
            w1, w2, wMid = LS1.value, LS2.value, 0.5 * ( w_xys1.value + w_xys2.value )
            f = math.log( w2 / wMid ) /  math.log( w2 / w1 )
            g = 1 - f
            if( len( LS1 ) < len( LS2 ) ) :
                lsMid = getLSMid( LS1, LS2 )
                for l, Cl in enumerate( lsMid ) : lsMid[l] = f * Cl + g * LS2[l]
            else :
                lsMid = getLSMid( LS2, LS1 )
                for l, Cl in enumerate( lsMid ) : lsMid[l] = f * LS1[l] + g * Cl
            w_xysMid = lsMid.toPointwiseLinear( accuracy )
            w_xysMid.value = wMid
            yMax = max( w_xys1.yMax( ), w_xysMid.yMax( ), w_xys2.yMax( ) )
            rEps = 0.1 * accuracy * yMax        # 0.1 is arbitary.
            ave = 0.5 * ( w_xys1 + w_xys2 )
            diff = ave - w_xysMid
            doMore = False
            for x, y in diff : doMore = doMore or ( ( abs( y ) > rEps ) and ( abs( y ) > accuracy * w_xysMid.getValue( x ) )  )
            if( doMore ) :
                logFill( n, LS1, w_xys1, LS2, w_xysMid, level + 1 )
                n.append( w_xysMid )
                logFill( n, LS1, w_xysMid, LS2, w_xys2, level + 1 )

        independent, dependent, qualifier = self.axes[0].interpolation.getInterpolationTokens( )
        if( independent not in [ axes.linearToken, axes.logToken ] ) : raise Exception( 'independent interpolation = %s not supported' % independent )
        if( dependent not in [ axes.linearToken, axes.logToken, axes.flatToken ] ) : raise Exception( 'dependent interpolation = %s not supported' % dependent )

        if( axes_ is None ) :
            axes__ = self.axes
            if( hasattr( axes__, '_getReferenceAxes' ) ) : axes__ = self.axes._getReferenceAxes( )
            axes_ = axes.axes( 3 )
            axes_[0] = axes__[0]
            axes_[0].setInterpolation( axes.interpolationXY( axes.linearToken, axes.linearToken, qualifier ) )
            axes_[1] = axes.axis( 'mu', 1, '', frame = axes__[1].getFrame( ), 
                interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken, qualifier ) )
            axes_[2] = axes.axis( 'P', 2, axes__[1].getUnit( ), frame = axes__[1].getFrame( ) )
        n = W_XYs.W_XYs( axes_ )
        ls_prior, w_prior = None, None
        for w_LegendreSeries in self :
            w_xys = w_LegendreSeries.toPointwiseLinear( accuracy )
            w_xys.value = w_LegendreSeries.value
            if( w_prior is not None ) :
                if( independent == axes.flatToken ) :
                    value = mathMisc.shiftFloatDownABit( w_LegendreSeries.value, accuracy )
                    n.append( w_prior.copy( value = value ) )
                elif( dependent == axes.logToken ) :
                    logFill( n, ls_prior, w_prior, w_LegendreSeries, w_xys )
            n.append( w_xys )
            ls_prior, w_prior = w_LegendreSeries, w_xys
        return( n )

    def toXML( self, tag = 'xData', indent = '', incrementalIndent = '  ' ) :

        return( '\n'.join( self.toXMLList( indent = indent, incrementalIndent = incrementalIndent ) ) )

    def toXMLList( self, tag = 'xData', indent = '', incrementalIndent = '  ' ) :

        if( hasattr( self, 'tag' ) ) : tag = self.tag
        indent2 = indent + incrementalIndent
        extraStr = ''
        if( self.index is not None ) : extraStr += ' index="%s"' % self.index
        if( self.value is not None ) : extraStr += ' value="%s"' % self.value
        xDataString = ''
        if( self.isPrimaryXData ) : xDataString = ' xData="%s"' % self.xData
        xmlString = [ '%s<%s%s%s>' % ( indent, tag, extraStr, xDataString ) ]
        xmlString += self.axes.toXMLList( indent = indent2 )
        for LegendreSeries_ in self : xmlString += LegendreSeries_.toXMLList( tag = self.axes[0].getLabel( ), indent = indent2, includeType = False )
        xmlString[-1] += '</%s>' % tag
        return( xmlString )

class V_W_XYs_LegendreSeries( ancestry ) :

    moniker = monikerV_W_XYs_LegendreSeries
    xData = monikerV_W_XYs_LegendreSeries

    def __init__( self, axes_, parent = None ) :

        ancestry.__init__( self, self.moniker, parent )
        self.axes = axes_.copy( parent = self )
        self.W_LegendreSeries_s = []

    def __len__( self ) :

        return( len( self.W_LegendreSeries_s ) )

    def __getitem__( self, index ) :

        return( self.W_LegendreSeries_s[index] )

    def __setitem__( self, index, W_LegendreSeries_ ) :

        if( not( isinstance( W_LegendreSeries_, W_XYs_LegendreSeries ) ) ) : 
            raise Exception( 'right-hand-side must be instance of W_XYs_LegendreSeries; it is %s' % brb.getType( W_LegendreSeries_ ) )
        if( type( W_LegendreSeries_.value ) != type( 1. ) ) : 
            raise Exception( 'W_LegendreSeries value must be a float; it is a %s' % brb.getType( W_LegendreSeries_.value ) )
        n = len( self )
        if( n < index ) : raise IndexError( 'index = %s while length is %s' % ( index, n ) )
        if( index < 0 ) : raise IndexError( 'index = %s' % index )
        W_LegendreSeries_ = W_LegendreSeries_.copy( index = index, value = W_LegendreSeries_.value, parent = self, isPrimaryXData = False )
        if( n == 0 ) :
            self.W_LegendreSeries_s.append( W_LegendreSeries_ )
        elif( n == index ) :
            if( W_LegendreSeries_.value <= self.W_LegendreSeries_s[-1].value ) : 
                raise Exception( 'W_LegendreSeries.value = %s is <= prior value = %s' % ( W_LegendreSeries_.value, self.W_LegendreSeries_s[-1].value ) )
            self.W_LegendreSeries_s.append( W_LegendreSeries_ )
        else :
            if( index > 0 ) :
                if( W_LegendreSeries_.value <= self.W_LegendreSeries_s[index - 1].value ) :
                    raise Exception( 'W_LegendreSeries.value = %s is <= prior value = %s' % ( W_LegendreSeries_.value, self.W_LegendreSeries_s[index - 1].value ) )
            if( index < ( n - 1 ) ) :
                if( W_LegendreSeries_.value >= self.W_LegendreSeries_s[index + 1].value ) :
                    raise Exception( 'W_LegendreSeries.value = %s is >= next value = %s' % ( W_LegendreSeries_.value, self.W_LegendreSeries_s[index + 1].value ) )
            self.LegendreSeries_s[index] = W_LegendreSeries_

    def append( self, W_LegendreSeries_ ) :

        self[len( self )] = W_LegendreSeries_

    def copy( self, parent = None ) :

        n = V_W_XYs_LegendreSeries( self.axes, parent = parent )
        for i, W_LegendreSeries_ in enumerate( self ) : n[i] = W_LegendreSeries_
        return( n )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( physicalQuantityWithUncertainty.valueOrPQ( self.W_LegendreSeries_s[0].value, unitFrom = self.axes[0].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( physicalQuantityWithUncertainty.valueOrPQ( self.W_LegendreSeries_s[-1].value, unitFrom = self.axes[0].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomainUnit( self ) :

        return( self.axes[0].getUnit( ) )

    def maxLegendreOrder( self ) :

        lMax = 0
        for W_LegendreSeries_s in self : lMax = max( lMax, W_LegendreSeries_s.maxLegendreOrder( ) )
        return( lMax )

    def toPointwiseLinear( self, accuracy, cls = None ) :

        import V_W_XYs

        independentVY, dependentVY, qualifierVY = self.axes[0].interpolation.getInterpolationTokens( )
        if( independentVY not in [ axes.linearToken ] ) : raise Exception( 'vy independent interpolation = %s not supported' % independentVY )
        if( dependentVY not in [ axes.linearToken ] ) : raise Exception( 'vy dependent interpolation = %s not supported' % dependentVY )

        independentWY, dependentWY, qualifierWY = self.axes[1].interpolation.getInterpolationTokens( )  # WY interpolation check in W_XYs_LegendreSeries.toPointwiseLinear.

        axes_ = axes.axes( dimension = 4 )
        axes_[0] = axes.axis( self.axes[0].getLabel( ), 0, self.axes[0].getUnit( ), frame = self.axes[0].getFrame( ), \
            interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken, qualifierVY ) )
        axes_[1] = axes.axis( self.axes[1].getLabel( ), 1, self.axes[1].getUnit( ), frame = self.axes[1].getFrame( ), \
            interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken, qualifierWY ) )
        axes_[2] = axes.axis( "mu", 2, "", frame = self.axes[2].getFrame( ), \
            interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
        axes_[3] = axes.axis( "P", 3, self.axes[2].getUnit( ), frame = self.axes[2].getFrame( ) )

        if( cls is None ) : cls = V_W_XYs.V_W_XYs
        pwl = cls( axes_ )
        axesW_XYs = axes.referenceAxes( pwl, 3 )
        for w_xys in self :
            n = w_xys.toPointwiseLinear( accuracy, axes_ = axesW_XYs )
            n.value = w_xys.value
            pwl.append( n )
        return( pwl )

    def toXML( self, tag = 'xData', indent = '', incrementalIndent = '  ' ) :

        return( '\n'.join( self.toXMLList( indent = indent, incrementalIndent = incrementalIndent ) ) )

    def toXMLList( self, tag = 'xData', genre = 'V_W_XYs_LegendreSeries', indent = '', incrementalIndent = '  ' ) :

        indent2 = indent + incrementalIndent
        indent3 = indent2 + incrementalIndent
        xmlString = [ '%s<%s type="%s" length="%s">' % ( indent, tag, genre, len( self ) ) ]
        xmlString += self.axes.toXMLList( indent = indent2 )
        xmlString.append( '%s<data>' % indent2 )
        for W_LegendreSeries_ in self : xmlString += W_LegendreSeries_.toXMLList( tag = self.axes[0].getLabel( ), indent = indent3, includeType = False )
        xmlString[-1] += '</data></%s>' % tag
        return( xmlString )

if( __name__ == '__main__' ) :

    c_ls = XYs_LegendreSeries( '1/eV', [ 0.5, 0.3, -0.2, 0.1, 0.01 ] )
    print len( c_ls )
    print c_ls
    print c_ls.toXML( )

    pw = c_ls.toPointwiseLinear( 1e-3 )
    print pw.axes
    print pw
    print pw.toXML( )

    axes_ = axes.defaultAxes( dimension = 3, labelsUnits = { 0 : [ 'x', 'eV' ], 2 : [ 'z', '1/eV' ] } )
    w_xys_ls = W_XYs_LegendreSeries( axes_ )
    w_xys_ls[0] = XYs_LegendreSeries( '1/eV', [ 0.5,  0.3, -0.2,  0.1,  0.01 ], value = 1. )
    w_xys_ls[1] = XYs_LegendreSeries( '1/eV', [ 0.5, -0.3,  0.2, -0.1, -0.01 ], value = 2. )
    w_xys_ls[2] = XYs_LegendreSeries( '1/eV', [ 0.5,  0.3,  0.2,  0.1,  0.01 ], value = 3. )
    w_xys_ls[3] = XYs_LegendreSeries( '1/eV', [ 0.5, -0.3, -0.2, -0.1, -0.01 ], value = 5. )
    print w_xys_ls.toXML( )

    i = w_xys_ls.invert( )
    print i.toXML( )
