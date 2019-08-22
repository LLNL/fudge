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

import math
from fudge.core.utilities import brb
from xData.ancestry import ancestry
import axes, XYs, W_XYs
from fudge.core.math import fudgemath
from pqu import PQU

__metaclass__ = type

monikerXYs_LegendreSeries = 'XYs_LegendreSeries'
monikerW_XYs_LegendreSeries = 'W_XYs_LegendreSeries'
monikerV_W_XYs_LegendreSeries = 'V_W_XYs_LegendreSeries'

noExtrapolationToken = 'noExtrapolation'
flatExtrapolationToken = 'flatExtrapolation'

try :
    from numericalFunctions import Legendre as Legendre_C
    maxLegendreOrder = Legendre_C.maxMaxOrder( )
except :
    maxLegendreOrder = 64

def Legendre( n, x, checkXRange = True ) :
    """Returns the value of the Legendre function of order n at x.  For n <= 10, use analytical form,
    for n > 10 use the recursive relationship. (This would be way faster in C or using numpy version)"""

    if( n < 0 ) : raise ValueError( "\nError, n = %d < 0" % n )
    if( checkXRange and ( abs( x ) > 1 ) ) : raise ValueError( "Legendre: |x| > 1; x = %g" % x )
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

class XYs_LegendreSeries( ancestry ) :
    """
    This class stores and manipulates an angular pdf (i.e. pdf(mu) where mu is the cos of the angle) 
    represented as Legendre coefficients. The pdf and Legendre coefficients are related by

        pdf(mu) = Sum_over_l_of ( l + 0.5 ) * C_l * P_l(mu)

    where the sum is from l = 0 to lMax, lMax is the highest Legendre coefficients in the instance, C_l 
    is the Legendre coefficient for Legendre order l and P_l(mu) is the Legendre polynomial of order l.
    """

    moniker = monikerXYs_LegendreSeries

    def __init__( self, coefficients, index = None, value = None, parent = None ) :

        ancestry.__init__( self )
        self.index = index
        if( value is not None ) : value = fudgemath.toFloat( value )
        self.value = value
        self.coefficients = map( float, coefficients )

    def __len__( self ) :
        """Returns the number of Legendre coefficients in the instance (i.e., lMax + 1)."""

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

    def __add__( self, other ) :
        """Returns a XYs_LegendreSeries that is the sum of self and other. Other must be of type XYs_LegendreSeries."""

        if( not( isinstance( other, XYs_LegendreSeries ) ) ) : raise TypeError( 'other of type "%s"' % type( other ) )
        c_l1, c_l2 = self, other
        if( len( self ) < len( other ) ) : c_l1, c_l2 = other, self
        c_ls = c_l1.copy( )
        for l, c_l in enumerate( c_l2 ) : c_ls.coefficients[l] += c_l
        return( c_ls )

    def __sub__( self, other ) :
        """Returns a XYs_LegendreSeries that is the difference of self and other. Other must be of type XYs_LegendreSeries."""

        if( not( isinstance( other, XYs_LegendreSeries ) ) ) : raise TypeError( 'other of type "%s"' % type( other ) )
        c_l1, c_l2 = self.coefficients, other.coefficients
        if( len( self ) < len( other ) ) : c_l1, c_l2 = c_l2, c_l1
        c_ls = c_l1.copy( )
        for l, c_l in enumerate( c_l2 ) : c_ls.coefficients[l] += c_l
        return( c_ls )

    def __mul__( self, value ) :
        """Multiplies each coefficient of self by value. Value must be convertible to a float."""

        value_ = float( value )
        c_ls = self.copy( )
        for l, c_l in enumerate( self ) : c_ls[l] *= value
        return( c_ls )

    def __rmul__( self, value ) :
        "returns self.__mul__( value )."

        return( self.__mul__( value ) )

    def __str__( self ) :
        """Returns a string representation of the Legendre coefficients of self."""

        return( ' '.join( [ "%g" % c_l for c_l in self.coefficients ] ) )

    def copy( self, index = None, value = None, parent = None ) :
        """Creates a new XYs_LegendreSeries that is a copy of self except for its parent. The new 
        instance's index and value members are changes if index or value arguments are not None 
        respectively. The new instance's parent determined by the parent argument."""

        if( index is None ) : index = self.index
        if( value is None ) : value = self.value
        n = XYs_LegendreSeries( self.coefficients, index = index, value = value, parent = parent )
        return( n )

    def getCoefficientSafely( self, l ) :
        """
        Returns the (l+1)^th Legendre coefficient. Returns 0 if l is greater than lMax. This is like
        __getitem__ but allows for l to be greater than lMax.
        """

        if( l >= len( self ) ) : return( 0. )
        return( self.coefficients[l] )

    def getValue( self, mu ) :
        """Using the Legendre coefficients, this method calculates f(mu) and returns it."""

        P = 0.
        for l, c_l in enumerate( self.coefficients ) : P += ( l + 0.5 ) * c_l * Legendre( l, mu, checkXRange = False ) 
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

    def toPointwise_withLinearXYs( self, accuracy, biSectionMax = 16 ) :
        """
        The method constructs the pdf(mu) versus mu and returns it as a XYs instance. The accuracy of the 
        reconstruction (hence the number of points in the returned XYs) is determined by the accuracy argument.
        """

        if( accuracy < 1e-6 ) : accuracy = 1e-6
        if( accuracy > 0.1 ) : accuracy = 0.1

        try :
            from numericalFunctions import Legendre
            L = Legendre_C.Series( self.coefficients )
            P = L.toPointwiseLinear( accuracy, biSectionMax = biSectionMax, checkForRoots = True )
        except :
            P, n = [], 100
            for i in xrange( n ) :
                mu = -1. + ( 2. * i ) / n
                P.append( [ mu, self.getValue( mu ) ] )
            P.append( [ 1., self.getValue( 1. ) ] )
        axes_ = axes.axes( )
        axes_[0] = axes.axis( 'mu', 0, '', interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
        try :
            axes__ = self.findAttributeInAncestry( 'axes' )
            unit = axes__[-1].getUnit( )
        except :
            unit = ''
        axes_[1] = axes.axis( 'P(mu)', 1, unit )
        P = XYs.XYs( axes_, P, accuracy )
        return( P.thin( accuracy = accuracy ) )

    def toXML( self, tag = 'xData', indent = '', incrementalIndent = '  ' ) :
        """This method returns the XML string representation of self."""

        return( '\n'.join( self.toXMLList( tag = tag, indent = indent, incrementalIndent = incrementalIndent ) ) )

    def toXMLList( self, tag = 'xData', indent = '', incrementalIndent = '  ', includeType = True ) :
        """This method returns self as a list of string which converts to the XML string representation of self via
        the python statement '\n'.join( xmlString )."""

        indent2 = indent + incrementalIndent
        attrs = ''
        if( self.value is not None ) : attrs += ' value="%s"' % self.value
        if( self.index is not None ) : attrs += ' index="%s"' % self.index
        if( includeType ) : attrs += ' type="%s"' % self.moniker
        xmlString = [ '%s<%s%s length="%s">' % ( indent, tag, attrs, len( self ) ) ]
        xmlString += [ PQU.toShortestString(c_l) for c_l in self.coefficients ]
        xmlString[-1] += '</%s>' % tag
        return( [ ' '.join( xmlString ) ] )

class W_XYs_LegendreSeries( ancestry ) :

    moniker = monikerW_XYs_LegendreSeries
    xData = monikerW_XYs_LegendreSeries

    def __init__( self, axes_, index = None, value = None, parent = None, isPrimaryXData = True ) :

        ancestry.__init__( self )
        self.index = index
        if( value is not None ) : value = fudgemath.toFloat( value )
        self.value = value
        self.axes = axes_.copy( parent = self )
        self.LegendreSeries_s = []
        self.isPrimaryXData = isPrimaryXData

    def __len__( self ) :

        return( len( self.LegendreSeries_s ) )

    def __getitem__( self, index ) :

        return( self.LegendreSeries_s[index] )

    def __setitem__( self, index, LegendreSeries_ ) :
        self._setitem_internal( index, LegendreSeries_, copy=True )

    def _setitem_internal( self, index, LegendreSeries_, copy=True ) :

        if( not( isinstance( LegendreSeries_, XYs_LegendreSeries ) ) ) : 
            raise TypeError( 'right-hand-side must be instance of XYs_LegendreSeries; it is %s' % brb.getType( LegendreSeries_ ) )
        n = len( self )
        if( n < index ) : raise IndexError( 'index = %s while length is %s' % ( index, n ) )
        if( index < 0 ) : raise IndexError( 'index = %s' % index )
        if copy:
            LegendreSeries_ = LegendreSeries_.copy( index = index, value = LegendreSeries_.value, parent = self )
        else:
            LegendreSeries_.setAncestor( self )
            LegendreSeries_.index = index
            LegendreSeries_.isPrimaryXData = False
        if( n == 0 ) :
            self.LegendreSeries_s.append( LegendreSeries_ )
        elif( n == index ) :
            if( LegendreSeries_.value <= self.LegendreSeries_s[-1].value ) : 
                raise ValueError( 'LegendreSeries.value = %s is <= prior value = %s' % ( LegendreSeries_.value, self.LegendreSeries_s[-1].value ) )
            self.LegendreSeries_s.append( LegendreSeries_ )
        else :
            if( index > 0 ) :
                if( LegendreSeries_.value <= self.LegendreSeries_s[index - 1].value ) :
                    raise ValueError( 'LegendreSeries.value = %s is <= prior value = %s' % ( LegendreSeries_.value, self.LegendreSeries_s[index - 1].value ) )
            if( index < ( n - 1 ) ) :
                if( LegendreSeries_.value >= self.LegendreSeries_s[index + 1].value ) :
                    raise ValueError( 'LegendreSeries.value = %s is >= next value = %s' % ( LegendreSeries_.value, self.LegendreSeries_s[index + 1].value ) )
            self.LegendreSeries_s[index] = LegendreSeries_

    def append( self, LegendreSeries_, copy=True ) :

        self._setitem_internal( len(self), LegendreSeries_, copy )

    def copy( self, parent = None, index = None, value = None, moniker = None, isPrimaryXData = True, standAlone = False ) :
            # Name is not used here but needed for now for regions.

        if( index is None ) : index = self.index
        if( value is None ) : value = self.value
        axes_ = self.axes
        if( standAlone ) : axes_ = self.axes.copy( standAlone = standAlone )
        n = W_XYs_LegendreSeries( axes_, index = index, value = value, parent = parent, isPrimaryXData = isPrimaryXData )
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
            s2 = ( w - LegendreSeries1.value ) / ( LegendreSeries2.value - LegendreSeries1.value )
        elif( independent == axes.logToken ) :
            s2 = math.log( w / LegendreSeries1.value ) / math.log( LegendreSeries2.value / LegendreSeries1.value )
        else :
            raise Exception( 'Unsupported interpolation = %s' % independent  )
        if( dependent == axes.linearToken ) :
            coefficients = []
            s1 = 1.0 - s2
            for l1, Cl_1 in enumerate( coefficients1 ) :
                coefficients.append( s2 * coefficients2[l1] + s1 * Cl_1 )
        elif( dependent == axes.flatToken ) :
            coefficients = coefficients1
        else : raise Exception( 'Unsupported interpolation = %s' % dependent  )
        return( XYs_LegendreSeries( coefficients, value = w ) )

    def getCoefficientSafely( self, l, accuracy=1e-6 ) :
        """
        A version of XYs_LegendreSeries.getCoefficientSafely. 
        Returns the (l+1)^th Legendre coefficient as an XYs object of [ Ein, C_l] pairs. 
        C_l is set to 0 if l is greater than lMax. 
        The accuracy argument gets passed into the XYs object
        """
        if( accuracy < 1e-6 ) : accuracy = 1e-6
        if( accuracy > 0.1 ) : accuracy = 0.1
        axes_ = self.axes.copy( standAlone = True )
        axes_[1].label = "C_L=%i" % l
        return XYs.XYs( axes_, [ [ x.value, x.getCoefficientSafely( l ) ] for x in self ], accuracy )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( PQU.valueOrPQ( self.LegendreSeries_s[0].value, unitFrom = self.axes[0].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( PQU.valueOrPQ( self.LegendreSeries_s[-1].value, unitFrom = self.axes[0].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnit( self ) :

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

        pwl = self.toPointwise_withLinearXYs( accuracy )
        n = W_XYs_LegendreSeries( axes__ )
        for e_in in pwl : n.append( self.getLegendreSeriesAtW( e_in.value ) )
        return( n )

    def toPointwise_withLinearXYs( self, accuracy, axes_ = None ) :

        from fudge.core.math import miscellaneous

        def logFill( n, LS1, w_xys1, LS2, w_xys2, level = 0 ) :

            def getLSMid( ls1, ls2 ) :

                lsMid = XYs_LegendreSeries( ls1.coefficients )
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
            w_xysMid = lsMid.toPointwise_withLinearXYs( accuracy )
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
            axes_[1] = axes.axis( 'mu', 1, '', interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken, qualifier ) )
            axes_[2] = axes.axis( 'P', 2, axes__[1].getUnit( ) )
        n = W_XYs.W_XYs( axes_ )
        ls_prior, w_prior = None, None
        for w_LegendreSeries in self :
            w_xys = w_LegendreSeries.toPointwise_withLinearXYs( accuracy )
            w_xys.value = w_LegendreSeries.value
            if( w_prior is not None ) :
                if( dependent == axes.flatToken ) :
                    value = miscellaneous.shiftFloatDownABit( w_LegendreSeries.value, accuracy )
                    n.append( w_prior.copy( value = value ) )
                elif( dependent == axes.logToken ) :
                    logFill( n, ls_prior, w_prior, w_LegendreSeries, w_xys )
            n.append( w_xys )
            ls_prior, w_prior = w_LegendreSeries, w_xys
        return( n )

    def toXML( self, tag = 'xData', indent = '', incrementalIndent = '  ' ) :

        return( '\n'.join( self.toXMLList( tag = tag, indent = indent, incrementalIndent = incrementalIndent ) ) )

    def toXMLList( self, tag = 'xData', indent = '', incrementalIndent = '  ' ) :

        if( hasattr( self, 'tag' ) ) : tag = self.tag
        indent2 = indent + incrementalIndent
        extraStr = ''
        if( self.index is not None ) : extraStr += ' index="%s"' % self.index
        if( self.value is not None ) : extraStr += ' value="%s"' % self.value
        xDataString = ''
        if( self.isPrimaryXData ) : xDataString = ' xData="%s"' % self.xData
        extraXMLAttributeString = ''
        if( hasattr( self, 'extraXMLAttributeString' ) ) : extraXMLAttributeString = ' ' + self.extraXMLAttributeString( )
        xmlString = [ '%s<%s%s%s%s>' % ( indent, tag, extraStr, xDataString, extraXMLAttributeString ) ]
        xmlString += self.axes.toXMLList( indent = indent2 )
        for LegendreSeries_ in self : xmlString += LegendreSeries_.toXMLList( tag = self.axes[0].getLabel( ), indent = indent2, includeType = False )
        xmlString[-1] += '</%s>' % tag
        return( xmlString )

class V_W_XYs_LegendreSeries( ancestry ) :

    moniker = monikerV_W_XYs_LegendreSeries
    xData = monikerV_W_XYs_LegendreSeries

    def __init__( self, axes_, parent = None ) :

        ancestry.__init__( self )
        self.axes = axes_.copy( parent = self )
        self.W_LegendreSeries_s = []

    def __len__( self ) :

        return( len( self.W_LegendreSeries_s ) )

    def __getitem__( self, index ) :

        return( self.W_LegendreSeries_s[index] )

    def __setitem__( self, index, W_LegendreSeries_ ) :

        self._setitem_internal( index, W_LegendreSeries_, copy=True )

    def _setitem_internal( self, index, W_LegendreSeries_, copy=True ) :

        if( not( isinstance( W_LegendreSeries_, W_XYs_LegendreSeries ) ) ) : 
            raise Exception( 'right-hand-side must be instance of W_XYs_LegendreSeries; it is %s' % brb.getType( W_LegendreSeries_ ) )
        n = len( self )
        if( n < index ) : raise IndexError( 'index = %s while length is %s' % ( index, n ) )
        if( index < 0 ) : raise IndexError( 'index = %s' % index )
        if copy:
            W_LegendreSeries_ = W_LegendreSeries_.copy( index = index, value = W_LegendreSeries_.value, parent = self, isPrimaryXData = False )
        else:
            W_LegendreSeries_.setAncestor( self )
            W_LegendreSeries_.index = index
            W_LegendreSeries_.isPrimaryXData = False
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

    def append( self, W_LegendreSeries_, copy=True ) :

        self._setitem_internal( len(self), W_LegendreSeries_, copy=copy )

    def copy( self, parent = None ) :

        n = V_W_XYs_LegendreSeries( self.axes, parent = parent )
        for i, W_LegendreSeries_ in enumerate( self ) : n[i] = W_LegendreSeries_
        return( n )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( PQU.valueOrPQ( self.W_LegendreSeries_s[0].value, unitFrom = self.axes[0].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( PQU.valueOrPQ( self.W_LegendreSeries_s[-1].value, unitFrom = self.axes[0].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnit( self ) :

        return( self.axes[0].getUnit( ) )

    def getW_XYs_LegendreSeriesAtV( self, v, extrapolation = noExtrapolationToken ) :
        """Returns the W_XYs_LegendreSeries for self at v. Currently, extrapolation is ignored."""

        eps = 1e-12
        smallerEps = 1e-2 * eps         # Must be greater than machine precision.

        if( self[0].value >= v ) : return( self[0] )
        if( self[-1].value <= v ) : return( self[-1] )
        for W_XYs_LegendreSeries2 in self :
            if( W_XYs_LegendreSeries2.value >= v ) : break
            W_XYs_LegendreSeries1 = W_XYs_LegendreSeries2
        if( W_XYs_LegendreSeries2.value == v ) : return( W_XYs_LegendreSeries2 )
        independent, dependent, qualifier = self.axes[0].interpolation.getInterpolationTokens( )
        if( independent != axes.linearToken ) : raise Exception( 'Interpolation %s not support for independent axis.' % independent )
        if(   dependent != axes.linearToken ) : raise Exception( 'Interpolation %s not support for dependent axis.' % dependent )

        def getEps( Eps ) :

            dE = Eps[-1] - Eps[0]
            Eps_unitBase = [ ( Ep - Eps[0] ) / dE for Ep in Eps[:-1] ]
            Eps_unitBase.append( 1. )           # Make sure last point is 1. First point should not be an issue.
            return( Eps, Eps_unitBase )

        def getLs( Ep, Eps, W_XYs_LegendreSeries_ ) :

            Ep_ = Ep * Eps[-1] + ( 1 - Ep ) * Eps[0]
            for Ep in Eps :
                if( abs( Ep - Ep_ ) < smallerEps * Ep ) :
                    Ep_ = Ep
                    break
            Ep_ = min( Eps[-1], max( Eps[0], Ep_ ) )
            return( W_XYs_LegendreSeries_.getLegendreSeriesAtW( Ep_ ) )

        Eps1, Eps1_unitBase = getEps( [ Ls.value for Ls in W_XYs_LegendreSeries1 ] )
        Eps2, Eps2_unitBase = getEps( [ Ls.value for Ls in W_XYs_LegendreSeries2 ] )
        for Ep1 in Eps1_unitBase :
            if( Ep1 not in Eps2_unitBase ) : Eps2_unitBase.append( Ep1 )
        Eps2_unitBase.sort( )
        Ep_p, iEsToDelete = -1, []
        for iE, Ep in enumerate( Eps2_unitBase ) :      # Remove any close points.
            if( Ep - Ep_p < 1e-12 * Ep ) : iEsToDelete.insert( 0, iE )
            Ep_p = Ep
        for iE in iEsToDelete : del Eps2_unitBase[iE]

        n1 = W_XYs_LegendreSeries( W_XYs_LegendreSeries1.axes, value = v )
        f1 = ( W_XYs_LegendreSeries2.value - v ) / ( W_XYs_LegendreSeries2.value - W_XYs_LegendreSeries1.value )
        EpMin, EpMax = f1 * Eps1[0] + ( 1 - f1 ) * Eps2[0], f1 * Eps1[-1] + ( 1 - f1 ) * Eps2[-1]
        f2 = ( 1 - f1 ) * ( Eps2[-1] - Eps2[0] ) / ( EpMax - EpMin )
        f1 *= ( Eps1[-1] - Eps1[0] ) / ( EpMax - EpMin )
        for Ep in Eps2_unitBase :
            XYs_LegendreSeries1 = getLs( Ep, Eps1, W_XYs_LegendreSeries1 )
            XYs_LegendreSeries2 = getLs( Ep, Eps2, W_XYs_LegendreSeries2 )
            n1.append( XYs_LegendreSeries( f1 * XYs_LegendreSeries1 + f2 * XYs_LegendreSeries2, value = ( 1 - Ep ) * EpMin + Ep * EpMax ) )
        return( n1 )

    def maxLegendreOrder( self ) :

        lMax = 0
        for W_LegendreSeries_s in self : lMax = max( lMax, W_LegendreSeries_s.maxLegendreOrder( ) )
        return( lMax )

    def toPointwise_withLinearXYs( self, accuracy, cls = None ) :

        import V_W_XYs

        independentVY, dependentVY, qualifierVY = self.axes[0].interpolation.getInterpolationTokens( )
        if( independentVY not in [ axes.linearToken ] ) : raise Exception( 'vy independent interpolation = %s not supported' % independentVY )
        if( dependentVY not in [ axes.linearToken ] ) : raise Exception( 'vy dependent interpolation = %s not supported' % dependentVY )

                                                        # WY interpolation checked in W_XYs_LegendreSeries.toPointwise_withLinearXYs.
        independentWY, dependentWY, qualifierWY = self.axes[1].interpolation.getInterpolationTokens( )

        axes_ = axes.axes( dimension = 4 )
        axes_[0] = axes.axis( self.axes[0].getLabel( ), 0, self.axes[0].getUnit( ), \
            interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken, qualifierVY ) )
        axes_[1] = axes.axis( self.axes[1].getLabel( ), 1, self.axes[1].getUnit( ), \
            interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken, qualifierWY ) )
        axes_[2] = axes.axis( "mu", 2, "", interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
        axes_[3] = axes.axis( "P", 3, self.axes[2].getUnit( ) )

        if( cls is None ) : cls = V_W_XYs.V_W_XYs
        pwl = cls( axes_, self.getProductFrame() )
        axesW_XYs = axes.referenceAxes( pwl, 3 )
        for w_xys in self :
            n = w_xys.toPointwise_withLinearXYs( accuracy, axes_ = axesW_XYs )
            n.value = w_xys.value
            pwl.append( n )
        return( pwl )

    def toXML( self, tag = 'xData', indent = '', incrementalIndent = '  ' ) :

        return( '\n'.join( self.toXMLList( tag = tag, indent = indent, incrementalIndent = incrementalIndent ) ) )

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

    c_ls = XYs_LegendreSeries( [ 0.5, 0.3, -0.2, 0.1, 0.01 ] )
    print len( c_ls )
    print c_ls
    print c_ls.toXML( )

    pw = c_ls.toPointwise_withLinearXYs( 1e-3 )
    print pw.axes
    print pw
    print pw.toXML( )
    pw.plot( )

    axes_ = axes.defaultAxes( dimension = 3, labelsUnits = { 0 : [ 'x', 'eV' ], 2 : [ 'z', '1/eV' ] } )
    w_xys_ls = W_XYs_LegendreSeries( axes_ )
    w_xys_ls[0] = XYs_LegendreSeries( [ 0.5,  0.3, -0.2,  0.1,  0.01 ], value = 1. )
    w_xys_ls[1] = XYs_LegendreSeries( [ 0.5, -0.3,  0.2, -0.1, -0.01 ], value = 2. )
    w_xys_ls[2] = XYs_LegendreSeries( [ 0.5,  0.3,  0.2,  0.1,  0.01 ], value = 3. )
    w_xys_ls[3] = XYs_LegendreSeries( [ 0.5, -0.3, -0.2, -0.1, -0.01 ], value = 5. )
    print w_xys_ls.toXML( )

    i = w_xys_ls.invert( )
    print i.toXML( )
