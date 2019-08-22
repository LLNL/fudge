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

__metaclass__ = type

import abc
import math

from pqu import PQU as PQUModule

import base as baseModule
import standards as standardsModule
import values as valuesModule
import axes as axesModule
import XYs as XYsModule
import uncertainties as uncertaintiesModule

from numericalFunctions import pointwiseXY_C as pointwiseXY_CModule
floatToShortestString = pointwiseXY_CModule.floatToShortestString

from numericalFunctions import Legendre as Legendre_C
maxLegendreOrder = Legendre_C.maxMaxOrder( )


def Legendre( n, mu, checkXRange = True ) :
    """
    Returns the value of the Legendre function of order n at mu using the recursive relationship.
    """

    if( n < 0 ) : raise ValueError( "\nError, n = %d < 0" % n )
    if( checkXRange and ( abs( mu ) > 1 ) ) : raise ValueError( "Legendre: |mu| > 1; mu = %g" % mu )
    Pn = 0.
    Pnp1 = 1.
    n_ = 0
    twoNp1 = 1
    while( n_ < n ) :
        Pnm1 = Pn
        Pn = Pnp1
        n_p1 = n_ + 1
        Pnp1 = ( twoNp1 * mu * Pn - n_ * Pnm1 ) / n_p1
        twoNp1 += 2
        n_ = n_p1
    return( Pnp1 )


class series( baseModule.xDataFunctional ) :
    """
    This class is the base class for storing a 1d function as a polynomial series. The function store the 
    polynomial coefficients and has methods for manipulations of the coefficients that are generic to
    all polynomial series (e.g., simple polynomial, Legendre).
    """

    dimension = 1
    __metaclass__ = abc.ABCMeta
    ancestryMembers = baseModule.xDataFunctional.ancestryMembers # + ( 'coefficients', )

    def __init__( self, coefficients, domainMin, domainMax, lowerIndex = 0, axes = None,
            index = None, valueType = standardsModule.types.float64Token, value = None, label = None, sep = ' ' ) :

        baseModule.xDataFunctional.__init__( self, self.moniker, axes, index = index, valueType = valueType,
                value = value, label = label )

        if( not( isinstance( sep, str ) ) ) : raise TypeError( 'sep must be a str instance' )
        if( len( sep ) != 1 ) : raise TypeError( 'sep length must be 1 not %d' % len( sep ) )
        self.__sep = sep

        self.__domainMin = float( domainMin )
        self.__domainMax = float( domainMax )
        self.__lowerIndex = int( lowerIndex )

        self.coefficients = map( float, coefficients )

    def __len__( self ) :
        """Returns the number of Legendre coefficients in the instance (e.g., for Legendre series it is lMax + 1)."""

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
        """Returns a series that is the sum of self and other. Other must be of type series."""

        try :
            other = float( other )
            c_ls = self.copy( )
            for l, c_l in enumerate( c_ls ) : c_ls.coefficients[l] += other
        except :
            self.checkSameSeriesType( other, 'add' )
            c_l1, c_l2 = self, other
            if( len( self ) < len( other ) ) : c_l1, c_l2 = other, self
            c_ls = c_l1.copy( )
            for l, c_l in enumerate( c_l2 ) : c_ls.coefficients[l] += c_l
        return( c_ls )

    def __sub__( self, other ) :
        """Returns a series that is the difference of self and other. Other must be of type series."""

        try :
            other = float( other )
            c_ls = self.copy( )
            for l, c_l in enumerate( c_ls ) : c_ls.coefficients[l] -= other
        except :
            self.checkSameSeriesType( other, 'subtract' )
            c_l1, c_l2 = self.coefficients, other.coefficients
            if( len( self ) < len( other ) ) : c_l1, c_l2 = c_l2, c_l1
            c_ls = c_l1.copy( )
            for l, c_l in enumerate( c_l2 ) : c_ls.coefficients[l] += c_l
        return( c_ls )

    def __mul__( self, value ) :
        """Returns a new series that is each coefficient of self multiplied by value. Value must be convertible to a float."""

        value_ = float( value )
        c_ls = self.copy( )
        for l, c_l in enumerate( self ) : c_ls[l] *= value
        return( c_ls )

    __rmul__ = __mul__

    def __str__( self ) :
        """Returns a string representation of the Legendre coefficients of self."""

        return( ' '.join( [ "%g" % c_l for c_l in self.coefficients ] ) )

    def evaluate(self, x):
        """Override this in derived classes if you have a faster approach"""
        total = 0.0
        for i,c in enumerate(self.coefficients):
            total += c * self.evaluateBasisFunction(x,i)
        return total

    def evaluateBasisFunction(self, x, i):
        raise NotImplementedError("Implement in derived classes")

    def checkSameSeriesType( self, other, operator ) :

        if( not( isinstance( other, series ) ) ) : raise TypeError( 'other of type "%s"' % type( other ) )
        if( self.moniker != other.moniker ) : raise TypeError( 'Cannot %s series %s to series %s' % ( operator, self.moniker, other.moniker ) )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        if( self.axes is None ) : print self.toXLink( )
        factors = self.axes.convertUnits( unitMap )
        if( factors[:2] !=[ 1., 1. ] ) :
            self.__domainMin *= factors[1]
            self.__domainMax *= factors[1]
            for l, c_l in enumerate( self ) : self[l] *= factors[0]
            factor = math.pow( factors[1], -self.__lowerIndex )
            for l, c_l in enumerate( self ) :
                if( ( l + self.__lowerIndex ) == 0 ) : factor = 1
                self[l] *= factor
                factor *= 1 / factors[1]
        self.fixValuePerUnitChange( factors )

        if( self.uncertainties ) : self.uncertainties.convertUnits( unitMap )

    def copy( self ) :
        """
        Creates a new series that is a copy of self. The new 
        instance's index and value members are changes if index or value arguments are not None 
        respectively.
        """

        return( self.returnAsClass( self, self.coefficients, index = self.index, value = self.value ) )

    __copy__ = copy
    __deepcopy__ = __copy__

    @property
    def domainMin( self ) :

        return( self.__domainMin )

    @property
    def domainMax( self ) :

        return( self.__domainMax )

    @property
    def lowerIndex( self ) :

        return( self.__lowerIndex )

    @property
    def upperIndex( self ) :

        return( self.__lowerIndex + len( self ) )

    @property
    def sep( self ) :

        return( self.__sep )

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    def getCoefficientSafely( self, l ) :
        """
        Returns the (l+1)^th Legendre coefficient. Returns 0 if l is greater than lMax. This is like
        __getitem__ but allows for l to be greater than lMax.
        """

        if( l >= len( self ) ) : return( 0. )
        return( self.coefficients[l] )

    def invert( self ) :
        """
        This method returns a series instance that is the mirror of self about mu = 0.
        That is, returns a Legendre series which represent self's pdf(-mu).
        """

        series_ = self.copy( )
        for l in xrange( 1, len( series_ ), 2 ) : series_.coefficients[l] *= -1
        return( series_ )

    def setData( self, data ) :

        self.coefficients = data

    def rangeMin( self, unitTo = None ) :

        raise NotImplementedError( )

    def rangeMax( self, unitTo = None ) :

        raise NotImplementedError( )

    def toString( self ) :

        str = [ "%d %16.8g" % ( l, coefficient ) for l, coefficient in enumerate( self ) ]
        return( '\n'.join( str ) )

    def toXML( self, indent = '', **kwargs ) :
        """This method returns the XML string representation of self."""

        return( '\n'.join( self.toXMLList( **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :
        """
        This method returns self as a list of strings which converts to the XML string representation of self via
        the python statement '\n'.join( XMLList ).
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        valueFormatter = kwargs.get( 'valueFormatter', floatToShortestString )
        significantDigits = kwargs.get( 'significantDigits', 15 )

        attributesStr = baseModule.xDataFunctional.attributesToXMLAttributeStr( self )
        if( None in self.fixedDomain( ) ) : attributesStr += ' domainMin="%s" domainMax="%s"' % (
                valueFormatter( self.domainMin, significantDigits = significantDigits ),
                valueFormatter( self.domainMax, significantDigits = significantDigits ) )
        if( self.lowerIndex != 0 ) : attributesStr += ' lowerIndex="%s"' % self.lowerIndex

        # FIXME: converting self.coefficients to values for printing. Should it be stored as values in the first place?
        coefs = valuesModule.values( self.coefficients, valueType = self.valueType, sep = self.__sep )

        XMLList = [ '%s<%s%s>' % ( indent, self.moniker, attributesStr) ]
        if( self.isPrimaryXData( ) and ( self.axes is not None ) ) :
            XMLList += self.axes.toXMLList( indent2, **kwargs )
            XMLList += coefs.toXMLList( indent2, **kwargs )
            if self.uncertainties:
                XMLList.append( "%s<uncertainties>" % indent2 )
                indent3 = indent2 + kwargs.get( 'incrementalIndent', '  ' )
                for uncertainty in self.uncertainties:
                    XMLList += uncertainty.toXMLList( indent3, **kwargs )
                XMLList[-1] += "</uncertainties>"
            XMLList[-1] += '</%s>' % self.moniker
            return XMLList
        else:
            XMLList += coefs.toXMLList( '', **kwargs )
            XMLList[-1] += '</%s>' % self.moniker
            return( [ ''.join( XMLList ) ] )

    @classmethod
    def returnAsClass( cls, self, coefficients, index = None, value = None ) :

        return( cls( coefficients, self.domainMin, self.domainMax, lowerIndex = self.lowerIndex, axes = self.axes,
            index = index, valueType = self.valueType, value = value, label = self.label, sep = self.__sep ) )

    @classmethod
    def parseXMLNode( cls, xDataElement, xPath, linkData, axes = None, **kwargs ) :
        """
        Translate a series XML element into its python xData class.
        """

        xPath.append( xDataElement.tag )

        domainMin, domainMax = cls.fixedDomain( )
        attrs = { 'label' : None, 'index' : None, 'value' : None, 'domainMin' : domainMin, 'domainMax' : domainMax, 'lowerIndex' : 0 }
        attributes = { 'label' : str, 'index' : int, 'value' : float, 'domainMin' : float, 'domainMax' : float, 'lowerIndex' : int }
        if xDataElement.find('axes') is not None:
            axes = axesModule.axes.parseXMLNode( xDataElement.find('axes'), xPath, linkData )
        for key, item in xDataElement.items( ) :
            if( key == 'axes' ) : continue
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key]( item )
        if( attrs['domainMin'] == None ) : raise ValueError( 'missing attribute "domainMin"' )
        if( attrs['domainMax'] == None ) : raise ValueError( 'missing attribute "domainMax"' )
        coefficients = map( float, xDataElement.find('values').text.split() )
        series = cls( coefficients = coefficients, axes = axes, **attrs )
        uncertElement = xDataElement.find(uncertaintiesModule.uncertainties.moniker)
        if uncertElement is not None:
            series.uncertainties = uncertaintiesModule.uncertainties.parseXMLNode(uncertElement, xPath, linkData)
        xPath.pop( )
        return( series )

    @classmethod
    def parseXMLString( cls, XMLString ) :

        from xml.etree import cElementTree

        return( cls.parseXMLNode( cElementTree.fromstring( XMLString ), [], [] ) )

    @staticmethod
    def fixedDomain( ) :

        return( None, None )


class LegendreSeries( series ) :
    """
    This class represent a Legendre series for a function f(mu) as:

    ..math::

        f(\mu) = \sum_L ( l + 0.5 ) * C_l * P_l(\mu)

    so
    ..math::
        C_l=\int_{-1}^1 d\mu P_l(\mu) f(\mu)

    where the sum is from l = 0 to lMax, lMax is the highest Legendre coefficients in the instance, C_l
    is the Legendre coefficient for Legendre order l and P_l(mu) is the Legendre polynomial of order l.
    This class stores the Legendre coefficients C_l.
    """

    moniker = 'Legendre'
    dimension = 1

    def __init__( self, coefficients, domainMin = -1, domainMax = 1, lowerIndex = 0, axes = None,
            index = None, valueType = standardsModule.types.float64Token, value = None, label = None, sep = ' ' ) :

        if( lowerIndex != 0 ) : raise ValueError( 'lowerIndex = %s must be 0' )
        if( ( domainMin != -1 ) or ( domainMax != 1 ) ) : raise ValueError( 'domain must be [-1, 1], not [%s, %s]' % ( domainMin, domainMax ) )
        series.__init__( self, coefficients, -1, 1, axes = axes, index = index, valueType = valueType, value = value, 
                label = label, sep = sep )

    def evaluate( self, mu ) :
        """Using the Legendre coefficients, this method calculates f(mu) and returns it."""

        P = 0.
        for l, c_l in enumerate( self.coefficients ) : P += ( l + 0.5 ) * c_l * Legendre( l, mu, checkXRange = False ) 
        return( P )

    def evaluateBasisFunction(self, mu, l):
        return ( l + 0.5 ) * Legendre( l, mu, checkXRange = False )

    def isIsotropic( self ) :
        """Returns True if self is isotropic."""

        for coefficient in self[1:] :
            if( coefficient != 0. ) : return( False )
        return( True )

    def normalize( self, insitu = False, dimension = 1 ) :
        """
        The dimension argument is ignored. Only here to be compatable with calling from XYsnd.normalize.
        """

        norm = self.coefficients[0]
        coefficients = [ coefficient / norm for coefficient in self ]
        if( insitu ) :
            self.setData( coefficients )
            return( self )
        return( self.returnAsClass( self, coefficients ) )

    def toLinearXYsClass( self ) :

        return( XYsModule.XYs1d )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        This method constructs the pdf(mu) versus mu and returns it as a XYs1d instance. The accuracy of the 
        reconstruction (hence the number of points in the returned XYs1d) is determined by the accuracy argument.
        """

        arguments = self.getArguments( kwargs, { 'accuracy' : 1e-3, 'biSectionMax' : 16 } )
        accuracy = arguments['accuracy']
        biSectionMax = arguments['biSectionMax']

        if( accuracy < 1e-6 ) : accuracy = 1e-6
        if( accuracy > 0.1 ) : accuracy = 0.1

        try :
            L = Legendre_C.Series( self.coefficients )
            P = L.toPointwiseLinear( accuracy, biSectionMax = biSectionMax, checkForRoots = True )
        except :
            P, n = [], 400
            for i in xrange( n ) :
                mu = -1. + ( 2. * i ) / n
                P.append( [ mu, self.evaluate( mu ) ] )
            P.append( [ 1., self.evaluate( 1. ) ] )
        axes = axesModule.axes( )
        unit = self.getAxisUnitSafely( 0 )
        axes[0] = axesModule.axis( 'P(mu)', 0, unit )
        axes[1] = axesModule.axis( 'mu', 1, '' )
        Pclass = self.toLinearXYsClass()
        P = Pclass( P, axes = axes )
        return( P.thin( accuracy = accuracy ) )

    @staticmethod
    def fixedDomain( ) :

        return( -1, 1 )


class polynomial1d( series ) :

    moniker = 'polynomial1d'
    dimension = 1

    def __init__( self, coefficients, domainMin, domainMax, lowerIndex = 0, axes = None,
            index = None, valueType = standardsModule.types.float64Token, value = None, label = None, sep = ' ' ) :

        series.__init__( self, coefficients, domainMin, domainMax, lowerIndex = lowerIndex, axes = axes, index = index, 
                valueType = valueType, value = value, label = label, sep = sep )

    def evaluate( self, x ) :
        """Using the polynomial coefficients, this method calculates p(x) and returns it."""

        P = 0.
        for c_l in reversed( self.coefficients ) : P = c_l  + x * P
        return( P )

    def evaluateBasisFunction(self, x, i):
        return pow(x,i)

    def toLinearXYsClass( self ) :

        return( XYsModule.XYs1d )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        This method constructs the y(x) versus x and returns it as a XYs1d instance. The accuracy of the 
        reconstruction (hence the number of points in the returned XYs1d) is determined by the accuracy argument.
        Currently, accuracy is not implemented.
        """

        arguments = self.getArguments( kwargs, { 'accuracy' : 1e-3, 'biSectionMax' : 16 } )
        accuracy = arguments['accuracy']
        biSectionMax = arguments['biSectionMax']

        if( accuracy < 1e-6 ) : accuracy = 1e-6
        if( accuracy > 0.1 ) : accuracy = 0.1

        P, n = [], 1000
        for i in xrange( n + 1 ) :
            x = ( ( n - i ) * self.domainMin + self.domainMax * i ) / n
            P.append( [ x, self.evaluate( x ) ] )
        axes = axesModule.axes( )
        yUnit = self.getAxisUnitSafely( 0 )
        xUnit = self.getAxisUnitSafely( 1 )
        axes[0] = axesModule.axis( 'y(x)', 0, yUnit )       # FIXME
        axes[1] = axesModule.axis( 'x', 1, xUnit )          # FIXME
        Pclass = self.toLinearXYsClass()
        P = Pclass( P, axes = axes )
        return( P.thin( accuracy = accuracy ) )

class linearSpline1d( series ) :
    """
    This class is a simple linear spline.  It basically wraps the XYs1d class.
    Linear interpolation uses the linear spline or the "hat" basis.  The first basis function looks like this::
        |\
        | \
        |  \
        ----
    The second basis function looks like this::
          /\
         /  \
        /    \
        ------
    and so on.  Together y(x) = \sum_i B_i(x) y_i such that y_i = y(x_i)
    """

    moniker = 'linearSpline1d'
    dimension = 1

    def __init__( self, xdata, ydata, axes,
            index = None, valueType = standardsModule.types.float64Token, value = None, label = None, sep = ' ' ) :
        if len(xdata) != len(ydata): raise ValueError("Number of x and y values not equal")

        series.__init__( self, ydata, xdata[0], xdata[-1], lowerIndex = 0, axes = axes, index = index,
                valueType = valueType, value = value, label = label, sep = sep )

        self.axes=axes
        self.basis = XYsModule.XYs1d(axes=self.axes, data=zip(xdata,ydata), interpolation=standardsModule.interpolation.linlinToken)

    def evaluateBasisFunction(self, x, i):

        self.basis[i] = (self.basis[i][0], 1.0)
        result=self.basis.evaluate(x)
        self.basis[i] = (self.basis[i][0], 0.0)
        return result

    def toLinearXYsClass( self ) :

        return( XYsModule.XYs1d )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        This method constructs the y(x) versus x and returns it as a XYs1d instance. Basically we just copy the basis
        function widget and put the y values back.
        """
        result = self.basis.copy()
        for i,c in enumerate(self.coefficients):
            result[i]=(result[i][0],c)
        return result
