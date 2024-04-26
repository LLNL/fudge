# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import abc
import math

from pqu import PQU as PQUModule

from numericalFunctions import pointwiseXY_C as pointwiseXY_CModule
floatToShortestString = pointwiseXY_CModule.floatToShortestString

from numericalFunctions import Legendre as Legendre_C
maxLegendreOrder = Legendre_C.maxMaxOrder( )

from . import enums as enumsModule
from . import base as baseModule
from . import values as valuesModule
from . import axes as axesModule
from . import XYs1d as XYs1dModule

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


class Series1d( baseModule.XDataFunctional ) :
    """
    This class is the base class for storing a 1d function as a polynomial series. The function store the 
    polynomial coefficients and has methods for manipulations of the coefficients that are generic to
    all polynomial series (e.g., simple polynomial, Legendre).
    """

    dimension = 1
    ancestryMembers = baseModule.XDataFunctional.ancestryMembers # + ( 'coefficients', )

    def __init__(self, coefficients, domainMin, domainMax, lowerIndex = 0, axes = None,
            index = None, valueType = enumsModule.ValueType.float64, outerDomainValue = None, label = None):

        baseModule.XDataFunctional.__init__(self, axes, index=index, valueType=valueType, outerDomainValue=outerDomainValue, label=label)

        self.__domainMin = float( domainMin )
        self.__domainMax = float( domainMax )
        self.__lowerIndex = int( lowerIndex )

        self.coefficients = coefficients

    def __len__( self ) :
        """Returns the number of coefficients in the instance (e.g., for Legendre series it is lMax + 1)."""

        return( len( self.__coefficients ) )

    def __getitem__( self, l ) :
        """Returns the (l+1)^th coefficient."""

        return( self.__coefficients[l] )

    def __setitem__( self, l, c_l ) :
        """Sets the (l+1)^th coefficient to c_l. l must be less than or equal to lMax."""

        if( l == len( self ) ) :
            self.__coefficients.append( float( c_l ) )
        else :
            self.__coefficients[l] = float( c_l ) 

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
            c_l1, c_l2 = self.__coefficients, other.coefficients
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

    def __div__(self, value):
        '''Returns a new series that is each coefficient of self divided by value. Value must be convertible to a float.'''

        value_ = float(value)
        c_ls = self.copy()
        for l, c_l in enumerate(self):
            c_ls[l] /= value
        return c_ls 

    def __eq__(self, other):
        """Returns true if self and other have the same coefficients, domain and lowerIndex."""

        return self.areDomainsMutual(other) and self.lowerIndex == other.lowerIndex and self.coefficients == other.coefficients

    def __str__(self):
        """Returns a string representation of the coefficients of self."""

        return ' '.join(["%g" % c_l for c_l in self.__coefficients])

    @property
    def coefficients(self):

        return self.__coefficients

    @coefficients.setter
    def coefficients(self, coefficients):

        self.__coefficients = list(map(float, coefficients))

    def areDomainsMutual(self, other):
        '''
        Returns True if domain mins are the same for *self* and *other*, and domain maxs are the same for *self* and *other*;
        otherwise, returns False.

        :param other:       Another xData 1d instance.
        '''

        return self.domainMin == other.domainMin and self.domainMax == other.domainMax

    def evaluate(self, mu):
        """Evaluates the Legendre series at the mu."""

        total = 0.0
        for i,c in enumerate(self.__coefficients):
            total += c * self.evaluateBasisFunction(mu, i)
        return total

    def evaluateBasisFunction(self, x, i):
        raise NotImplementedError("Implement in derived classes")

    def checkSameSeriesType( self, other, operator ) :

        if( not( isinstance( other, Series1d ) ) ) : raise TypeError( 'other of type "%s"' % type( other ) )
        if( self.moniker != other.moniker ) : raise TypeError( 'Cannot %s series %s to series %s' % ( operator, self.moniker, other.moniker ) )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        if len(self.axes) == 0: return
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

        if( self.uncertainty ) : self.uncertainty.convertUnits( unitMap )

    def copy( self ) :
        """
        Creates a new series that is a copy of self. The new 
        instance's index and value members are changes if index or value arguments are not None 
        respectively.
        """

        return( self.returnAsClass( self, self.__coefficients, index = self.index, outerDomainValue = self.outerDomainValue ) )

    __copy__ = copy

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

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    def getCoefficientSafely( self, l ) :
        """
        Returns the (l+1)^th coefficient. Returns 0 if l is greater than lMax. This is like
        __getitem__ but allows for l to be greater than lMax.
        """

        if( l >= len( self ) ) : return( 0. )
        return( self.__coefficients[l] )

    def hasData(self):
        '''
        Returns True if the length of self's coefficients is greater than 0.
        '''

        return len(self.__coefficients) > 0

    def invert( self ) :
        """
        This method returns a series instance that is the mirror of self about x1 = 0 (i.e., negates all odd coefficients).
        For example, for a Legendre series return equilavent of pdf(-mu).
        """

        series_ = self.copy( )
        for l in range( 1, len( series_ ), 2 ) : series_.coefficients[l] *= -1
        return( series_ )

    def setData( self, data ) :

        self.__coefficients = data

    @property
    def rangeMin( self ) :

        raise NotImplementedError( )

    @property
    def rangeMax( self ) :

        raise NotImplementedError( )

    def toString( self ) :

        str = [ "%d %16.8g" % ( l, coefficient ) for l, coefficient in enumerate( self.__coefficients ) ]
        return( '\n'.join( str ) )

    def toXML_strList(self, indent = '', **kwargs):
        """
        This method returns self as a list of strings.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        valueFormatter = kwargs.get('valueFormatter', floatToShortestString)
        significantDigits = kwargs.get('significantDigits', 15)

        attributesStr = baseModule.XDataFunctional.attributesToXMLAttributeStr(self)
        if None in self.fixedDomain(): attributesStr += ' domainMin="%s" domainMax="%s"' % (
                valueFormatter(self.domainMin, significantDigits = significantDigits),
                valueFormatter(self.domainMax, significantDigits = significantDigits) )
        if self.lowerIndex != 0: attributesStr += ' lowerIndex="%s"' % self.lowerIndex

        # FIXME2: converting self.__coefficients to values for printing. Should it be stored as values in the first place?
        coefs = valuesModule.Values(self.__coefficients, valueType = self.valueType)

        XML_strList = [ '%s<%s%s>' % ( indent, self.moniker, attributesStr) ]
        if self.isPrimaryXData() and len(self.axes) > 0:
            XML_strList += self.axes.toXML_strList(indent2, **kwargs)
            XML_strList += coefs.toXML_strList(indent2, **kwargs)
            if self.uncertainty: XML_strList += self.uncertainty.toXML_strList(indent2, **kwargs)
            XML_strList[-1] += '</%s>' % self.moniker
        else:
            XML_strList += coefs.toXML_strList('', **kwargs)
            XML_strList[-1] += '</%s>' % self.moniker
            XML_strList = [ ''.join(XML_strList) ]

        return XML_strList

    @classmethod
    def returnAsClass( cls, self, coefficients, index = None, outerDomainValue = None ) :

        return cls(coefficients, self.domainMin, self.domainMax, lowerIndex = self.lowerIndex, axes = self.axes, index = index, 
                valueType = self.valueType, outerDomainValue = outerDomainValue, label = self.label)

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Translate a series XML element into its python xData class.
        """

        attributes, extraAttributes = baseModule.XDataFunctional.parseBareNodeCommonAttributes(node, xPath)     # parseBareNodeCommonAttributes adds to xPath.

        domainMin, domainMax = cls.fixedDomain()
        if domainMin is None: domainMin = float(extraAttributes['domainMin'])
        extraAttributes.pop('domainMin', None)
        if domainMax is None: domainMax = float(extraAttributes['domainMax'])
        extraAttributes.pop('domainMax', None)

        series = cls([], domainMin=domainMin, domainMax=domainMax, **attributes)

        extraNodes = baseModule.XDataFunctional.parseNodeStandardChildren(series, node, xPath, linkData, **kwargs)

        if len(extraNodes) == 1:
            values = extraNodes.pop()
            values = valuesModule.Values.parseNodeUsingClass(values, xPath, linkData, **kwargs)
            series.coefficients = values                                    # FIXME store self.__coefficients as values instance instead of list?

        if len(extraNodes) > 0: raise Exception('Invalid nodes: %s.' % (', '.join([extraNode.tag for extraNode in extraNodes])))

        xPath.pop()

        return series

    @staticmethod
    def fixedDomain( ) :

        return( None, None )

class LegendreSeries( Series1d ) :
    r"""
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

    def __init__(self, coefficients, domainMin = -1, domainMax = 1, lowerIndex = 0, axes = None, index = None, 
            valueType = enumsModule.ValueType.float64, outerDomainValue = None, label = None):

        if( lowerIndex != 0 ) : raise ValueError( 'lowerIndex = %s must be 0' )
        if( ( domainMin != -1 ) or ( domainMax != 1 ) ) : raise ValueError( 'domain must be [-1, 1], not [%s, %s]' % ( domainMin, domainMax ) )
        Series1d.__init__( self, coefficients, -1, 1, axes = axes, index = index, valueType = valueType, outerDomainValue = outerDomainValue, 
                label = label)

    def __iter__( self ) :

        length = len( self )
        for i1 in range( length ) :
            yield self[i1]

    def evaluate( self, mu ) :
        """Using the Legendre coefficients, this method calculates f(mu) and returns it."""

        P = 0.
        for l, c_l in enumerate( self.coefficients ) : P += ( l + 0.5 ) * c_l * Legendre( l, mu, checkXRange = False ) 
        return( P )

    def evaluateBasisFunction(self, mu, l):
        return ( l + 0.5 ) * Legendre( l, mu, checkXRange = False )

    def integrate(self, domainMin=None, domainMax=None):

        if domainMin is None:
            domainMin = -1.0
        if domainMax is None:
            domainMax =  1.0

        sign = 1
        if domainMin > domainMax:
            domainMin, domainMax, sign = domainMax, domainMin, -1
            
        if domainMin < -1.0:
            raise ValueError('domainMin must be greater than or equal to -1: it is %e' % domainMin )
        if domainMax >  1.0:
            raise ValueError('domainMax must be greater than or equal to 1: it is %e' % domainMax )
        if domainMin == domainMax:
            return 0.0
        if domainMin == -1.0 and domainMax == 1.0:
            return self.coefficients[0]

        maxOrder = len(self.coefficients)
        if maxOrder == 0:
            return 0.0

        integral = 0.5 * (domainMax - domainMin) * self.coefficients[0]

        if maxOrder > 1:
            integral += 0.75 * (domainMin + domainMax) * (domainMax - domainMin) * self.coefficients[1]

            P_l_m1_1 = domainMin
            P_l_m1_2 = domainMax
            P_l_1 = 0.5 * (3.0 * domainMin * domainMin - 1.0)
            P_l_2 = 0.5 * (3.0 * domainMax * domainMax - 1.0)
            for order in range(2, maxOrder):
                P_l_p1_1 = Legendre(order + 1, domainMin)
                P_l_p1_2 = Legendre(order + 1, domainMax)
                integral += 0.5 * (P_l_p1_2 - P_l_p1_1 + P_l_m1_1 - P_l_m1_2) * self.coefficients[order]
                P_l_m1_1 = P_l_1
                P_l_m1_2 = P_l_2
                P_l_1 = P_l_p1_1
                P_l_2 = P_l_p1_2
        
        return sign * integral

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

        return( XYs1dModule.XYs1d )

    def toPointwiseLinear( self, **kwargs ) :
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
            for i in range( n ) :
                mu = -1. + ( 2. * i ) / n
                P.append( [ mu, self.evaluate( mu ) ] )
            P.append( [ 1., self.evaluate( 1. ) ] )
        axes = axesModule.Axes(2)
        unit = self.getAxisUnitSafely( 0 )
        axes[0] = axesModule.Axis( 'P(mu)', 0, unit )
        axes[1] = axesModule.Axis( 'mu', 1, '' )
        Pclass = self.toLinearXYsClass()
        P = Pclass( P, axes = axes, outerDomainValue = self.outerDomainValue )
        return( P.thin( accuracy = accuracy ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """Calls self's toPointwiseLinear method."""

        return( self.toPointwiseLinear( **kwargs ) )

    @staticmethod
    def fixedDomain( ) :

        return( -1, 1 )

    @staticmethod
    def fromXYs1d( xys1d, maxOrder ) :

        if( xys1d.domainMin != -1 ) : raise ValueError( "Domain min must be -1. It is %e" % xys1d.domainMin )
        if( xys1d.domainMax !=  1 ) : raise ValueError( "Domain max must be 1. It is %e" % xys1d.domainMax )

        axes = xys1d.axes
        return( LegendreSeries( Legendre_C.from_pointwiseXY_C( xys1d, maxOrder ), axes = axes ) )

class Polynomial1d( Series1d ) :

    moniker = 'polynomial1d'
    dimension = 1

    def __init__(self, coefficients, domainMin, domainMax, lowerIndex = 0, axes = None, index = None, 
            valueType = enumsModule.ValueType.float64, outerDomainValue = None, label = None):

        Series1d.__init__( self, coefficients, domainMin, domainMax, lowerIndex = lowerIndex, axes = axes, index = index, 
                valueType = valueType, outerDomainValue = outerDomainValue, label = label)

    def __findExtrema( self ):

        from numpy.polynomial.polynomial import Polynomial
        p1 = Polynomial(self.coefficients, domain=(self.domainMin, self.domainMax),
            window=(self.domainMin, self.domainMax))
        pprimeroots = p1.deriv().roots()
        return [self.domainMin] + list(pprimeroots) + [self.domainMax]

    def evaluate( self, x ) :
        """Using the polynomial coefficients, this method calculates p(x) and returns it."""

        P = 0.
        for c_l in reversed( self.coefficients ) : P = c_l  + x * P
        return( P )

    def evaluateBasisFunction(self, x, i):
        return pow(x,i)

    @property
    def rangeMin( self ):

        extrema = self.__findExtrema()
        return min([self.evaluate(x) for x in extrema])

    @property
    def rangeMax( self ):

        extrema = self.__findExtrema()
        return max([self.evaluate(x) for x in extrema])

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """
        Sets *domainMin* and *domainMax* per the arguments.
        """

        OldDomainMin = self.domainMin
        OldDomainMax = self.domainMax

        domainMin = max(domainMin, self.domainMin)
        domainMax = min(domainMax, self.domainMax)
        if fixToDomain == enumsModule.FixDomain.lower:
            self.__domainMin = domainMin
        elif fixToDomain == enumsModule.FixDomain.upper:
            self.__domainMax = domainMax
        else:
            self.__domainMin = domainMin
            self.__domainMax = domainMax

        if OldDomainMin == self.domainMin and OldDomainMax == self.domainMax: return 0
        return 1

    def toLinearXYsClass( self ) :

        return( XYs1dModule.XYs1d )

    def toPointwiseLinear( self, **kwargs ) :
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
        for i in range( n + 1 ) :
            x = ( ( n - i ) * self.domainMin + self.domainMax * i ) / n
            P.append( [ x, self.evaluate( x ) ] )
        Pclass = self.toLinearXYsClass()
        P = Pclass( P, axes = self.axes, outerDomainValue = self.outerDomainValue )
        return( P.thin( accuracy = accuracy ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """Calls self' toPointwiseLinear method."""

        return( self.toPointwiseLinear( **kwargs ) )

class LinearSpline1d( Series1d ) :
    r"""
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

    def __init__(self, xdata, ydata, axes, index = None, valueType = enumsModule.ValueType.float64, outerDomainValue= None, label = None):

        if( len( xdata ) != len( ydata ) ) : raise ValueError( "Number of x and y values not equal." )

        Series1d.__init__(self, ydata, xdata[0], xdata[-1], lowerIndex = 0, axes = axes, index = index, valueType = valueType, 
                outerDomainValue = outerDomainValue, label = label)

        self.axes = axes
        self.basis = XYs1dModule.XYs1d(axes=self.axes, data=list(zip(xdata, ydata)), interpolation=enumsModule.Interpolation.linlin)

    def evaluateBasisFunction( self, x, i ) :

        self.basis[i] = (self.basis[i][0], 1.0)
        result=self.basis.evaluate(x)
        self.basis[i] = (self.basis[i][0], 0.0)
        return result

    def toLinearXYsClass( self ) :

        return( XYs1dModule.XYs1d )

    def toPointwiseLinear( self, **kwargs ) :
        """
        This method constructs the y(x) versus x and returns it as a XYs1d instance. Basically we just copy the basis
        function widget and put the y values back.
        """

        result = self.basis.copy()
        for index, coefficient in enumerate( self.coefficients ) :
            result[index] = ( result[index][0], coefficient )

        return result

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """Calls self's toPointwiseLinear method."""

        return( self.toPointwiseLinear( **kwargs ) )

