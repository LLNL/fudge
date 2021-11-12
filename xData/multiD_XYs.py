# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the XYsnd classes for n > 1. 
"""

__metaclass__ = type

"""
Missing methods

    copyDataToGridWsAndXsAndYs
    def getValue( self, *values ) :
    setFromList
    setFromW_XYs
    thin
    toPointwise_withLinearXYs
    toString
    plot
"""

import abc

from pqu import PQU as PQUModule

from . import formatVersion as formatVersionModule
from . import standards as standardsModule
from . import base as baseModule
from . import axes as axesModule
from . import XYs as XYsModule
from . import regions as regionsModule
from . import series1d as series1dModule
from . import xs_pdf_cdf as xs_pdf_cdfMoudle
from . import uncertainties as uncertaintiesModule

def flatInterpolationToLinearPoint( lowerDomain, upperDomain, domain, epsilon ) :

    if( domain < 0 ) :
        domain *= ( 1.0 - epsilon )
    elif( domain > 0 ) :
        domain *= ( 1.0 + epsilon )
    else :
        domain = epsilon

    if( epsilon != 0.0 ) :
        if( ( domain <= lowerDomain ) or ( domain >= upperDomain ) ) : domain = 0.5 * ( lowerDomain + upperDomain )

    return( domain )
            
class XYsnd( baseModule.xDataFunctional ) :

    ancestryMembers = ( '[functionals', )

    def __init__( self, interpolation = standardsModule.interpolation.linlinToken, axes = None,
            index = None, valueType = standardsModule.types.float64Token, outerDomainValue = None, label = None, 
            interpolationQualifier = standardsModule.interpolation.noneQualifierToken ) :
        """
        Abstract base class constructor for XYsnd class.
        """

        baseModule.xDataFunctional.__init__( self, self.moniker, axes, index = index, valueType = valueType,
                outerDomainValue = outerDomainValue, label = label )

        if( not( isinstance( interpolation, str ) ) ) : raise TypeError( 'interpolation must be a string' )
        self.interpolation = interpolation

        if( not( isinstance( interpolationQualifier, str ) ) ) : raise TypeError( 'interpolation qualifier must be a string' )
        self.interpolationQualifier = interpolationQualifier

        self.__functionals = []

    def __len__( self ) :

        return( len( self.__functionals ) )

    def __getitem__( self, index ) :

        return( self.__functionals[index] )

    def __setitem__( self, index, functional ) :

        index1, functional1 = self._set_insertCommon( index, functional.outerDomainValue, functional )
        if( index1 is not None ) :
            if( ( index1 > 0 ) and ( functional1.outerDomainValue <= self.__functionals[index1-1].outerDomainValue ) ) :
                raise ValueError( 'functional.outerDomainValue = %s is <= prior functional.outerDomainValue = %s' % ( functional1.outerDomainValue , self.__functionals[index1-1].outerDomainValue ) )
            if( ( index1 < ( len( self ) - 1 ) ) and ( functional1.outerDomainValue >= self.__functionals[index1+1].outerDomainValue ) ) :
                raise ValueError( 'functional.outerDomainValue = %s is >= next functional.outerDomainValue = %s' % ( functional1.outerDomainValue, self.__functionals[index1+1].outerDomainValue ) )
            self.__functionals[index1] = functional1

    def  __mul__( self, other ) :

        try :
            value = float( other )
        except :
            raise ValueError( "Other must be a number." )

        functionNd = self.__class__( interpolation = self.interpolation, axes = self.axes.copy( ), index = self.index, outerDomainValue = self.outerDomainValue,
                label = self.label, interpolationQualifier = self.interpolationQualifier )
        for function in self.__functionals : functionNd.append( function * value )

        return( functionNd )

    __rmul__ = __mul__

    @property
    def domainMin( self ) :

        return( self.__functionals[0].outerDomainValue )

    @property
    def domainMax( self ) :

        return( self.__functionals[-1].outerDomainValue )

    @property
    def domainUnit( self ) :

        return( self.getAxisUnitSafely( self.dimension ) )

    @property
    def domainGrid( self ) :

        return( [ functional.outerDomainValue for functional in self ] )

    @property
    def rangeMin( self ) :

        return( min( [ func.rangeMin for func in self ] ) )

    @property
    def rangeMax( self ) :

        return( max( [ func.rangeMax for func in self ] ) )

    @property
    def rangeUnit( self ) :

        return( self.getAxisUnitSafely( 0 ) )

    @property
    def functionals( self ) :
        """Returns self's __functionals."""

        return( self.__functionals )

    @property
    def functionNdsName( self ) :
        """Returns the node name for the child "function#ds"."""

        return( "function%dds" % ( self.dimension - 1 ) )

    def append( self, functional ) :

        self.insert( len( self ), functional )

    def insert( self, index, functional, outerDomainValue = None ) :
        """
        Inserts functional at index. If outerDomainValue is None, outerDomainValue is take from the outerDomainValue of functional.
        """

        if( outerDomainValue is None ) : outerDomainValue = functional.outerDomainValue
        index1, functional1 = self._set_insertCommon( index, outerDomainValue, functional )
        if( index1 is not None ) :
            if( ( index1 > 0 ) and ( outerDomainValue <= self.__functionals[index1-1].outerDomainValue ) ) :
                raise Exception( 'outerDomainValue = %s is <= prior functionals.outerDomainValue = %s' % ( outerDomainValue, self.__functionals[index1-1].outerDomainValue ) )
            if( outerDomainValue >= self.__functionals[index1].outerDomainValue ) :
                raise Exception( 'outerDomainValue = %s is >= next functionals.outerDomainValue = %s. index = %d' % ( outerDomainValue, self.__functionals[index1].outerDomainValue, index1 ) )
            self.__functionals.insert( index1, functional1 )

    def insertAtValue( self, functional, outerDomainValue = None ) :
        """
        Inserts functional at the appropriate index for outerDomainValue. The inserted functional instance will have outerDomainValue outerDomainValue, 
        even if functional as a outerDomainValue.
        """

        if( outerDomainValue is None ) : outerDomainValue = functional.outerDomainValue
        outerDomainValue = float( outerDomainValue )
        index = -1               # Set in case self is empty and next line does not set index or functional.
        for index, functional1 in enumerate( self ) :
            if( functional1.outerDomainValue >= outerDomainValue ) : break
        if( index == -1 ) :
            index = 0
        else :
            if( functional1.outerDomainValue == outerDomainValue ) :
                del self.__functionals[index]
            elif( functional1.outerDomainValue < outerDomainValue ) :     # Happens when outerDomainValue is greater than last items outerDomainValue.
                index += 1
        self.insert( index, functional, outerDomainValue = outerDomainValue )

    def pop( self, index ):

        return self.__functionals.pop( index )

    def _set_insertCommon( self, index, outerDomainValue, functional ) :
        """For internal use only."""

        if( not( isinstance( functional, self.allowedSubElements( ) ) ) ) :
            raise TypeError( 'Invalid class "%s" for insertion into "%s".' % ( functional.__class__, self.__class__ ) )
        outerDomainValue = float( outerDomainValue )
        if( not( isinstance( functional, baseModule.xDataFunctional ) ) ) :
            raise TypeError( 'right-hand-side must be instance of xDataFunctional' )
        if( functional.dimension != ( self.dimension - 1 ) ) :
            raise Exception( 'functional dimension = %d not one less than self diemension = %d'
                             % ( functional.dimension, self.dimension ) )
        n1 = len( self )
        if( n1 < index ) : raise IndexError( 'index = %s while length is %s' % ( index, n1 ) )
        index1 = index
        if( index1 < 0 ) : index1 += n1
        if( index1 < 0 ) : raise IndexError( 'index = %s' % index )
        functional.setAncestor( self )
        if( n1 == 0 ) :
            self.__functionals.append( functional )
            return( None, None )
        elif( n1 == index1 ) :
            if( outerDomainValue <= self.__functionals[-1].outerDomainValue ) :
                raise Exception( 'outerDomainValue = %s is <= prior functional.outerDomainValue = %s' % ( outerDomainValue, self.__functionals[-1].outerDomainValue ) )
            self.__functionals.append( functional )
            return( None, None )
        return( ( index1, functional ) )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        for functional in self : functional.convertUnits( unitMap )
        factors = self.axes.convertUnits( unitMap )
        self.fixValuePerUnitChange( factors )

    def copy( self ) :

        axes = self.axes
        if( axes is not None ) : axes = axes.copy( )
        multid_xys = self.__class__( interpolation = self.interpolation, index = self.index,
                outerDomainValue = self.outerDomainValue, axes = axes, interpolationQualifier = self.interpolationQualifier )
        for i1, functional in enumerate( self ) : multid_xys[i1] = functional.copy( )
        return( multid_xys )

    __copy__ = copy

    def copyDataToNestedLists( self ) :

        return( [ [ subData.outerDomainValue, subData.copyDataToNestedLists( ) ] for subData in self ] )

    def evaluate( self, domainValue, extrapolation = standardsModule.noExtrapolationToken, epsilon = 0, interpolationQualifier = None ) :
        """
        Evaluates the function at the domain point domainValue.
        Interpolation is used if domainValue is between two sub-functions. However, if one of the
        sub-functions is within domainValue * epsilon of domainValue then that sub-function is returned.
        If both sub-functions are within domainValue * epsilon of domainValue, the closest is returned.
        """

        if( interpolationQualifier is None ) : interpolationQualifier = self.interpolationQualifier
        outerDomainValue = baseModule.getDomainValue2( domainValue )

        if( extrapolation not in standardsModule.validExtrapolations ) :
            raise ValueError( 'Invalid extrapolation outerDomainValue = "%s"' % extrapolation )
        position, function1, function2, frac, interpolation, interpolationQualifier2 = self.getBoundingSubFunctions( domainValue )
        if( position is None ) : raise Exception( "No data to interpolate" )

        if( frac <= epsilon ) :             # If close to first point pick it.
            function = function1.copy( )
        elif( abs( 1 - frac ) <= epsilon ) :     # If close to second point pick it.
            function = function2.copy( )
        else :
            if( position in ( '=', '<', '>' ) ) :
                if( position != '=' ) :
                    if( extrapolation != standardsModule.flatExtrapolationToken ) :
                        index = { '<' : 0, '>' : -1 }[position]
                        raise Exception( "evaluation point = %s %s than %s" % 
                                ( outerDomainValue, { '<' : 'less', '>' : 'greater' }[position], self[index].outerDomainValue ) )
                function = function1.copy( )
            else :
                if( not( isinstance( function1, XYsModule.XYs1d ) ) ) :      # FIXME, accuracy, lowerEps and upperEps should not be hardwired.
                    if( hasattr( function1, 'toPointwiseLinear' ) ) :
                        function1 = function1.toPointwiseLinear( accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6 )
                    else :
                        function1 = function1.toPointwise_withLinearXYs( accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6 )
                if( not( isinstance( function2, XYsModule.XYs1d ) ) ) :
                    if( hasattr( function1, 'toPointwiseLinear' ) ) :
                        function2 = function2.toPointwiseLinear( accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6 )
                    else :
                        function2 = function2.toPointwise_withLinearXYs( accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6 )
                if( interpolationQualifier == standardsModule.interpolation.unitBaseToken ) :
                    if( function1.dimension == 1 ) :
                        xy = XYsModule.pointwiseXY_C.unitbaseInterpolate( outerDomainValue, function1.outerDomainValue, function1.nf_pointwiseXY,
                                function2.outerDomainValue, function2.nf_pointwiseXY, 1 )
                    elif( function1.dimension == 2 ) :
                        frac1 = 1.0 - frac
                        frac2 = frac
                        function1 = function1.copy( )
                        function2 = function2.copy( )
                        EPrime1_1 = function1.domainMin
                        EPrime2_1 = function1.domainMax
                        EPrime1_2 = function2.domainMin
                        EPrime2_2 = function2.domainMax
                        EPrime1 = frac1 * EPrime1_1 + frac2 * EPrime1_2
                        EPrime2 = frac1 * EPrime2_1 + frac2 * EPrime2_2

                        energyPrimes = set( )

                        for function1d in function1 :
                            frac = ( EPrime2_1 - function1d.outerDomainValue ) / ( EPrime2_1 - EPrime1_1 )
                            function1d.outerDomainValue = frac * EPrime1 + ( 1.0 - frac ) * EPrime2
                            energyPrimes.add( function1d.outerDomainValue )

                        for function1d in function2 :
                            frac = ( EPrime2_2 - function1d.outerDomainValue ) / ( EPrime2_2 - EPrime1_2 )
                            function1d.outerDomainValue = frac * EPrime1 + ( 1.0 - frac ) * EPrime2
                            energyPrimes.add( function1d.outerDomainValue )

                        energyPrimes = sorted( list( energyPrimes ) )
                        function = function1.__class__( outerDomainValue = outerDomainValue )
                        scale1 = ( EPrime2_1 - EPrime1_1 ) / ( EPrime2 - EPrime1 )
                        scale2 = ( EPrime2_2 - EPrime1_2 ) / ( EPrime2 - EPrime1 )
                        for energyPrime in energyPrimes :
                            function1d1 = function1.evaluate( energyPrime )
                            function1d2 = function2.evaluate( energyPrime )
                            function1d = scale1 * frac1 * function1d1 + scale2 * frac2 * function1d2
                            function.append( function1d )
                        return( function )
                    else :
                        raise ValueError( "Unitbase interpolate of %d dimensional function not supported." % function1.dimension )
                elif( interpolationQualifier == standardsModule.interpolation.unitBaseUnscaledToken ) :
                    xy = XYsModule.pointwiseXY_C.unitbaseInterpolate( outerDomainValue, function1.outerDomainValue, function1.nf_pointwiseXY,
                                                                             function2.outerDomainValue, function2.nf_pointwiseXY, 0 )
                else :
                    xy = ( 1.0 - frac ) * function1 + frac * function2

                try :
                    interpolation = xy.interpolation
                except :
                    interpolation = xy.getInterpolation( )
                function = function1.returnAsClass( function1, xy, outerDomainValue = outerDomainValue, interpolation = interpolation )

        function.outerDomainValue = outerDomainValue

        return( function )

    def findInstancesOfClassInChildren( self, cls, level = 9999 ) :
        """
        Finds all instances of class *cls* in self's children, grand-children, etc.
        """

        foundInstances = []
        level -= 1
        if( level < 0 ) : return( foundInstances )
        for functional in self :
            if( isinstance( functional, cls ) ) : foundInstances.append( functional )
            foundInstances += functional.findInstancesOfClassInChildren( cls, level = level )

        return( foundInstances )

    def integrate( self, **limits ):
        """
        Integrate a XYsnd function. Supports limits for each axis.
        Example:
        >XYsnd.integrate( energy_in = ('1e-5 eV', '10 eV'), energy_out = ('1 keV', '10 keV') )

        :param limits: dictionary containing limits for each independent axis (keyed by axis label or index).
        If an independent axis is missing from the dictionary, integrate over the entire domain of that axis.

        :return: float or PQU
        """

        domainMin, domainMax = None, None
        if( len( limits ) > 0 ) :
            if self.axes[-1].label in limits :
                domainMin, domainMax = limits.pop( self.axes[-1].label )
            elif self.axes[-1].index in limits :
                domainMin, domainMax = limits.pop( self.axes[-1].index )

        xys_ = []
        for functional in self :
            if isinstance( functional, ( XYsModule.XYs1d, series1dModule.series ) ) :
                xys_.append( [ functional.outerDomainValue, functional.integrate( domainMin = domainMin, domainMax = domainMax ) ] )
            elif isinstance( functional, ( XYsnd, regionsModule.regions ) ) :
                xys_.append( [ functional.outerDomainValue, functional.integrate( **limits ) ] )
            else :
                raise TypeError( "Unsupported class for integration: %s" % type( functional ) )
        if( type( xys_[0][1] ) == float ) :
            yUnit = "1"
        else :
            yUnit = xys_[0][1].getUnitSymbol( )
        xys = [ [ x, float( y ) ] for x, y in xys_ ]

        unit = self.getAxisUnitSafely( self.dimension )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        value = float( XYsModule.XYs1d( xys, interpolation = self.interpolation ).integrate( domainMin, domainMax ) )
        if( yUnit == "1" ) :
            return( PQUModule.PQU( value, "" ) )
        else :
            return( PQUModule.PQU( value, baseModule.processUnits( unit, yUnit, '*' ) ) )

    def interpolateAtValue( self, value, unitBase = False, extrapolation = standardsModule.noExtrapolationToken ) :
        """
        Returns a functional with dimension one less than self that is the interpolation of self at value. 
        If value is outside the domain of self and extrapolation is 'noExtrapolationToken' a raise is executed. Otherwise,
        a flat interpolated functional is returned.  If unitBase is True, then unit base interpolation is performed on 
        the lowest dimension and the dependent data.  This method is deprecated (see evaluate).
        """

        if( extrapolation not in standardsModule.validExtrapolations ) : raise ValueError( 'Invalid extrapolation value = "%s"' % extrapolation )
        if( len( self ) == 0 ) : raise Exception( "No data to interpolate" )
        if( value < self[0].outerDomainValue ) :
            if( extrapolation == standardsModule.flatExtrapolationToken ) :
                function = self[0].copy( )
                function.outerDomainValue = value
                return( function )
            else :
                raise Exception( "Interpolation point = %s less than %s" % ( value, self[0].outerDomainValue ) )
        if( value > self[-1].outerDomainValue ) :
            if( extrapolation == standardsModule.flatExtrapolationToken ) :
                function = self[-1].copy( )
                function.outerDomainValue = value
                return( function )
            else :
                raise Exception( "Interpolation point = %s greater than %s" % ( value, self[-1].outerDomainValue ) )
        for index, functional2 in enumerate( self ) :
            if( functional2.outerDomainValue >= value ) : break
        if( value == functional2.outerDomainValue ) :
            function = functional2.copy( )
            function.outerDomainValue = value
            return( function )
        functional1 = self[index-1]
# FIXME: following logic only works if functional1 and functional2 are both XYs1d:
        if( unitBase ) :
            xy = XYsModule.pointwiseXY_C.unitbaseInterpolate( value, functional1.outerDomainValue, functional1.nf_pointwiseXY, 
                    functional2.outerDomainValue, functional2.nf_pointwiseXY, 1 )
        else :
            f = ( functional2.outerDomainValue - value ) / ( functional2.outerDomainValue - functional1.outerDomainValue )
            xy = f * functional1 + ( 1. - f ) * functional2
        xyp = functional1.returnAsClass( functional1, xy, outerDomainValue = value )
        return( xyp )

    def getBoundingSubFunctions( self, value ) :
        """
        Returns the tuple flag, functional1, functional2, frac, interpolation and interpolationQualifier.

        Flag is one of 
            +-------+---------------------------------------------------------------------------+
            | None  | no data,                                                                  |
            +-------+---------------------------------------------------------------------------+
            | '<'   | value below domainMin,                                                    |
            +-------+---------------------------------------------------------------------------+
            | '>'   | value above domainMax,                                                    |
            +-------+---------------------------------------------------------------------------+
            | '='   | value at functional1 or                                                   |
            +-------+---------------------------------------------------------------------------+
            | ''    | functional1.outerDomainValue <= value < functional2.outerDomainValue.     |
            +-------+---------------------------------------------------------------------------+

        If flag is None then functional1, functional2 and frac are also None.  If flag is not '' then functional2 is None.
        """

        interpolation = self.interpolation
        interpolationQualifier = self.interpolationQualifier
        if( len( self ) == 0 ) :
            return( None, None, None, None, interpolation, interpolationQualifier )
        elif( len( self ) == 1 ) :
            if( value == self[0].outerDomainValue ) :
                return( '=', self[0], None, 0.0, interpolation, interpolationQualifier )
            symbol = '<'
            if( value > self[0].outerDomainValue ) : symbol = '>'
            return( symbol, self[0], None, 0.0, interpolation, interpolationQualifier )
        elif( value < self[0].outerDomainValue ) :
            frac = ( self[0].outerDomainValue - value ) / max( abs( value ), abs( self[0].outerDomainValue ) )
            return( '<', self[0], None, frac, interpolation, interpolationQualifier )
        elif( value > self[-1].outerDomainValue ) :
            frac = ( value - self[-1].outerDomainValue ) / max( abs( value ), abs( self[-1].outerDomainValue ) )
            return( '>', self[-1], None, frac, interpolation, interpolationQualifier )

        for index, functional2 in enumerate( self ) :
            if( functional2.outerDomainValue >= value ) : break
            functional1 = functional2
        if( value == functional2.outerDomainValue ) : return( '=', functional2, None, 0, interpolation, interpolationQualifier )
        frac = ( value - functional1.outerDomainValue ) / ( functional2.outerDomainValue - functional1.outerDomainValue )
        return( '', functional1, functional2, frac, interpolation, interpolationQualifier )

    def normalize( self, insitu = True, dimension = None ) :

        selfsDimension = self.dimension
        if( dimension is None ) : dimension = selfsDimension
        if( dimension < 0 ) : dimension += selfsDimension
        if( dimension < 1 ) : raise Exception( 'Dimension %d out of range, must be greater than 1' % dimension )
        multid_xys = self
        if( not( insitu ) ) : multid_xys = self.copy( )
        if( dimension == 0 ) : return( multid_xys )
        if( dimension >= selfsDimension ) :
            multid_xys.scaleDependent( 1. / multid_xys.integrate( ), insitu = True )
        else :
            for functional in multid_xys.__functionals : functional.normalize( insitu = True, dimension = dimension )
        return( multid_xys )

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    def domainSlice( self, domainMin = None, domainMax = None, fill = 1, dullEps = 0. ) :
        """
        Returns a new instance with self sliced between ``domainMin`` and ``domainMax``.

        :param domainMin:   [optional] the lower x-value of the slice, default is domain minimum of self,
        :param domainMax:   [optional] the upper x-value of the slice, default is domain maximum of self,
        :param fill:        [optional] if True, points are added at domainMin and domainMax if they are not in self,
                                       else only existing points in the range [domainMin, domainMax] are included.
        :param dullEps:     [optional] (Currently not implemented) the lower and upper points are dulled, default is 0.
        """

        if( domainMin is None ) : domainMin = self.domainMin
        domainMin = max( domainMin, self.domainMin )
        if( domainMax is None ) : domainMax = self.domainMax
        domainMax = max( domainMax, self.domainMax )

        newMultiD = self.__class__( interpolation = self.interpolation, axes = self.axes.copy( ), index = self.index, valueType = self.valueType, 
                outerDomainValue = self.outerDomainValue, label = self.label, interpolationQualifier = self.interpolationQualifier )

        domainGrid = [ tmp.outerDomainValue for tmp in self ]
        for idx1, val in enumerate( domainGrid ) :
            if( val >= domainMin ) : break
        for idx2, val in enumerate( domainGrid ) :
            if( val >= domainMax ) : break

        if( domainGrid[idx1] == domainMin ) :
            newMultiD.append( self[idx1].copy( ) )
            idx1 += 1
        else :
            newMultiD.append( self.evaluate( domainMin ) )

        for idx in range( idx1, idx2 ) : newMultiD.append( self[idx].copy( ) )

        if( domainGrid[idx2] == domainMax ) :
            newMultiD.append( self[idx2].copy( ) )
        else:
            newMultiD.append( self.evaluate( domainMax ) )

        return( newMultiD )

    def rangeUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.rangeUnit ).getValueAs( unitTo ) )

    def scaleDependent( self, value, insitu = False ) :

        multid_xys = self
        if( not( insitu ) ) : multid_xys = self.copy( )
        for functional in multid_xys : functional.scaleDependent( value, insitu = True )

    def toPointwiseLinear( self, **kwargs ) :

        if( self.interpolation not in [ standardsModule.interpolation.flatToken, standardsModule.interpolation.linlinToken ] ) :
            raise TypeError( 'Unsupported interpolation = "%s".' % self.interpolation )

        flatInterpolation = self.interpolation == standardsModule.interpolation.flatToken
        if( flatInterpolation ) :
            lowerEps = kwargs.get( 'lowerEps', 0.0 )
            upperEps = kwargs.get( 'upperEps', 0.0 )
            lowerEps = abs( lowerEps )
            upperEps = abs( upperEps )
            if( ( lowerEps <= 0.0 ) and ( upperEps <= 0.0 ) ) :
                raise ValueError( 'For "%s" interpolation, one or both of lowerEps and/or upperEps must be greater than zero.' % self.interpolation )
        lowerFlatInterpolation = flatInterpolation and ( lowerEps > 0 )
        upperFlatInterpolation = flatInterpolation and ( upperEps > 0 )

        pointwiseLinear = self.toLinearXYsClass( )
        pointwiseLinear =  pointwiseLinear( interpolation = standardsModule.interpolation.linlinToken, axes = self.axes, index = self.index, 
                valueType = self.valueType, outerDomainValue = self.outerDomainValue, interpolationQualifier = self.interpolationQualifier )

        lastIndex = len( self.__functionals ) - 1
        for index, subFunction in enumerate( self ) :
            if( hasattr( subFunction, 'toPointwiseLinear' ) ) :
                subFunctionPointwiseLinear = subFunction.toPointwiseLinear( **kwargs )
            else :
                subFunctionPointwiseLinear = subFunction.toPointwiseLinear( **kwargs )

            if( flatInterpolation and ( index > 0 ) ) :
                lowerDomainValue = priorSumPointwiseLinear.outerDomainValue
                upperDomainValue = subFunction.outerDomainValue
                if( index == lastIndex ) : lowerEps = 0.0

                priorSumPointwiseLinear.outerDomainValue = flatInterpolationToLinearPoint( lowerDomainValue, upperDomainValue, upperDomainValue, -lowerEps )
                pointwiseLinear.append( priorSumPointwiseLinear )

                if( index < lastIndex ) :
                    nextDomainValue = self[index+1].outerDomainValue
                    subFunctionPointwiseLinear.outerDomainValue = flatInterpolationToLinearPoint( upperDomainValue, nextDomainValue, upperDomainValue, upperEps )
                    pointwiseLinear.append( subFunctionPointwiseLinear )
                    priorSumPointwiseLinear = subFunctionPointwiseLinear.copy( )
            else :
                pointwiseLinear.append( subFunctionPointwiseLinear )
                if( flatInterpolation ) : priorSumPointwiseLinear = subFunctionPointwiseLinear.copy( )

        return( pointwiseLinear )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        arguments = self.getArguments( kwargs, { 'cls' : None } )
        cls = arguments['cls']
        kwargs.pop( 'cls', None )

        if( cls is None ) : cls = self.__class__
        newMultiD = cls( interpolation = self.interpolation, axes = self.axes, index = self.index, valueType = self.valueType, 
                outerDomainValue = self.outerDomainValue, label = self.label, interpolationQualifier = self.interpolationQualifier )
        for subsec in self:
            newPW = subsec.toPointwise_withLinearXYs( cls = subsec.toLinearXYsClass( ), **kwargs )
            newPW.outerDomainValue = subsec.outerDomainValue
            newMultiD.append( newPW )

        return newMultiD

    def toXMLList( self, indent = '', **kwargs ) :

        formatVersion = kwargs.get( 'formatVersion', formatVersionModule.default )

        indent2 = indent  + kwargs.get( 'incrementalIndent', '  ' )
        indent3 = indent2 + kwargs.get( 'incrementalIndent', '  ' )
        if( formatVersion == formatVersionModule.version_1_10 ) : indent3 = indent2

        outline = kwargs.get( 'outline', False )
        if( len( self ) < 6 ) : outline = False

        attributeStr = baseModule.xDataFunctional.attributesToXMLAttributeStr( self )
        if( self.interpolation != standardsModule.interpolation.linlinToken ) :
            attributeStr += ' interpolation="%s"' % self.interpolation
        if( self.interpolationQualifier != standardsModule.interpolation.noneQualifierToken ) :
            attributeStr += ' interpolationQualifier="%s"' % self.interpolationQualifier

        XMLList = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ] 
        if( self.isPrimaryXData( ) ) :
            if( self.axes is not None ) : XMLList += self.axes.toXMLList( indent2 )
        if( 'oneLine' not in kwargs ) :
            if( self.dimension == 2 ) : kwargs['oneLine'] = True

        if( formatVersion != formatVersionModule.version_1_10 ) : XMLList.append( '%s<%s>' % ( indent2, self.functionNdsName ) )
        if( outline ) :
            XMLList += self.__functionals[0].toXMLList( indent3, **kwargs )
            XMLList += self.__functionals[1].toXMLList( indent3, **kwargs )
            XMLList += [ '%s    ... ' % indent3 ]
            XMLList += self.__functionals[-2].toXMLList( indent3, **kwargs )
            XMLList += self.__functionals[-1].toXMLList( indent3, **kwargs )
        else :
            for functional in self.__functionals : XMLList += functional.toXMLList( indent3, **kwargs )
        if( formatVersion != formatVersionModule.version_1_10 ) : XMLList[-1] += "</%s>" % ( self.functionNdsName )
        if( self.uncertainty ) : XMLList += self.uncertainty.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @classmethod
    def parseXMLNode( cls, xDataElement, xPath, linkData, axes = None ) :
        """
        Translates XYsnd XML into the python XYsnd xData class.
        """

        xmlAttr = False
        for attrName in ( 'outerDomainValue', 'label' ) :
            if xDataElement.get(attrName) is not None:
                xmlAttr = True
                xPath.append( '%s[@%s="%s"]' % (xDataElement.tag, attrName, xDataElement.get(attrName) ) )
        if( not xmlAttr ) : xPath.append( xDataElement.tag )

        allowedSubElements = cls.allowedSubElements( )

        attrs = {      'index' : None, 'valueType' : None, 'outerDomainValue' : None,  'label' : None,
                'interpolation' : standardsModule.interpolation.linlinToken, 'interpolationQualifier' : standardsModule.interpolation.noneQualifierToken }
        attributes = { 'index' : int,  'valueType' : str,  'outerDomainValue' : float, 'label' : str,
                'interpolation' : str,                              'interpolationQualifier' : str }
        for key, item in list( xDataElement.items( ) ) :
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key]( item )

        multid_xys = cls( axes = axes, **attrs )

        functionElements = []                   # This support GNDS 1.10 and 2.0
        functions = xDataElement.find( multid_xys.functionNdsName )
        if( functions is not None ) :
            for child in functions : functionElements.append( child )

        for child in xDataElement :
            if( child.tag == axesModule.axes.moniker ) :
                multid_xys.axes = axesModule.axes.parseXMLNode( child, xPath, linkData )
            elif( child.tag == uncertaintiesModule.uncertainty.moniker ) :
                multid_xys.uncertainty = uncertaintiesModule.uncertainty.parseXMLNode( child, xPath, linkData )
            elif( child.tag == multid_xys.functionNdsName ) :
                continue
            else :
                if( functions is not None ) : raise Exception( 'Unsupported child name = "%s".' % child.tag )
                functionElements.append( child ) 

        for child in functionElements :
            subElementClass = None
            for subElement in allowedSubElements :
                if( subElement.moniker == child.tag ) :
                    subElementClass = subElement
                    break
            if( subElementClass is None ) : raise TypeError( 'unknown sub-element "%s" in element "%s"' % ( child.tag, cls.moniker ) )
            xdata = subElementClass.parseXMLNode( child, xPath = xPath, linkData = linkData, axes = multid_xys.axes )
            multid_xys.append( xdata )

        xPath.pop( )
        return( multid_xys )

    @classmethod
    def defaultAxes( cls, labelsUnits ) :
        """
        :param labelsUnits: dictionary of form {
                 0:('dependent label','dependent unit'),
                 1:('1st independent label','1st independent unit'),
                 2:('2nd independent label','2nd independent unit'), ... }
        :return: new axes instance
        """

        return( axesModule.axes( rank = cls.dimension + 1, labelsUnits = labelsUnits ) )

class XYs2d( XYsnd ) :

    moniker = 'XYs2d'
    dimension = 2

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYsModule.XYs1d, series1dModule.series, regionsModule.regions1d, xs_pdf_cdfMoudle.xs_pdf_cdf1d ) )

class XYs3d( XYsnd ) :

    moniker = 'XYs3d'
    dimension = 3

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, regionsModule.regions2d ) )
