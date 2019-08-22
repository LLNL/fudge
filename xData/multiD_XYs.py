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

import standards as standardsModule
import base as baseModule
import axes as axesModule
import XYs as XYsModule
import regions as regionsModule
import series1d as series1dModule
import uncertainties as uncertaintiesModule
from pqu import PQU

class XYsnd( baseModule.xDataFunctional ) :

    __metaclass__ = abc.ABCMeta
    ancestryMembers = ( '[functionals', )

    def __init__( self, interpolation = standardsModule.interpolation.linlinToken, axes = None,
            index = None, valueType = standardsModule.types.float64Token, value = None, label = None, 
            interpolationQualifier = standardsModule.interpolation.noneQualifierToken ) :
        """
        Abstract base class constructor for XYsnd class.
        """

        baseModule.xDataFunctional.__init__( self, self.moniker, axes, index = index, valueType = valueType,
                value = value, label = label )

        if( not( isinstance( interpolation, str ) ) ) : raise TypeError( 'interpolation must be a string' )
        self.interpolation = interpolation

        if( not( isinstance( interpolationQualifier, str ) ) ) : raise TypeError( 'interpolation qualifier must be a string' )
        self.interpolationQualifier = interpolationQualifier

        self.functionals = []

    def __len__( self ) :

        return( len( self.functionals ) )

    def __getitem__( self, index ) :

        return( self.functionals[index] )

    def __setitem__( self, index, functional ) :

        index_, functional_ = self._set_insertCommon( index, functional.value, functional )
        if( index_ is not None ) :
            if( ( index_ > 0 ) and ( functional_.value <= self.functionals[index_-1].value ) ) :
                raise ValueError( 'functional.value = %s is <= prior functional.value = %s' % ( functional_.value, self.functionals[index_-1].value ) )
            if( ( index_ < ( len( self ) - 1 ) ) and ( functional_.value >= self.functionals[index_+1].value ) ) :
                raise ValueError( 'functional.value = %s is >= next functional.value = %s' % ( functional_.value, self.functionals[index_+1].value ) )
            self.functionals[index_] = functional_

    def append( self, functional ) :

        self.insert( len( self ), functional )

    def insert( self, index, functional, value = None ) :
        """
        Inserts functional at index. If value is None, value is take from the value of functional.
        """

        if( value is None ) : value = functional.value
        index_, functional_ = self._set_insertCommon( index, value, functional )
        if( index_ is not None ) :
            if( ( index_ > 0 ) and ( value <= self.functionals[index_-1].value ) ) :
                raise Exception( 'value = %s is <= prior functionals.value = %s' % ( value, self.functionals[index_-1].value ) )
            if( value >= self.functionals[index_].value ) :
                raise Exception( 'value = %s is >= next functionals.value = %s. index = %d' % ( value, self.functionals[index_].value, index_ ) )
            self.functionals.insert( index_, functional_ )

    def insertAtValue( self, functional, value = None ) :
        """
        Inserts functional at the appropriate index for value. The inserted functional instance will have value value, 
        even if functional as a value.
        """

        if( value is None ) : value = functional.value
        value = float( value )
        index = -1               # Set in case self is empty and next line does not set index or functional.
        for index, functional_ in enumerate( self ) :
            if( functional_.value >= value ) : break
        if( index == -1 ) :
            index = 0
        else :
            if( functional_.value == value ) :
                del self.functionals[index]
            elif( functional_.value < value ) :     # Happens when value is greater than last items value.
                index += 1
        self.insert( index, functional, value = value )

    def pop( self, index ):

        return self.functionals.pop( index )

    def _set_insertCommon( self, index, value, functional ) :
        """For internal use only."""

        if( not( isinstance( functional, self.allowedSubElements( ) ) ) ) :
            raise TypeError( 'Invalid class for insertion: %s' % functional.__class__ )
        value = float( value )
        if( not( isinstance( functional, baseModule.xDataFunctional ) ) ) :
            raise TypeError( 'right-hand-side must be instance of xDataFunctional' )
        if( functional.dimension != ( self.dimension - 1 ) ) :
            raise Exception( 'functional dimension = %d not one less than self diemension = %d'
                             % ( functional.dimension, self.dimension ) )
        n1 = len( self )
        if( n1 < index ) : raise IndexError( 'index = %s while length is %s' % ( index, n1 ) )
        index_ = index
        if( index_ < 0 ) : index_ += n1
        if( index_ < 0 ) : raise IndexError( 'index = %s' % index )
        functional.setAncestor( self )
        if( n1 == 0 ) :
            self.functionals.append( functional )
            return( None, None )
        elif( n1 == index_ ) :
            if( value <= self.functionals[-1].value ) :
                raise Exception( 'value = %s is <= prior functional.value = %s' % ( value, self.functionals[-1].value ) )
            self.functionals.append( functional )
            return( None, None )
        return( ( index_, functional ) )

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
                value = self.value, axes = axes, interpolationQualifier = self.interpolationQualifier )
        for i1, functional in enumerate( self ) : multid_xys[i1] = functional.copy( )
        return( multid_xys )

    __copy__ = copy
    __deepcopy__ = __copy__

    def copyDataToNestedLists( self ) :

        return( [ [ subData.value, subData.copyDataToNestedLists( ) ] for subData in self ] )

    def evaluate( self, domainValue, extrapolation = standardsModule.noExtrapolationToken, epsilon = 0 ) :
        """
        Evaluates the function at the domain point domainValue.
        Interpolation is used if domainValue is between two sub-functions. However, if one of the
        sub-functions is within domainValue * epsilon of domainValue then that sub-function is returned.
        If both sub-functions are within domainValue * epsilon of domainValue, the closes is returned.
        """

        value = baseModule.getDomainValue2( domainValue )

        if( extrapolation not in standardsModule.validExtrapolations ) :
            raise ValueError( 'Invalid extrapolation value = "%s"' % extrapolation )
        position, function1, function2, frac = self.getBoundingSubFunctions( domainValue )
        if( position is None ) : raise Exception( "No data to interpolate" )

        if( frac <= epsilon ) :             # If close to first point pick it.
            function = function1.copy( )
            function.value = value
        if( ( 1 - frac ) <= epsilon ) :     # If close to second point pick it.
            function = function2.copy( )
            function.value = value
        else :
            if( position in ( '=', '<', '>' ) ) :
                if( position != '=' ) :
                    if( extrapolation != standardsModule.flatExtrapolationToken ) :
                        raise Exception( "evaluation point = %s %s than %s" % 
                                ( value, { '<' : 'less', '>' : 'greater' }[position], self[0].value ) )
                function = function1.copy( )
                function.value = value
            else :
                if( not( isinstance( function1, XYsModule.XYs1d ) ) ) :      # FIXME, accuracy, lowerEps and upperEps should not be hardwired.
                    function1 = function1.toPointwise_withLinearXYs( accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6 )
                if( not( isinstance( function2, XYsModule.XYs1d ) ) ) :
                    function2 = function2.toPointwise_withLinearXYs( accuracy = 1e-4, lowerEps = 1e-6, upperEps = 1e-6 )
                if( self.interpolationQualifier == standardsModule.interpolation.unitBaseToken ) :
                    xy = XYsModule.pointwiseXY_C.unitbaseInterpolate( value, function1.value, function1,
                                                                             function2.value, function2, 1 )
                elif( self.interpolationQualifier == standardsModule.interpolation.unitBaseUnscaledToken ) :
                    xy = XYsModule.pointwiseXY_C.unitbaseInterpolate( value, function1.value, function1,
                                                                             function2.value, function2, 0 )
                else :
                    f = ( function2.value - value ) / ( function2.value - function1.value )
                    xy = f * function1 + ( 1. - f ) * function2
                function = function1.returnAsClass( function1, xy, value = value )
        return( function )

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
        if self.axes[-1].label in limits :
            domainMin, domainMax = limits.pop( self.axes[-1].label )
        elif self.axes[-1].index in limits :
            domainMin, domainMax = limits.pop( self.axes[-1].index )

        xys_ = []
        for functional in self :
            if isinstance( functional, ( XYsModule.XYs1d, series1dModule.series ) ) :
                xys_.append( [functional.value, functional.integrate( domainMin = domainMin, domainMax = domainMax ) ] )
            elif isinstance( functional, ( XYsnd, regionsModule.regions ) ) :
                xys_.append( [ functional.value, functional.integrate( **limits ) ] )
            else :
                raise TypeError( "Unsupported class for integration: %s" % type( functional ) )
        yUnit = xys_[0][1].getUnitSymbol( )
        xys = [ [ x, float( y ) ] for x, y in xys_ ]

        unit = self.getAxisUnitSafely( self.dimension )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        value = float( XYsModule.XYs1d( xys, interpolation = self.interpolation ).integrate( domainMin, domainMax ) )
        return( PQU.PQU( value, baseModule.processUnits( unit, yUnit, '*' ) ) )

    def interpolateAtValue( self, value, unitBase = False, extrapolation = standardsModule.noExtrapolationToken ) :
        """
        Returns a functional with dimension one less than self that is the interpolation of self at value. 
        If value is outside the domain of self and extrapolation is 'noExtrapolationToken' a raise is executed. Otherwise,
        a flat interpolated functional is returned.  If unitBase is True, then unit base interpolation is performed on 
        the lowest dimension and the dependent data.  This method is deprecated (see evaluate).
        """

        if( extrapolation not in standardsModule.validExtrapolations ) : raise ValueError( 'Invalid extrapolation value = "%s"' % extrapolation )
        if( len( self ) == 0 ) : raise Exception( "No data to interpolate" )
        if( value < self[0].value ) :
            if( extrapolation == standardsModule.flatExtrapolationToken ) :
                function = self[0].copy( )
                function.value = value
                return( function )
            else :
                raise Exception( "Interpolation point = %s less than %s" % ( value, self[0].value ) )
        if( value > self[-1].value ) :
            if( extrapolation == standardsModule.flatExtrapolationToken ) :
                function = self[-1].copy( )
                function.value = value
                return( function )
            else :
                raise Exception( "Interpolation point = %s greater than %s" % ( value, self[-1].value ) )
        for index, functional2 in enumerate( self ) :
            if( functional2.value >= value ) : break
        if( value == functional2.value ) :
            function = functional2.copy( )
            function.value = value
            return( function )
        functional1 = self[index-1]
# FIXME: following logic only works if functional1 and functional2 are both XYs1d:
        if( unitBase ) :
            xy = XYsModule.pointwiseXY_C.unitbaseInterpolate( value, functional1.value, functional1,
                                                              functional2.value, functional2, 1 )
        else :
            f = ( functional2.value - value ) / ( functional2.value - functional1.value )
            xy = f * functional1 + ( 1. - f ) * functional2
        xyp = functional1.returnAsClass( functional1, xy, value = value )
        return( xyp )

    def getBoundingSubFunctions( self, value ) :

        if( len( self ) == 0 ) : return( None, None, None, None )
        if( value < self[0].value ) :
            frac = ( self[0].value - value ) / max( abs( value ), abs( self[0].value ) )
            return( '<', self[0], None, frac )
        if( value > self[-1].value ) :
            frac = ( value - self[-1].value ) / max( abs( value ), abs( self[-1].value ) )
            return( '>', self[-1], None, frac )
        for index, functional2 in enumerate( self ) :
            if( functional2.value >= value ) : break
            functional1 = functional2
        if( value == functional2.value ) : return( '=', functional2, None, 0 )
        frac = ( value - functional1.value ) / ( functional2.value - functional1.value )
        return( '', functional1, functional2, frac )

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
            for functional in multid_xys.functionals : functional.normalize( insitu = True, dimension = dimension )
        return( multid_xys )

    @property
    def domainMin( self ) :

        return( self.functionals[0].value )

    @property
    def domainMax( self ) :

        return( self.functionals[-1].value )

    @property
    def domainUnit( self ) :

        return( self.getAxisUnitSafely( self.dimension ) )

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQU.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    @property
    def domainGrid( self ) :

        return( [ functional.value for functional in self ] )

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

        newMultiD = self.__class__( interpolation = self.interpolation,
            axes = self.axes.copy( ), index = self.index, valueType = self.valueType, value = self.value,
            label = self.label, interpolationQualifier = self.interpolationQualifier )

        domainGrid = [ tmp.value for tmp in self ]
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

    @property
    def rangeMin( self ) :

        return( min( [ func.rangeMin for func in self ] ) )

    @property
    def rangeMax( self ) :

        return( max( [ func.rangeMax for func in self ] ) )

    @property
    def rangeUnit( self ) :

        return( self.getAxisUnitSafely( 0 ) )

    def rangeUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQU.PQU( '1 ' + self.rangeUnit ).getValueAs( unitTo ) )

    def scaleDependent( self, value, insitu = False ) :

        multid_xys = self
        if( not( insitu ) ) : multid_xys = self.copy( )
        for functional in multid_xys : functional.scaleDependent( value, insitu = True )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        arguments = self.getArguments( kwargs, { 'cls' : None } )
        cls = arguments['cls']
        kwargs.pop( 'cls', None )

        if( cls is None ) : cls = self.__class__
        newMultiD = cls( interpolation = self.interpolation, axes = self.axes, index = self.index,
                    valueType = self.valueType, value = self.value, label = self.label, 
                    interpolationQualifier = self.interpolationQualifier )
        for subsec in self:
            newPW = subsec.toPointwise_withLinearXYs( **kwargs )
            newPW.value = subsec.value
            newMultiD.append( newPW )

        return newMultiD

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
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
        if( outline ) :
            XMLList += self.functionals[0].toXMLList( indent2, **kwargs )
            XMLList += self.functionals[1].toXMLList( indent2, **kwargs )
            XMLList += [ '%s    ... ' % indent2 ]
            XMLList += self.functionals[-2].toXMLList( indent2, **kwargs )
            XMLList += self.functionals[-1].toXMLList( indent2, **kwargs )
        else :
            for functional in self.functionals : XMLList += functional.toXMLList( indent2, **kwargs )
        if( self.uncertainties ) : XMLList += self.uncertainties.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @classmethod
    def parseXMLNode( cls, xDataElement, xPath, linkData, axes = None ) :
        """
        Translates XYsnd XML into the python XYsnd xData class.
        """

        xmlAttr = False
        for attrName in ('value','label'):
            if xDataElement.get(attrName) is not None:
                xmlAttr = True
                xPath.append( '%s[@%s="%s"]' % (xDataElement.tag, attrName, xDataElement.get(attrName) ) )
        if( not xmlAttr ) : xPath.append( xDataElement.tag )

        allowedSubElements = cls.allowedSubElements( )

        attrs = {      'index' : None, 'valueType' : None, 'value' : None,  'label' : None,
                'interpolation' : standardsModule.interpolation.linlinToken, 'interpolationQualifier' : standardsModule.interpolation.noneQualifierToken }
        attributes = { 'index' : int,  'valueType' : str,  'value' : float, 'label' : str,
                'interpolation' : str,                              'interpolationQualifier' : str }
        for key, item in xDataElement.items( ) :
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key]( item )

        for subElement in xDataElement :
            if( subElement.tag == axesModule.axes.moniker ) :
                if( axes is not None ) : raise Exception( 'Multiple "axes" elements present' )
                axes = axesModule.axes.parseXMLNode( subElement, xPath, linkData )

        multid_xys = cls( axes = axes, **attrs )
        uncertainties = None
        for child in xDataElement :
            if( child.tag == axesModule.axes.moniker ) :
                continue
            elif( child.tag == uncertaintiesModule.uncertainties.moniker ) :
                uncertainties = uncertaintiesModule.uncertainties.parseXMLNode( child, xPath, linkData )
                continue
            else :
                subElementClass = None
                for subElement in allowedSubElements :
                    if( subElement.moniker == child.tag ) :
                        subElementClass = subElement
                        break
                if( subElementClass is None ) : raise TypeError( 'unknown sub-element "%s" in element "%s"' % ( child.tag, cls.moniker ) )
                xdata = subElementClass.parseXMLNode( child, xPath = xPath, linkData = linkData, axes = axes )
                multid_xys.append( xdata )
        if uncertainties is not None: multid_xys.uncertainties = uncertainties

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

        return( ( XYsModule.XYs1d, series1dModule.series, regionsModule.regions1d ) )

class XYs3d( XYsnd ) :

    moniker = 'XYs3d'
    dimension = 3

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs2d, regionsModule.regions2d ) )
