# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__metaclass__ = type

from numericalFunctions import pointwiseXY_C as pointwiseXY_CModule
floatToShortestString = pointwiseXY_CModule.floatToShortestString

from . import base as baseModule
from . import axes as axesModule

class constant( baseModule.xDataFunctional )  :

    def __init__( self, _value, domainMin, domainMax, axes = None, label = None ) :

        baseModule.xDataFunctional.__init__( self, self.moniker, label = label, axes = axes )

        self.value = _value

        if( isinstance( domainMin, int ) ) : domainMin = float( domainMin )
        if( not( isinstance( domainMin, float ) ) ) : TypeError( 'domainMin not a float instance' )
        self.__domainMin = domainMin

        if( isinstance( domainMax, int ) ) : domainMax = float( domainMax )
        if( not( isinstance( domainMax, float ) ) ) : TypeError( 'domainMax not a float instance' )
        self.__domainMax = domainMax

    def copy( self ) :

        axes = self.axes
        if( axes is not None ) : axes = self.axes.copy( )
        return( self.__class__( self.value, self.domainMin, self.domainMax, axes = axes, label = self.label ) )

    __copy__ = copy

    @property
    def value( self ) :

        return( self.__value )

    @value.setter
    def value( self, _value ) :

        if( isinstance( _value, int ) ) : _value = float( _value )
        if( not( isinstance( _value, float ) ) ) : TypeError( 'value not a float instance' )
        self.__value = _value

    @property
    def domainMin( self ) :

        return( self.__domainMin )

    @property
    def domainMax( self ) :

        return( self.__domainMax )

    @property
    def domainUnit( self ) :

        return( self.getAxisUnitSafely( self.dimension ) )

    @property
    def rangeMin( self ) :

        return( self.__value )

    rangeMax = rangeMin

    @property
    def rangeUnit( self ) :

        return( self.getAxisUnitSafely( 0 ) )

    def fixDomainPerUnitChange( self, factors ) :

        self.__domainMin *= factors[self.dimension]
        self.__domainMax *= factors[self.dimension]

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        valueFormatter = kwargs.get( 'valueFormatter', floatToShortestString )
        significantDigits = kwargs.get( 'significantDigits', 15 )

        attributeStr = baseModule.xDataCoreMembers.attributesToXMLAttributeStr( self )
        XMLList = [ '%s<%s%s value="%s" domainMin="%s" domainMax="%s">' % ( indent, self.moniker, attributeStr, 
                valueFormatter( self.__value, significantDigits = significantDigits ),
                valueFormatter( self.domainMin, significantDigits = significantDigits ),
                valueFormatter( self.domainMax, significantDigits = significantDigits ) ) ]

        if( self.isPrimaryXData( ) ) :
            if( self.axes is not None ) : XMLList += self.axes.toXMLList( indent = indent2, **kwargs )

        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @classmethod
    def parseXMLNode( cls, xDataElement, xPath, linkData, axes = None, **kwargs ) :

        attrs =      { 'value' :  None, 'domainMin' :  None, 'domainMax' :  None, 'label' : None }
        attributes = { 'value' : float, 'domainMin' : float, 'domainMax' : float, 'label' : str }
        for key, item in list( xDataElement.items( ) ) :
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key]( item )

        uncertainty = None
        for subElement in xDataElement :
            if( subElement.tag == 'axes' ) :
                axes = axesModule.axes.parseXMLNode( subElement, xPath, linkData )
            elif( subElement.tag == 'uncertainty' ) :
                from . import uncertainties as uncertaintiesModule
                uncertainty = uncertaintiesModule.uncertainty.parseXMLNode( subElement, xPath, linkData )
            else :
                raise TypeError( 'sub-element "%s" not valid' % subElement.tag )

        _value = attrs.pop( 'value' )
        newConstant = cls( _value, axes = axes, **attrs )
        newConstant.uncertainty = uncertainty
        return newConstant

    @staticmethod
    def parseXMLString( XMLString ) :

        from xml.etree import cElementTree

        return( constant.parseXMLNode( cElementTree.fromstring( XMLString ), xPath = [], linkData = {} ) )

class constant1d( constant ) :

    moniker = 'constant1d'
    dimension = 1

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        factors = self.axes.convertUnits( unitMap )
        self.value *= factors[0]
        self.fixDomainPerUnitChange( factors )
        self.fixValuePerUnitChange( factors )

    def evaluate( self, x ) :

        return( self.value )
