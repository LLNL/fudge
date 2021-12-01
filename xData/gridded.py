# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__metaclass__ = type

from . import standards as standardsModule
from . import base as baseModule
from . import axes as axesModule
from . import xDataArray as arrayModule

class gridded( baseModule.xDataFunctional ) :

    ancestryMembers = baseModule.xDataFunctional.ancestryMembers + ( 'array', )

    def __init__( self, axes, array, index = None, valueType = standardsModule.types.float64Token, outerDomainValue = None, label = None ) :

        for axis in axes :
            if( axis.index == 0 ) :
                if( not( isinstance( axis, axesModule.axis ) ) ) : raise Exception( 'dependent axis must not have a grid' )
            else :
                if( not( isinstance( axis, axesModule.grid ) ) ) : raise Exception( 'independent axis must have a grid' )
        if( not( isinstance( array, arrayModule.arrayBase ) ) ) : raise TypeError( 'array not an array instance' )

        baseModule.xDataFunctional.__init__( self, self.moniker, axes, index = index, valueType = valueType,
                outerDomainValue = outerDomainValue, label = label )

        if( self.dimension != len( array ) ) : ValueError( 'self.dimension = %d != len( array ) = %d' % 
                ( self.dimension, len( array ) ) )
        self.array = array
        self.array.setAncestor( self )

    def __len__( self ) :
        """Returns the number of regions in self."""

        return( len( self.array ) )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        if( self.axes is None ) : return
        factors = self.axes.convertUnits( unitMap )
        self.array.offsetScaleValues( 0.0, factors[0] )
        self.fixValuePerUnitChange( factors )

    def copy( self ):

        return self.__class__( self.axes.copy( ), self.array.copy( ), index = self.index, valueType = self.valueType,
            outerDomainValue = self.outerDomainValue, label = self.label )

    __copy__ = copy

    def getGrid( self, index ) :

        if( 0 < index <= len( self ) ) : return( self.axes[index].values )
        raise IndexError( 'index = %d out of range [1,%s]' % (index, len( self ) ) )

    def getGridUnit( self, index ) :

        return( self.axes[index].unit )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = baseModule.xDataFunctional.attributesToXMLAttributeStr( self )
        XMLList = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
        if( self.axes is not None ) : XMLList += self.axes.toXMLList( indent2, **kwargs )
        XMLList += self.array.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData, **kwargs ) :

        xPath.append( element.tag )

        axes = None
        if element.find('axes') is not None:
            axes = axesModule.axes.parseXMLNode( element.find('axes'), xPath, linkData )
        array = arrayModule.arrayBase.parseXMLNode( element[1], xPath, linkData )
        index = int( element.get('index') ) if 'index' in list( element.keys( ) ) else None
        outerDomainValue = float( element.get( 'outerDomainValue' ) ) if 'outerDomainValue' in list( element.keys( ) ) else None
        Grid = cls( axes = axes, array = array, index = index, outerDomainValue = outerDomainValue, label = element.get( 'label' ) )
        xPath.pop( )
        return Grid

class gridded1d( gridded ) :

    moniker = 'gridded1d'
    dimension = 1

class gridded2d( gridded ) :

    moniker = 'gridded2d'
    dimension = 2

class gridded3d( gridded ) :

    moniker = 'gridded3d'
    dimension = 3
