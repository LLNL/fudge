# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from . import enums as enumsModule
from . import base as baseModule
from . import axes as axesModule
from . import xDataArray as arrayModule

class Gridded( baseModule.XDataFunctional ) :

    ancestryMembers = baseModule.XDataFunctional.ancestryMembers + ( 'array', )

    def __init__( self, axes, array, index = None, valueType = enumsModule.ValueType.float64, outerDomainValue = None, label = None ) :

        for axis in axes[:self.dimension+1]:
            if( axis.index == 0 ) :
                if( not( isinstance( axis, axesModule.Axis ) ) ) : raise Exception( 'dependent axis must not have a grid' )
            else :
                if( not( isinstance( axis, axesModule.Grid ) ) ) : raise Exception( 'independent axis must have a grid' )
        if( not( isinstance( array, arrayModule.ArrayBase ) ) ) : raise TypeError( 'array not an array instance' )

        baseModule.XDataFunctional.__init__(self, axes, index=index, valueType=valueType, outerDomainValue=outerDomainValue, label=label)

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

    def toXML_strList(self, indent = '', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        attributeStr = baseModule.XDataFunctional.attributesToXMLAttributeStr(self)
        XML_strList = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
        if self.axes is not None: XML_strList += self.axes.toXML_strList(indent2, **kwargs)
        XML_strList += self.array.toXML_strList(indent2, **kwargs)
        XML_strList[-1] += '</%s>' % self.moniker

        return XML_strList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        attributes, extraAttributes = baseModule.XDataFunctional.parseBareNodeCommonAttributes(node, xPath)     # parseBareNodeCommonAttributes adds to xPath.
        if len(extraAttributes) > 0: raise Exception('Invalid attributes: %s.' % ( ', '.join(list(extraAttributes.keys())) ))

        axes = kwargs.get('axes', None)
        if node.find('axes') is not None: axes = axesModule.Axes.parseNodeUsingClass(node.find('axes'), xPath, linkData, **kwargs)

        array = arrayModule.ArrayBase.parseNodeUsingClass(node[1], xPath, linkData, **kwargs)

        instance = cls(axes, array, **attributes)

        xPath.pop()

        return instance

class Gridded1d( Gridded ) :

    moniker = 'gridded1d'
    dimension = 1

    def constructVector(self):
        """
        Generate an xData.vector.Vector object from an xData.gridded.gridded1d instance.

        :returns: instance of xData.vector.Vector.
        """
        return self.array.constructVector()

class Gridded2d( Gridded ) :

    moniker = 'gridded2d'
    dimension = 2

class Gridded3d( Gridded ) :

    moniker = 'gridded3d'
    dimension = 3

    def constructMatrix(self, thirdIndex):
        """
        Generate an xData.matrix.Matrix object from an xData.gridded.gridded3d instance.

        :param thirdIndex: Index along the third dimension from which to extract 2d-array.
        :returns: instance of xData.matrix.Matrix.
        """
        return self.array.constructMatrix(thirdIndex)
