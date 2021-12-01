# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the Ys1d class. 
"""

__metaclass__ = type

from pqu import PQU as PQUModule

from . import standards as standardsModule
from . import base as baseModule
from . import axes as axesModule
from . import values as valuesModule
from . import XYs as XYsModule

class Ys1d( baseModule.xDataFunctional ) :

    moniker = 'Ys1d'
    dimension = 1

    ancestryMembers = ( 'Ys', )

    def __init__( self, Ys, interpolation = standardsModule.interpolation.linlinToken, axes = None,
            index = None, valueType = standardsModule.types.float64Token, outerDomainValue = None, label = None ) :

        baseModule.xDataFunctional.__init__( self, self.moniker, axes, index = index, valueType = valueType,
                outerDomainValue = outerDomainValue, label = label )

        if( not( isinstance( interpolation, str ) ) ) : raise TypeError( 'interpolation must be a string' )
        self.interpolation = interpolation

        if( not( isinstance( Ys, valuesModule.values ) ) ) : raise TypeError( 'Ys must be an instance of values.values.' )
        self.__Ys = Ys
        self.__Ys.setAncestor( self )

    def __len__( self ) :

        return( len( self.__Ys ) )

    def __getitem__( self, index ) :

        return( self.__Ys[index] )

    def __add__( self, other ) :

        if( len( self.__Ys ) == 0 ) : return( other.copy( ) )
        if( self.Ys.size != other.Ys.size ) : raise Exception( 'self.Ys.size = %d != other.Ys.size = %d' % ( self.Ys.size, other.Ys.size ) )

        if( self.__Ys.start <= other.Ys.start ) :
            ys1d_1 = self.copy( )
            ys1d_2 = other
        else :
            ys1d_1 = other.copy( )
            ys1d_2 = self 

        offset = ys1d_2.Ys.start - ys1d_1.Ys.start
        values = [ y for y in ys1d_1.Ys.values ]
        for i1, y2 in enumerate( ys1d_2.Ys ) : values[i1+offset] += y2
        ys1d_1.Ys.values = values

        return( ys1d_1 )

    @property
    def Ys( self ) :

        return( self.__Ys )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        factors = self.axes.convertUnits( unitMap )
        yFactor = factors[0]
        if yFactor != 1:
            self.__Ys.values = [ yFactor * value for value in self.__Ys ]
        self.fixValuePerUnitChange( factors )

    def copy( self ) :

        axes = self.axes
        if( axes is not None ) : axes = axes.copy( )
        Ys = self.__class__( self.__Ys.copy( ), interpolation = self.interpolation, axes = axes, 
                index = self.index, outerDomainValue = self.outerDomainValue, label = self.label )
        return( Ys )

    __copy__ = copy

    def evaluate( self, domainValue, extrapolation = standardsModule.noExtrapolationToken, epsilon = 0 ) :

        pass

    @property
    def domainMin( self ) :

        return( self.axes[1].domainMin )

    @property
    def domainMax( self ) :

        return( self.axes[1].domainMax )

    @property
    def domainUnit( self ) :

        return( self.axes[1].domainUnit )

    def domainUnitConversionFactor( self, unitTo ) :

        return( self.axes[1].domainUnitConversionFactor )

    @property
    def domainGrid( self ) :

        return( self.axes[1].domainGrid )

    @property
    def rangeMin( self ) :

        return( min( self.__Ys.values ) )

    @property
    def rangeMax( self ) :

        return( max( self.__Ys.values ) )

    @property
    def rangeUnit( self ) :

        return( self.getAxisUnitSafely( 0 ) )

    def rangeUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.rangeUnit ).getValueAs( unitTo ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        cls = kwargs.pop( 'cls', XYsModule.XYs1d )

        xys = [ self.domainGrid[self.__Ys.start:self.__Ys.end], self.__Ys ]
        return( cls( data = xys, dataForm = 'xsandys', interpolation = self.interpolation, axes = self.axes, index = self.index,
                    valueType = self.valueType, outerDomainValue = self.outerDomainValue, label = self.label ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        outline = kwargs.get( 'outline', False )
        if( len( self ) < 6 ) : outline = False

        attributeStr = baseModule.xDataFunctional.attributesToXMLAttributeStr( self )
        if( self.interpolation != standardsModule.interpolation.linlinToken ) :
            attributeStr += ' interpolation="%s"' % self.interpolation

        XMLList = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ] 
        if( self.isPrimaryXData( ) ) :
            if( self.axes is not None ) : XMLList += self.axes.toXMLList( indent2 )
        XMLList += self.__Ys.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker

        return( XMLList )

    @classmethod
    def parseXMLNode( cls, xDataElement, xPath, linkData, axes = None ) :
        """
        Translates XML Ys1d into a Ys1d instance.
        """

        xmlAttr = False
        for attrName in ( 'label', 'outerDomainValue' ) :
            if xDataElement.get(attrName) is not None:
                xmlAttr = True
                xPath.append( '%s[@%s="%s"]' % (xDataElement.tag, attrName, xDataElement.get(attrName) ) )
                break
        if( not xmlAttr ) : xPath.append( xDataElement.tag )

        attrs = {      'interpolation' : standardsModule.interpolation.linlinToken, 'label' : None, 'index' : None, 'outerDomainValue' : None  }
        attributes = { 'interpolation' : str,                                       'label' : str,  'index' : int,  'outerDomainValue' : float }
        for key, item in list( xDataElement.items( ) ) :
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key]( item )

        axes = None
        for subElement in xDataElement :
            if( subElement.tag == axesModule.axes.moniker ) :
                axes = axesModule.axes.parseXMLNode( subElement, xPath, linkData )
            elif( subElement.tag == valuesModule.values.moniker ) :
                Ys = valuesModule.values.parseXMLNode( subElement, xPath, linkData )
            else :
                raise TypeError( 'sub-element "%s" not valid' % subElement.tag )

        ys1d = cls( Ys, axes = axes, **attrs )

        xPath.pop( )
        return( ys1d )

    @staticmethod
    def defaultAxes( labelsUnits = None ) :
        """
        :param labelsUnits: dictionary of form { 0 : ( 'dependent label',   'dependent unit' ),
                                                 1 : ( 'independent label', 'independent unit' ) }
        :return: new axes instance
        """

        return( axesModule.axes( rank = 2, labelsUnits = labelsUnits ) )
