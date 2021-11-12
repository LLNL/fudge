# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__metaclass__ = type

import string

from pqu import PQU as PQUModule

from . import ancestry as ancestryModule
from . import link as linkModule
from . import values as valuesModule

noneGridToken = 'none'
pointsGridToken = 'points'
boundariesGridToken = 'boundaries'
parametersGridToken = 'parameters'
linkGridToken = 'link'

normalPDF = 'normal'
lognormalPDF = 'log-normal'

class axis( ancestryModule.ancestry ) :

    moniker = 'axis'
    ancestryMembers = ( '', )

    def __init__( self, label, index, unit ) :
        """
        Returns a new instance of axis.
        """

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label = "%s" is not a string' % label )
        self.__label = label.strip( )

        self.__index = int( index )

        self.unit = unit

    def __str__( self ) :

        return( 'label="%s", index="%s", unit="%s"' % ( self.label, self.__index, self.unit ) )

    def __eq__( self, other ) :

# BRB: why is this method defined?
        return( isinstance( other, axis ) and ( self.label == other.label ) )

    def __ne__( self, other ) :

        return( not( self.__eq__( other ) ) )

    @property
    def keyName( self ) :

        return( 'index' )

    @property
    def keyValue( self ) :

        return( self.index )

    def convertUnits( self, unitMap ) :

        unit, factor = PQUModule.convertUnits( self.unit, unitMap )
        self.unit = unit
        return( factor )

    def copy( self, unresolvedLinks ) :
        """Returns a new instance that is a copy of self."""

        return( axis( self.label, self.index, self.unit ) )

    __copy__ = copy

    @property
    def index( self ) :

        return( self.__index )

    @index.setter
    def index( self, value ) :

        self.__index = value

    @property
    def label( self ) :

        return( self.__label )

    @label.setter
    def label( self, value ) :

        self.__label = value

    @property
    def unit( self ) :

        return( self.__unit )

    @unit.setter
    def unit( self, value ) :
        """Sets self's unit. Only checks that unit is a string. If unit is None, it is set to an empty string (i.e., '')."""

        if( value is None ) : value = ''
        if( not( isinstance( value, str ) ) ) : raise TypeError( 'unit type "%s" is not a string' % type( value ) )
        self.__unit = value.strip( )

    def divideUnit( self, other ) :
        """
        Returns the unit obtained by the division of self.unit by other.unit. Other must be an axis based instance.
        """

        pqu = PQUModule.PQU( 1, self.unit ) / PQUModule.PQU( 1, other.unit )
        return( str( pqu.unit ) )

    def multiplyUnit( self, other ) :
        """
        Returns the unit obtained by the product of self.unit times other.unit. Other must be an axis based instance.
        """

        pqu = PQUModule.PQU( 1, self.unit ) * PQUModule.PQU( 1, other.unit )
        return( str( pqu.unit ) )

    def plotLabel( self ) :

        label = self.label
        if( label == '' ) : label = 'unknown'
        if( self.unit != '' ) : label += ' (%s)' % self.unit
        return( label )

    def toXML( self, indent = '', **kwargs ) :

        XMLStr = '%s<%s index="%d" label="%s" unit="%s"/>' % ( indent, self.moniker, self.index, self.label, self.unit )
        return( XMLStr )

    def toXMLList( self, indent = '', **kwargs ) :

        return( [ self.toXML( indent = indent, **kwargs ) ] )

    def unitConversionFactor( self, newUnit ) :
        """Returns as a float the factor needed to convert self's unit to newUnit. If units are not compatible, a raise is executed."""

        return( PQUModule.PQU( 1., self.unit ).getValueAs( newUnit ) )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( '%s[@index="%s"]' % ( axis.moniker, element.get( 'index' ) ) )

        _axis = axis( element.get( 'label' ), element.get( 'index' ), element.get( 'unit' ) )

        xPath.pop()
        return( _axis )

class grid( axis ) :

    moniker = 'grid'
    ancestryMembers = ( 'values', )

    def __init__( self, label, index, unit, style, values, uncertainty = None, pdf = normalPDF, interpolation = None ) :
        """
        Returns a new instance of grid.
        """

        axis.__init__( self, label, index, unit )

        if( style == linkGridToken ) :
            if( not isinstance( values , linkModule.link ) ):
                raise TypeError( "style = 'link' not consistent with grid '%s'" % values.moniker )
        else :
            if( style not in [ pointsGridToken, boundariesGridToken, parametersGridToken ] ) :
                raise ValueError( 'style = %s not supported' % style )
            if( not( isinstance( values, valuesModule.values ) ) ) : raise TypeError( 'grid not values instance.' )

        self.__style = style
        self.__values = values
        self.values.setAncestor( self )

        self.interpolation = interpolation

            # BRB: uncertainty needs work.
        self.uncertainty = uncertainty
        if( not( isinstance( pdf, str ) ) ) : raise TypeError( 'pdf must be a string' )
        self.pdf = pdf

    @property
    def style( self ) :

        return( self.__style )

    @property
    def values( self ) :

        return( self.__values )

    @property
    def domainMin( self ) :

        return( self.values[0] )

    @property
    def domainMax( self ) :

        return( self.values[-1] )

    @property
    def domainUnit( self ) :

        return( self.unit )

    def domainUnitConversionFactor( self, unitTo ) :

        return( self.unitConversionFactor( unitTo ) )

    @property
    def domainGrid( self ) :

        return( [ value for value in self.values ] )

    def convertToUnit( self, unit ) :

        factor = self.unitConversionFactor( unit )
        self.unit = unit
        if self.style==linkGridToken: return
        self.__values = valuesModule.values( [ factor * value for value in self.values ] )

    def convertUnits( self, unitMap ) :

        factor = axis.convertUnits( self, unitMap )
        if( factor != 1 ) :
            if isinstance( self.__values, linkModule.link ):
                pass
            else:
                self.__values.offsetScaleValues( 0, factor )
        return( factor )

    def copy( self, unresolvedLinks ) :             # FIXME, unresolvedLinks is a kludge until links are handled in a better way.
        """Returns a new grid instance that is a copy of self."""

        _grid = grid( self.label, self.index, self.unit, self.style, self.values.copy( ), uncertainty = self.uncertainty, 
                pdf = self.pdf, interpolation = self.interpolation )
        if( isinstance( self.values, linkModule.link ) ) : unresolvedLinks.append( _grid.values )
        return( _grid )

    __copy__ = copy

    def getIndexOfValue(self,v):
        """
        Get the index of the value in values where x would fit
        :param v:
        :return:
        """
        for ival,val in enumerate(self.values[:-1]):
            if v >= val and v <= self.values[ival+1]: return ival
        return None

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ' style="%s"' % self.style
        if( self.interpolation is not None ) : attributeStr += ' interpolation="%s"' % self.interpolation
        if( self.uncertainty is not None ) :
            attributeStr += ' uncertainty="%s"' % self.uncertainty
            attributeStr += ' pdf="%s"' % self.pdf
        XMLStrList = [ '%s<%s index="%d" label="%s" unit="%s"%s>' % ( indent, self.moniker, self.index, self.label, self.unit, attributeStr ) ]
        XMLStrList += self.values.toXMLList( indent2, **kwargs )
        XMLStrList[-1] += '</%s>' % self.moniker
        return( XMLStrList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( '%s[@index="%s"]' % ( grid.moniker, element.get( 'index' ) ) )

        moniker = None
        href = None
        for key in list( element.keys( ) ) :
            if( 'href' == key[-4:] ) :
                xPath.pop( )
                return( linkModule.link2.parseXMLNode( element, xPath, linkData ) )

        style = element.get( 'style' )
        gridClass = {
            'link': linkModule.link,
            'boundaries': valuesModule.values,
            'parameters': valuesModule.values,
            'points': valuesModule.values,
        }.get( style )
        if( gridClass is None ) : raise Exception( "grid style '%s' not yet supported" % style )

        gridData = gridClass.parseXMLNode( element[0], xPath, linkData )
        _grid = grid( element.get( 'label' ), element.get( 'index' ), element.get( 'unit' ), style, gridData,
                interpolation = element.get( 'interpolation' ) )

        xPath.pop( )
        return( _grid )

class axes( ancestryModule.ancestry ) :

    moniker = 'axes'
    ancestryMembers = ( '[axes', )

    def __init__( self, rank = None, labelsUnits = None ) :
        """
        Constructor for ``axes`` class. For example::

            _axes = axes( labelsUnits = { 0 : ( 'crossSection' , 'b' ), 1 : ( 'energy_in', 'eV' ) } )
        """

        ancestryModule.ancestry.__init__( self )

        if( labelsUnits is None ) : labelsUnits = {}
        if( rank is None ) :
            rank = 2
            if( len( labelsUnits ) > 0 ) : rank = len( labelsUnits )
        rank = int( rank )
        if( not( 0 < rank < 26 ) ) : raise Exception( 'rank = %d must be in the range [1, 25]' % rank )

        self.axes = []
        abcsOffset = string.ascii_lowercase.index( 'y' )
        for index in range( rank ) :
            label, unit = string.ascii_lowercase[abcsOffset-index], ''
            if( index in labelsUnits ) : label, unit = labelsUnits[index]
            self.axes.append( axis( label, index, unit ) )
            self.axes[-1].setAncestor( self, 'index' )

    def __eq__( self, other ) :

        if( isinstance( other, referenceAxes ) ) : return( other.__eq__( self ) )
        if( isinstance( other, axes ) and ( len( self ) == len( other ) ) ) :
            for index, _axis in enumerate( self.axes ) :
                if( _axis != other[index] ) : return( False )
            return( True )
        return( False )

    def __ne__( self, other ) :

        return( not( self.__eq__( other ) ) )

    def __len__( self ) :

        return( len( self.axes ) )

    def __getitem__( self, index ) :

        return( self.axes[index] )

    def __setitem__( self, index, axisOrGrid ) :

        if( not( isinstance( axisOrGrid, ( axis, grid, linkModule.link2 ) ) ) ) : raise TypeError( 'axisOrGrid is not an instance of axis or grid' )
        rank = len( self )
        index = int( index )
        if( index < 0 ) : index += rank 
        if( not( 0 <= index < rank ) ) : raise IndexError( "index = %s out of range for self of rank %s" % ( index, rank ) )
        self.axes[index] = axisOrGrid
        axisOrGrid.index = index
        self.axes[index].setAncestor( self, 'index' )

    def __str__( self ) :

        l = [ str( axis ) for axis in self ]
        return( '\n'.join( l ) )

    def checkRank( self, rank ) :

        if( len( self ) != rank ) : raise Exception( "self's rank = %s != %s" % ( len( self ), rank ) )

    def convertUnits( self, unitMap ) :
        """
        Converts each axis units.
        unitMap is a dictionary of mapping old units to new units (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        factors = []
        for axis in self :
            if( isinstance( axis, linkModule.link2 ) ) : continue
            factors.append( axis.convertUnits( unitMap ) )
        return( factors )

    def copy( self ) :

        unresolvedLinks = []
        newAxes = axes( rank = len( self ) )
        for index, axis in enumerate( self ) : newAxes[index] = axis.copy( unresolvedLinks )
        for object in unresolvedLinks : object.link = object.follow( object )
        return( newAxes )

    __copy__ = copy

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlAxisStringList = []
        for axis in self : xmlAxisStringList.append( axis.toXML( indent = indent2, **kwargs ) )
        XMLList += reversed( xmlAxisStringList )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @staticmethod
    def parseXMLNode( axesElement, xPath, linkData ) :
        """Parse XML element with tag '<axes>'."""

        xPath.append( axesElement.tag )
        if( axesElement.tag == axes.moniker ) :
            _axes = axes( rank = len( axesElement ) )
            for child in axesElement :
                childClass = { axis.moniker : axis, grid.moniker : grid }.get( child.tag )
                if childClass is None:
                    raise TypeError("Unexpected child element '%s' encountered in axes" % child.tag)
                index = child.get( "index" )
                _axes[index] = childClass.parseXMLNode( child, xPath, linkData )
        else :
            raise Exception( 'Invalid tag "%s" for axes' % ( axesElement.tag ) )
        xPath.pop()
        return( _axes )

    @staticmethod
    def parseXMLString( axisString, xPath, linkData ) :

        from xml.etree import cElementTree
        return( axes.parseXMLNode( cElementTree.fromstring( axisString ), xPath = xPath, linkData = linkData ) )

class referenceAxes( ancestryModule.ancestry ) :
    """
    A referenceAxes links to an axes or another referenceAxes instance, although the final link must always be
    an axes instance. All references to a referenceAxes's axis's are de-referenced to the linked axes or referenceAxes
    instance.  A referenceAxes does not write its self to an XML file; but, instead, only reside in Python instances.

    Unlike an axes instance, a referenceAxes does not allow one to change members of the linked axes.
    """

    moniker = 'referenceAxes'
    ancestryMembers = ( '', )

    def __init__( self, axes ) :
        """
        Constructor for ``referenceAxes`` class. For example::

            _axes = referenceAxes( axes )
        """

        ancestryModule.ancestry.__init__( self )

        self.__axes = axes

    def __getitem__( self, index ) :

        return( self.__axes[index] )

    def __str__( self ) :

        return( self.__axes.__str__( ) )

    def __eq__( self, other ) :

        return( self.__axes.__eq__( other ) )

    def __ne__( self, other ) :

        return( self.__axes.__ne__( other ) )

    def __len__( self ) :

        return( len( self.__axes ) )

    def copy( self, unresolvedLinks ) :

        return( referenceAxes( self.__axes ) )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        return( [] )
