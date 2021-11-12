# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__metaclass__ = type

from pqu import PQU as PQUModule

from . import ancestry as ancestryModule

from .uncertainty.physicalQuantity import uncertainty as uncertaintyModule

class physicalQuantity( ancestryModule.ancestry ) :

    moniker = 'physicalQuantity'

    def __init__( self, value, unit, label = None ) :

        ancestryModule.ancestry.__init__( self )
        self.__PQ = PQUModule.PQU( value, unit )

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a str instance.' )
        self.__label = label

        self.__uncertainty = None

    def __str__( self ) :

        return( self.toString( ) )

    def __float__( self ) :

        return( float( self.__PQ ) )

    def __eq__( self, other ) :
        """
        If other is not a physicalQuantity instance, False is returned. Otherwise, calls self's compare with epsilonFactor = 5 and test return value.
        """

        if( not isinstance( other, physicalQuantity ) ) : return( False )
        if( not( self.__PQ.isCompatible( other.__PQ.unit ) ) ) : return( False )
        return( self.compare( other, 5 ) == 0 )

    def __ne__( self, other ) :
        """
        If other is not a physicalQuantity instance, False is returned. Otherwise, calls self's compare with epsilonFactor = 5 and test return value.
        """

        if( not isinstance( other, physicalQuantity ) ) : return( True )
        if( not( self.__PQ.isCompatible( other.__PQ.unit ) ) ) : return( True )
        return( self.compare( other, 5 ) != 0 )

    def __lt__( self, other ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a physicalQuantity instance."""

        return( self.compare( other, 5 ) < 0 )

    def __le__( self, other ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a physicalQuantity instance."""

        return( self.compare( other, 5 ) <= 0 )

    def __gt__( self, other ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a physicalQuantity instance."""

        return( self.compare( other, 5 ) > 0 )

    def __ge__( self, other ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a physicalQuantity instance."""

        return( self.compare( other, 5 ) >= 0 )

    def compare( self, other, epsilonFactor = 0 ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a physicalQuantity instance."""

        if( not isinstance( other, physicalQuantity ) ) : raise TypeError( 'other not instance of physicalQuantity' )
        return( self.__PQ.compare( other.__PQ, epsilonFactor ) )

    @property
    def key( self ) :

        return( self.__label )

    @property
    def label( self ) :

        return( self.__label )

    @property
    def value( self ) :

        return( float( self ) )

    @property
    def unit( self ) :

        return( self.__PQ.unit )

    @property
    def uncertainty( self ) :

        return( self.__uncertainty )

    @uncertainty.setter
    def uncertainty( self, _uncertainty ) :

        if( _uncertainty is not None ) :
            if( not( isinstance( _uncertainty, uncertaintyModule.uncertainty ) ) ) : raise TypeError( 'Invalid uncertainty instance.' )

        self.__uncertainty = _uncertainty
        if( self.__uncertainty is not None ) : self.__uncertainty.setAncestor( self )

    def convertToUnit( self, unit ) :

        self.__PQ.convertToUnit( unit )

    def convertUnits( self, unitMap ) :

        unit, factor = PQUModule.convertUnits( self.unit, unitMap )
        self.__PQ.convertToUnit( unit )
        if( self.__uncertainty is not None ) : self.__uncertainty.parentConvertingUnits( [ factor ] )

    def copy( self ) :

        cls = self.__class__( self.value, self.unit, self.label )
        if( self.__uncertainty is not None ) : cls.uncertainty = self.__uncertainty.copy( )
        return( cls )

    __copy__ = copy

    def copyToUnit( self, unit ) :

        _copy = self.copy( )
        _copy.convertToUnit( unit )
        return( _copy )

    def getValueAs( self, unit ) :

        return( self.__PQ.getValueAs( unit ) )

    def toString( self, significantDigits = None, keepPeriod = True ) :

        return( self.__PQ.toString( significantDigits = significantDigits, keepPeriod = keepPeriod ) )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        label = ''
        if( self.label is not None ) : label = ' label="%s"' % self.label

        unit = ''
        if( not( self.unit.isDimensionless( ) ) ) : unit = ' unit="%s"' % self.unit

        ending = '>'
        if( self.__uncertainty is None ) : ending = '/>'
        XMLStringList = [ '%s<%s%s value="%s"%s%s' % ( indent, self.moniker, label, PQUModule.floatToShortestString( self.value, 12 ), unit, ending ) ]

        if( ending == '>' ) :
            if( self.__uncertainty is not None ) : XMLStringList += self.__uncertainty.toXMLList( indent = indent2, **kwargs )
            XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        value = element.get( 'value' )
        unit = element.get( 'unit' )
        label = element.get( 'label', None )
        _cls = cls( value, unit, label )

        for child in element :
            if( child.tag == uncertaintyModule.uncertainty.moniker ) :
                _cls.uncertainty = uncertaintyModule.uncertainty.parseXMLNodeAsClass( child, xPath, linkData )
            else :
                raise ValueError( 'child element with tag "%s" not allowed' % child.tag )

        xPath.pop( ) 
        return( _cls )
