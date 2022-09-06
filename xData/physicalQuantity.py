# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from pqu import PQU as PQUModule

from .uncertainty.physicalQuantity import uncertainty as uncertaintyModule

class PhysicalQuantity( ancestryModule.AncestryIO ) :

    moniker = 'physicalQuantity'

    def __init__( self, value, unit, label = None ) :

        ancestryModule.AncestryIO.__init__( self )
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
        If other is not a PhysicalQuantity instance, False is returned. Otherwise, calls self's compare with epsilonFactor = 5 and test return value.
        """

        if( not isinstance( other, PhysicalQuantity ) ) : return( False )
        if( not( self.__PQ.isCompatible( other.__PQ.unit ) ) ) : return( False )
        return( self.compare( other, 5 ) == 0 )

    def __ne__( self, other ) :
        """
        If other is not a PhysicalQuantity instance, False is returned. Otherwise, calls self's compare with epsilonFactor = 5 and test return value.
        """

        if( not isinstance( other, PhysicalQuantity ) ) : return( True )
        if( not( self.__PQ.isCompatible( other.__PQ.unit ) ) ) : return( True )
        return( self.compare( other, 5 ) != 0 )

    def __lt__( self, other ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a PhysicalQuantity instance."""

        return( self.compare( other, 5 ) < 0 )

    def __le__( self, other ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a PhysicalQuantity instance."""

        return( self.compare( other, 5 ) <= 0 )

    def __gt__( self, other ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a PhysicalQuantity instance."""

        return( self.compare( other, 5 ) > 0 )

    def __ge__( self, other ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a PhysicalQuantity instance."""

        return( self.compare( other, 5 ) >= 0 )

    def compare( self, other, epsilonFactor = 0 ) :
        """Calls self's compare with epsilonFactor = 5 and test return value. Other must be a PhysicalQuantity instance."""

        if( not isinstance( other, PhysicalQuantity ) ) : raise TypeError( 'other not instance of PhysicalQuantity' )
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
            if( not( isinstance( _uncertainty, uncertaintyModule.Uncertainty ) ) ) : raise TypeError( 'Invalid uncertainty instance.' )

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

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        moniker = kwargs.get('moniker', self.moniker)

        label = ''
        if self.label is not None: label = ' label="%s"' % self.label

        unit = ''
        if not self.unit.isDimensionless(): unit = ' unit="%s"' % self.unit

        ending = '>'
        if self.__uncertainty is None: ending = '/>'
        XML_strList = [ '%s<%s%s value="%s"%s%s' % ( indent, moniker, label, PQUModule.floatToShortestString(self.value, 12), unit, ending ) ]

        if ending == '>':
            if self.__uncertainty is not None: XML_strList += self.__uncertainty.toXML_strList(indent = indent2, **kwargs)
            XML_strList[-1] += '</%s>' % moniker

        return XML_strList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        value = node.get('value')
        unit = node.get('unit')
        label = node.get('label', None)
        instance = cls(value, unit, label)

        for child in node:
            if child.tag == uncertaintyModule.Uncertainty.moniker:
                instance.uncertainty = uncertaintyModule.Uncertainty.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else:
                raise ValueError('Invalid child node with tag "%s".' % child.tag)

        xPath.pop() 

        return instance
