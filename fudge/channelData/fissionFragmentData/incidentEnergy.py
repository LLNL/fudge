# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import ancestry as ancestryModule

from PoPs.fissionFragmentData import nuclides as nuclidesModule

from fudge import suites as suitesModule
from . import yields as yieldsModule

from PoPs import suite as suitePoPsModule
from PoPs.quantities import quantity as quantityModule

class double( quantityModule.double ) :

    pass

class energy( suitePoPsModule.sortedSuite ) :

    moniker = 'energy'

    def __init__( self ) :

        suitePoPsModule.sortedSuite.__init__( self, allowedClasses = ( double, ) )

class incidentEnergy( ancestryModule.ancestry ) :

    moniker = 'incidentEnergy'

    def __init__( self, label, nuclides ) :

        ancestryModule.ancestry.__init__( self )

        self.__label = label

        self.__energy = energy( )
        self.__energy.setAncestor( self )

        self.__nuclides = nuclidesModule.nuclides( nuclides )
        self.__nuclides.setAncestor( self )

        self.values = None
        self.uncertainty = None

    @property
    def label( self ) :

        return( self.__label )

    @property
    def energy( self ) :

        return( self.__energy )

    @property
    def nuclides( self ) :

        return( self.__nuclides )

    @property
    def values( self ) :

        return( self.__values )

    @values.setter
    def values( self, _values ) :

        if( _values is None ) :
            self.__values = _values
            return

        if( not( isinstance( _values, yieldsModule.values ) ) ) : raise TypeError( 'Invalid values instance.' )
        self.__values = _values
        self.__values.setAncestor( self )

    @property
    def uncertainty( self ) :

        return( self.__uncertainty )

    @uncertainty.setter
    def uncertainty( self, _uncertainty ) :

        if( _uncertainty is None ) :
            self.__uncertainty = _uncertainty
            return

        if( not( isinstance( _uncertainty, yieldsModule.uncertainty ) ) ) : raise TypeError( 'Invalid uncertainty instance.' )
        self.__uncertainty = _uncertainty
        self.__uncertainty.setAncestor( self )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s">' % ( indent, self.moniker, self.__label ) ]
        XMLStringList += self.__energy.toXMLList( indent2, **kwargs )
        XMLStringList += self.__nuclides.toXMLList( indent2, **kwargs )
        if( self.__values is not None ) : XMLStringList += self.__values.toXMLList( indent2, **kwargs )
        if( self.__uncertainty is not None ) : XMLStringList += self.__uncertainty.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):

        xPath.append(element.tag)
        nuclidesList = element.find(nuclidesModule.nuclides.moniker).text.strip().split()
        IE = cls(element.get('label'), nuclidesList)
        for child in element:
            if child.tag == energy.moniker:
                IE.energy.parseXMLNode(child, xPath, linkData)
            elif child.tag == yieldsModule.values.moniker:
                IE.values = yieldsModule.values.parseXMLNode(child, xPath, linkData)
            elif child.tag == yieldsModule.uncertainty.moniker:
                IE.uncertainty = yieldsModule.uncertainty.parseXMLNode(child, xPath, linkData)
            elif child.tag == nuclidesModule.nuclides.moniker:
                pass    # already read
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return IE

class suite( suitesModule.suite ) :

    moniker = 'incidentEnergies'

    def __init__( self ) :

        suitesModule.suite.__init__( self, allowedClasses = ( incidentEnergy, ) )
