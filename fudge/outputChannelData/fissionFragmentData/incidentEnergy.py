# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from PoPs.fissionFragmentData import yields as yieldsModule

from fudge import suites as suitesModule

from PoPs import suite as suitePoPsModule
from PoPs.quantities import quantity as quantityModule

class Double( quantityModule.Double ) :

    pass

class Energy( suitePoPsModule.SortedSuite ) :

    moniker = 'energy'

    def __init__( self ) :

        suitePoPsModule.SortedSuite.__init__( self, allowedClasses = ( Double, ) )

class IncidentEnergy( ancestryModule.AncestryIO ) :

    moniker = 'incidentEnergy'

    def __init__( self, label ) :

        ancestryModule.AncestryIO.__init__( self )

        self.__label = label

        self.__energy = Energy()
        self.__energy.setAncestor(self)

        self.__yields = yieldsModule.Yields()
        self.__yields.setAncestor(self)

    @property
    def label( self ) :

        return( self.__label )

    @property
    def energy( self ) :

        return( self.__energy )

    @property
    def yields(self):

        return self.__yields

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s">' % ( indent, self.moniker, self.__label ) ]
        XMLStringList += self.__energy.toXML_strList( indent2, **kwargs )
        XMLStringList += self.__yields.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)
        IE = cls(node.get('label'))
        for child in node:
            if child.tag == Energy.moniker:
                IE.energy.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == yieldsModule.Yields.moniker:
                IE.yields.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, node.tag))

        xPath.pop()
        return IE

class Suite( suitesModule.Suite ) :

    moniker = 'incidentEnergies'

    def __init__( self ) :

        suitesModule.Suite.__init__( self, allowedClasses = ( IncidentEnergy, ) )
