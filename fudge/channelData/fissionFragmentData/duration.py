# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import ancestry as ancestryModule

from PoPs.fissionFragmentData import time as timeModule

from fudge import suites as suitesModule

from . import incidentEnergy as incidentEnergyModule
from . import QMatrix as QMatrixModule

class duration( ancestryModule.ancestry ) :

    moniker = 'duration'

    def __init__( self, label ) :

        ancestryModule.ancestry.__init__( self )

        self.__label = label

        self.__time = timeModule.suite( )
        self.__time.setAncestor( self )

        self.__incidentEnergies = incidentEnergyModule.suite( )
        self.__incidentEnergies.setAncestor( self )

        self.__QMatrix = QMatrixModule.suite( )
        self.__QMatrix.setAncestor( self )

    @property
    def label( self ) :

        return( self.__label )

    @property
    def time( self ) :

        return( self.__time )

    @property
    def incidentEnergies( self ) :

        return( self.__incidentEnergies )

    @property
    def QMatrix( self ) :

        return( self.__QMatrix )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s">' % ( indent, self.moniker, self.__label ) ]
        XMLStringList += self.time.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.incidentEnergies.toXMLList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):

        xPath.append(element.tag)
        DUR = cls(element.get('label'))
        for child in element:
            if child.tag == timeModule.suite.moniker:
                DUR.time.parseXMLNode(child, xPath, linkData)
            elif child.tag == incidentEnergyModule.suite.moniker:
                DUR.incidentEnergies.parseXMLNode(child, xPath, linkData)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return DUR

class suite( suitesModule.suite ) :

    moniker = 'durations'

    def __init__( self ) :

        suitesModule.suite.__init__( self, allowedClasses = ( duration, ) )
