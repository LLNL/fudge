# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from PoPs.fissionFragmentData import time as timeModule

from fudge import suites as suitesModule

from . import incidentEnergy as incidentEnergyModule
from . import QMatrix as QMatrixModule

class ElapsedTime( ancestryModule.AncestryIO ) :

    moniker = 'elapsedTime'

    def __init__( self, label ) :

        ancestryModule.AncestryIO.__init__( self )

        self.__label = label

        self.__time = timeModule.Suite( )
        self.__time.setAncestor( self )

        self.__incidentEnergies = incidentEnergyModule.Suite( )
        self.__incidentEnergies.setAncestor( self )

        self.__QMatrix = QMatrixModule.Suite( )
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

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s">' % ( indent, self.moniker, self.__label ) ]
        XMLStringList += self.time.toXML_strList( indent = indent2, **kwargs )
        XMLStringList += self.incidentEnergies.toXML_strList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        instance = cls(element.get('label'))
        for child in element:
            if child.tag == timeModule.Suite.moniker:
                instance.time.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == incidentEnergyModule.Suite.moniker:
                instance.incidentEnergies.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return instance

class Suite( suitesModule.Suite ) :

    moniker = 'elapsedTimes'

    def __init__( self ) :

        suitesModule.Suite.__init__( self, allowedClasses = ( ElapsedTime, ) )
