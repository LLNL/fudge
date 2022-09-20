# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the 'atomic' class, used to store properties of atoms.
The primary role of atomic data in PoPs is to store excited and partially-ionized atomic configurations,
along with their binding energy and decay radiation.

Each chemicalElement in a PoPs database may contain an <atomic> section.
"""

from LUPY import ancestry as ancestryModule

from . import configuration as configurationModule

class Atomic( ancestryModule.AncestryIO):

    moniker = 'atomic'

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__( self )

        self.__configurations = configurationModule.Configurations( )
        self.__configurations.setAncestor( self )

    @property
    def configurations( self ) :

        return( self.__configurations )

    def convertUnits( self, unitMap ) :
        """See documentation in PoPs.database.convertUnits"""

        self.__configurations.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of self
        """

        _atomic = Atomic( )
        for configuration in self__.configurations : _atomic.configurations.add( configuration.copy( ) )
        return( _atomic )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__configurations.toXML_strList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        instance = cls()

        for child in node:
            if child.tag == configurationModule.Configurations.moniker:
                instance.configurations.parseNode(child, xPath, linkData, **kwargs)
            else :
                raise ValueError('Invalid child = "%s" for %s' % (child.tag, instance.moniker))

        xPath.pop()

        return instance
