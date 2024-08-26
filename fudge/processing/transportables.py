# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains classes for storing the data in a transportable node.

This module contains the following classes:

    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Class                                 | Description                                                                       |
    +=======================================+===================================================================================+
    | Transportable                         | This class stores a transportable instanse for a particle.                        |
    +---------------------------------------+-----------------------------------------------------------------------------------+
    | Transportables                        | This class stores a list of :py:class:`Transportable` instances.                  |
    +---------------------------------------+-----------------------------------------------------------------------------------+
"""

from LUPY import ancestry as ancestryModule

from fudge import enums as enumsModule
from fudge import suites as suitesModule

from . import group as groupModule

class Transportable(ancestryModule.AncestryIO):
    """
    This class stores a product's conservation flag and its multi-group boundaries.
    The conversation flag must be a :py:class:`enumsModule.Conserve` instance.
    The conversation flag indicates the weight used for the outgoing energy integral 
    when calculating the product's transfer matrices.
    The multi-group boundaries are stored in a groupModule.Group instance.

    +---------------+-------------------------------------------------------------------+
    | Member        | Description                                                       |
    +===============+===================================================================+
    | particle      | The GNDS PoPs id of the product whose data are defined.           |
    +---------------+-------------------------------------------------------------------+
    | conserve      | The conversation flag for the product's multi-group data.         |
    +---------------+-------------------------------------------------------------------+
    | group         | The multi-group boundaries.                                       |
    +---------------+-------------------------------------------------------------------+
    """

    moniker = 'transportable'

    def __init__(self, particle, conserve, group):
        """
        :param particle:    The GNDS PoPs id of the product whose data are defined.
        :param conserve:    The conversation flag for the product's multi-group data. 
        :param group:       The multi-group boundaries.
        """

        ancestryModule.AncestryIO.__init__(self)

        if( not( isinstance( particle, str ) ) ) : raise TypeError( 'particle must only be a string instance.' )
        self.__particle = particle

        self.__conserve = enumsModule.Conserve.checkEnumOrString(conserve)

        if( not( isinstance( group, groupModule.Group ) ) ) : raise TypeError( 'Group boundaries must only be a grid instance.' )
        self.__group = group

    @property
    def label( self ) :
        """This function returns the paritlce's GNDS PoPs id."""

        return( self.particle )

    @property
    def particle( self ) :
        """This function returns the paritlce's GNDS PoPs id."""

        return( self.__particle )

    @property
    def conserve( self ) :
        """This function returns the paritlce's conversation flag."""

        return( self.__conserve )

    @property
    def group( self ) :
        """This function returns a reference to the multi-group boundaries."""

        return( self.__group )

    def convertUnits(self, unitMap):
        """
        Converts unit per *unitMap*.

        :param unitMap:         A python dictionary with the keys being the current units and the values being the new units.
        """

        self.__group.convertUnits(unitMap)

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a python list of str instances representing the XML lines of *self*.
    
        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings. 
    
        :return:                Python list of str instances.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s label="%s" conserve="%s">' % ( indent, self.moniker, self.label, self.conserve ) ]
        xmlStringList += self.group.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    Dictionary that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *node*.
        """

        xPath.append('%s[@label="%s"]' % (node.tag, node.get('label')))

        label = node.get('label')
        group = groupModule.Group.parseNodeUsingClass(node.find(groupModule.Group.moniker), xPath, linkData, **kwargs)
        conserve = node.get('conserve')

        xPath.pop()

        return cls( label, conserve, group )

class Transportables( suitesModule.Suite ) :
    """
    This class stores a list of :py:class:`Transportable` instances.
    """

    moniker = 'transportables'

    def __init__( self ) :

        suitesModule.Suite.__init__( self, ( Transportable, ), allow_href = True )

    def parseNode(self, node, xPath, linkData, **kwargs):
        """
        This is a temporary kludge and should be in suitesModule.Suite.parseNode.

        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    Dictionary that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        """
# FIXME

        href = node.get('href')
        if href is None:
            suitesModule.Suite.parseNode(self, node, xPath, linkData, **kwargs)
        else:
            self.set_href(href)
