# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from .. import suite as suiteModule

"""
This module contains the class for the GNDS documentation child node collaborations.

This module contains the following classes:

    +---------------------------+-----------------------------------------------------------------------------------+
    | Class                     | Description                                                                       |
    +===========================+===================================================================================+
    | Collaboration             | This class represents a GNDS documentation/collaborations/collaboration node.     |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Collaborations            | This is the suite class for the GNDS documentation/collaborations node.           |
    +---------------------------+-----------------------------------------------------------------------------------+
"""

from LUPY import ancestry as ancestryModule

class Collaboration(ancestryModule.AncestryIO):
    """
    This class represent a GNDS documentation/collaborations/collaboration node.

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | name          | The name of the collaboration.                                |
    +---------------+---------------------------------------------------------------+
    | href          | TBD.                                                          |
    +---------------+---------------------------------------------------------------+
    """

    moniker = 'collaboration'
    keyName = 'label'

    def __init__(self, name, href=''):
        """
        :param name:    The nanme of the collaboration.
        :param href:    TBD.
        """

        ancestryModule.AncestryIO.__init__(self)

        self.__name = name
        self.__href = href

    @property
    def name( self ) :
        """
        This method returns self's label.

        :returns:       A python str.
        """

        return( self.__name )

    @property
    def href( self ) :
        """
        This method returns self's href.

        :returns:       A python str.
        """

        return( self.__href )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        href = ''
        if( len( self.__href ) > 0 ) : href = ' href="%s"' % self.__href

        return( [ '%s<%s name="%s"%s/>' % ( indent, self.moniker, self.__name, href ) ] )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :returns:           An instance of *cls* representing *node*.
        """

        name = node.get( 'name' )
        href = node.get( 'href' )

        return cls( name, href)

class Collaborations( suiteModule.Suite ) :
    """
    This is the suite class for the GNDS documentation/collaborations node.
    """

    moniker = 'collaborations'
    suiteName = 'label'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ Collaboration ] )
