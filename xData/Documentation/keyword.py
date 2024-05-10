# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the GNDS documentation child nodes keywords and keywork classes.

This module contains the following classes:

    +---------------------------+-----------------------------------------------------------------------------------+
    | Class                     | Description                                                                       |
    +===========================+===================================================================================+
    | Keyword                   | This class represents a GNDS documentation/keywords/keyword node.                 |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Keywords                  | This is the suite class for the GNDS documentation/keywords node.                 |
    +---------------------------+-----------------------------------------------------------------------------------+
"""

from LUPY import ancestry as ancestryModule

from .. import suite as suiteModule
from .. import text as textModule

class Keyword( textModule.Text ) :
    """
    A class representing a GNDS documentation/keywords/keyword node.

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | label         | The unique label for the keyword.                             |
    +---------------+---------------------------------------------------------------+
    | type          | TBD.                                                          |
    +---------------+---------------------------------------------------------------+
    | text          | TBD.                                                          |
    +---------------+---------------------------------------------------------------+
    """

    moniker = 'keyword'
    keyName = 'type'

    def __init__( self, label, type, text ) :
        """
        :param label:   The unique label for the keyword.
        :param type:    TBD.
        :param text:    TBD.
        """

        textModule.Text.__init__( self, text, label = label )

        self.__type = type

    @property
    def type( self ) :
        """
        This method returns the type of self.
        """

        return( self.__type )

    def XML_extraAttributes( self, **kwargs ) :
        """
        This methods returns the XML attributes for *self* as a single python str.

        :kwargs:        This argument is not used.

        :returns:       A python str.
        """

        if( self.__type == '' ) : return ''

        return ' type="%s"' % self.__type

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

        label = node.get( 'label' )
        type = node.get( 'type' )

        return cls(label, type, None)

class Keywords( suiteModule.Suite ) :
    """
    This is the suite class for the GNDS documentation/keywords node.
    """

    moniker = 'keywords'
    suiteName = 'type'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ Keyword ] )
