# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains class for representing a GNDS documentation/bibliography and its entries.

This module contains the following classes:

    +---------------------------+-----------------------------------------------------------------------------------+
    | Class                     | Description                                                                       |
    +===========================+===================================================================================+
    | Bibitem                   | This class represents entries for the GNDS documentation/bibliography node.       |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Bibliography              | This is the suite class for the GNDS documentation/bibliography node.             |
    +---------------------------+-----------------------------------------------------------------------------------+
"""

from .. import text as textModule
from .. import suite as suiteModule

class Bibitem( textModule.Text ) :
    """
    This class represents entries for the GNDS documentation/bibliography node.

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | xref          | TBD.                                                          |
    +---------------+---------------------------------------------------------------+
    | text          | TBD.                                                          |
    +---------------+---------------------------------------------------------------+
    """

    moniker = 'bibitem'
    keyName = 'label'

    def __init__( self, xref, text = None ) :
        """
        :param xref:    TBD.
        :param text:    TBD.
        """

        textModule.Text.__init__( self, text=text)

        self.__xref = textModule.raiseIfNotString( xref, 'xref' )

    @property
    def xref( self ) :
        """
        This method return the xref for self.
        """

        return ( self.__xref )

    def XML_extraAttributes( self, **kwargs ) :
        """
        This methods returns the XML attributes for *self* as a single python str.

        :kwargs:        This argument is not used.

        :returns:       A python str.
        """

        if( self.__xref == '' ) : return ''

        return ' xref="%s"' % self.__xref

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

        xref= node.get('xref')

        return cls(xref)

class Bibliography( suiteModule.Suite ) :
    """
    This is the suite class for the GNDS documentaion/bibliography node.
    """

    moniker = 'bibliography'
    suiteName = 'label'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ Bibitem ] )
