# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the GNDS documentation child node endfCompatible class.

This module contains the following classes:

    +---------------------------+-----------------------------------------------------------------------------------+
    | Class                     | Description                                                                       |
    +===========================+===================================================================================+
    | Acknowledgement           | This class represents a GNDS documentation/acknowledgements/acknowledgement node. |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Acknowledgements          | This is the suite class for the GNDS documentation/acknowledgements node.         |
    +---------------------------+-----------------------------------------------------------------------------------+
"""

from .. import suite as suiteModule
from .. import text as textModule

class Acknowledgement( textModule.Text ) :
    """
    This class represents a GNDS documentation/acknowledgements/acknowledgement node.

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | label         | The unique label for the acknowledgement.                     |
    +---------------+---------------------------------------------------------------+
    """

    moniker = 'acknowledgement'
    keyName = 'label'

    def __init__( self, _label ) :

        textModule.Text.__init__( self )

        self.__label = textModule.raiseIfNotString( _label, 'label' )

    @property
    def label( self ) :
        """
        This method returns the label.

        :returns:   A python str.
        """

        return( self.__label )

    def XML_extraAttributes( self, **kwargs ) :
        """
        This methods returns the XML attributes for *self* as a single python str.

        :kwargs:        This argument is not used.

        :returns:       A python str.
        """

        if( self.__label == '' ) : return ''

        return ' label="%s"' % self.__label

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

        label = node.get( 'label', '' )

        return cls(label)

class Acknowledgements( suiteModule.Suite ) :
    """
    This is the suite class for the GNDS documentation/acknowledgements node.
    """

    moniker = 'acknowledgements'
    suiteName = 'label'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ Acknowledgement ] )
