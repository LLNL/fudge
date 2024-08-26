# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the standard uncertainty class that is a Gaussian distribution (pdf) about the mean with a specified standard deviation.

This module contains the following classes:

    +---------------------------+-----------------------------------------------------------------------------------------------+
    | Class                     | Description                                                                                   |
    +===========================+===============================================================================================+
    | Standard                  | This is the class for storing an uncertainty as a Gaussian distribution with a specified      |
    |                           | standard deviation.                                                                           |
    +---------------------------+-----------------------------------------------------------------------------------------------+
"""

from . import uncertainty as uncertaintyModule

class Standard( uncertaintyModule.Base ) :

    """
    This class represents an uncertainty as a Gaussian distribution with a specified standard deviation.

    The following table list the primary members of this class:

    +-----------------------+---------------------------------------------------------------------------+
    | Member                | Description                                                               |
    +=======================+===========================================================================+
    | value                 | This is the standard deviation of the Gaussian distribution.              |
    +-----------------------+---------------------------------------------------------------------------+
    """

    moniker = "standard"

    def __init__( self, value ) :
        """
        :param value:       Ths standard deviation of the Gaussian distributed uncertainty.
        """

        uncertaintyModule.Base.__init__( self )

        if not isinstance(value, uncertaintyModule.Quantity): raise TypeError( 'Invalid quantity for %s uncertainty' % self.moniker )
        self.__value = value
        self.__value.setAncestor( self )

    @property
    def value( self ) :
        """
        This method return the *value* member of *self*.

        :returns:       An instance of :py:class:`uncertaintyModule.Quantity`.
        """

        return( self.__value )

    def copy( self ) :
        """
        This method returns a copy of *this*.

        :returns:       An instance of :py:class:`Standard`.
        """

        return( self.__class__( self.__value.copy( ) ) )

    def parentConvertingUnits( self, factors ) :
        """
        This method is call by the parent with the conversion factors used to change units.

        :param factors:     A list of python floats.
        """

        self.__value.parentConvertingUnits( factors )

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XML_strList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XML_strList += self.value.toXML_strList(indent = indent2, **kwargs)
        XML_strList[-1] += "</%s>" % self.moniker

        return XML_strList

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

        xPath.append(node.tag)

        child = node[0]
        if child.tag == uncertaintyModule.Double.moniker:
            uncertainty = uncertaintyModule.Double.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        else :
            raise Exception('Invalid child node = "%s".' % child.tag)

        instance = Standard(uncertainty)

        xPath.pop()

        return instance
