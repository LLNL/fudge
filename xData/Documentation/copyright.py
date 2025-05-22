# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the class for the GNDS documentation/copyright node.

This module contains the following classes:

    +---------------------------+-------------------------------------------------------------------------------+
    | Class                     | Description                                                                   |
    +===========================+===============================================================================+
    | Copyright                 | This class represents a GNDS documentation/copyright node.                    |
    +---------------------------+-------------------------------------------------------------------------------+
"""

from .. import text as textModule


class Copyright(textModule.Text):
    """
    A class representing a GNDS documentation/copyright node.

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | href          | A url to the copyright.                                       |
    +---------------+---------------------------------------------------------------+
    """

    moniker = 'copyright'

    def __init__(self, href=''):
        """
        :param href:    URL to the copyright.
        """

        textModule.Text.__init__(self)

        self.href = href

    @property
    def href(self):
        """
        This method returns self's href instance.

        :returns:       A python str.
        """

        return self.__href

    @href.setter
    def href(self, value):
        """
        This method sets the href member of *self* to *value*.

        :param value:       The new href value for *self*.
        """

        if not isinstance(value, str): raise TypeError('href must be a str instance.')
        self.__href = value

    def XML_extraAttributes(self, **kwargs):
        """
        This method returns the XML attributes for *self* as a single python str.

        :kwargs:        This argument is not used.

        :returns:       A python str.
        """

        if self.href == '': return ''

        return ' href="%s"' % self.href

    def parseNode(self, node, xPath, linkData, **kwargs):
        """
        This method fills *self* by parsing the data in *node*.

        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      dictionary of extra arguments controlling how *self* is converted to a list of XML strings.
        """

        textModule.Text.parseNode(self, node, xPath, linkData, **kwargs)
        xPath.append(node.tag)

        self.__href = node.get('href', '')

        xPath.pop()
