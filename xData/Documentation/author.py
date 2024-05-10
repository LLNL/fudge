# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the classes representing the GNDS documentation nodes author and authors.

This module contains the following classes:

    +---------------------------+-----------------------------------------------------------------------------------+
    | Class                     | Description                                                                       |
    +===========================+===================================================================================+
    | Author                    | This class represents a GNDS documentation/authors/author node.                   |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Authors                   | This is the suite class for the GNDS documentation/authors node.                  |
    +---------------------------+-----------------------------------------------------------------------------------+
"""

from .. import suite as suiteModule
from . import abstractClasses as abstractClassesModule

class Author(abstractClassesModule.AuthorAbstract):
    """
    This class represents a GNDS documentation/authors/author node.
    """

    moniker = 'author'
    keyName = 'name'

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

        name = node.get('name')
        orcid = node.get('orcid', '')
        email = node.get('email', '')

        return cls(name, orcid, email)

class Authors(suiteModule.Suite):
    """
    This is the suite class for the GNDS documentation/authors node.
    """

    moniker = 'authors'
    suiteName = 'name'

    def __init__(self):

        suiteModule.Suite.__init__(self, [ Author ])
