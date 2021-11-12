# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the classes representing the GNDS documentation nodes author and authors.
"""

from .. import suite as suiteModule
from . import abstractClasses as abstractClassesModule

class Author(abstractClassesModule.AuthorAbstract):

    moniker = 'author'

    @staticmethod
    def parseConstructBareNodeInstance(node, xPath, linkData, **kwargs):

        name = node.get('name')
        orcid = node.get('orcid', '')
        email = node.get('email', '')

        return Author(name, orcid, email)

class Authors(suiteModule.Suite):

    moniker = 'authors'

    def __init__(self):

        suiteModule.Suite.__init__(self, [ Author ])
