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
from .. import text as textModule
from . import abstractClasses as abstractClassesModule

class ContributorType:

    allowed = ( 'ContactPerson', 'DataCollector', 'DataCurator', 'DataManager', 'Distributor', 'Editor', 'HostingInstitution',
                'Producer', 'ProjectLeader', 'ProjectManager', 'ProjectMember', 'RegistrationAgency', 'RegistrationAuthoriy',
                'RelatedPerson', 'Researcher', 'ResearchGroup', 'RightsHolder', 'Sponsor', 'Supervisor', 'WorkPackageLeader', 'Other' )

class Contributor(abstractClassesModule.AuthorAbstract):

    moniker = 'contributor'
    keyName = 'name'

    def __init__(self, name, orcid, email, contributorType):

        abstractClassesModule.AuthorAbstract.__init__(self, name, orcid, email)

        self.__contributorType = abstractClassesModule.raiseIfNotInList(contributorType, ContributorType.allowed, 'contributorType')

    @property
    def contributorType(self):
        """."""

        return self.__contributorType

    def XML_extraAttributes(self, **kwargs):

        attributes = abstractClassesModule.AuthorAbstract.XML_extraAttributes(self, **kwargs)
        attributes += ' contributorType="%s"' % self.__contributorType
        return attributes

    @staticmethod
    def parseConstructBareNodeInstance(node, xPath, linkData, **kwargs):

        name = node.get('name')
        orcid = node.get('orcid', '')
        email = node.get('email', '')
        contributorType = node.get('contributorType', '')

        return Contributor(name, orcid, email, contributorType)

class Contributors(suiteModule.Suite):

    moniker = 'contributors'

    def __init__(self):

        suiteModule.Suite.__init__(self, [ Contributor ])

    def toXML(self, indent = '', **kwargs):

        return '\n'.join(self.toXMLList(**kwargs))
