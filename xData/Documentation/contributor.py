# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the classes representing the GNDS documentation nodes author and authors.
"""

from LUPY import enums as enumsModule

from .. import suite as suiteModule
from .. import text as textModule
from . import abstractClasses as abstractClassesModule

class ContributorType(enumsModule.Enum):

    contactPerson = 'ContactPerson'
    dataCollector = 'DataCollector'
    dataCurator = 'DataCurator'
    dataManager = 'DataManager'
    distributor = 'Distributor'
    editor = 'Editor'
    hostingInstitution = 'HostingInstitution'
    producer = 'Producer'
    projectLeader = 'ProjectLeader'
    projectManager = 'ProjectManager'
    projectMember = 'ProjectMember'
    registrationAgency = 'RegistrationAgency'
    registrationAuthoriy = 'RegistrationAuthoriy'
    relatedPerson = 'RelatedPerson'
    researcher = 'Researcher'
    researchGroup = 'ResearchGroup'
    rightsHolder = 'RightsHolder'
    sponsor = 'Sponsor'
    supervisor = 'Supervisor'
    workPackageLeader = 'WorkPackageLeader'
    other  = 'Other'

class Contributor(abstractClassesModule.AuthorAbstract):

    moniker = 'contributor'
    keyName = 'name'

    def __init__(self, name, orcid, email, contributorType):

        abstractClassesModule.AuthorAbstract.__init__(self, name, orcid, email)

        self.__contributorType = ContributorType.checkEnumOrString(contributorType)

    @property
    def contributorType(self):
        """."""

        return self.__contributorType

    def XML_extraAttributes(self, **kwargs):

        attributes = abstractClassesModule.AuthorAbstract.XML_extraAttributes(self, **kwargs)
        attributes += ' contributorType="%s"' % self.__contributorType
        return attributes

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        name = node.get('name')
        orcid = node.get('orcid', '')
        email = node.get('email', '')
        contributorType = node.get('contributorType', '')

        return cls(name, orcid, email, contributorType)

class Contributors(suiteModule.Suite):

    moniker = 'contributors'
    suiteName = 'name'

    def __init__(self):

        suiteModule.Suite.__init__(self, [ Contributor ])

    def toXML(self, indent = '', **kwargs):

        return '\n'.join(self.toXML_strList(**kwargs))
