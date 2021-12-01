# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the classes representing the GNDS documentation nodes author and authors.
"""

from .. import ancestry as ancestryModule
from .. import suite as suiteModule
from .. import text as textModule
from . import abstractClasses as abstractClassesModule

class RelationType:

    allowed = ( 'IsCitedBy', 'Cites', 'IsSupplementTo', 'IsSupplementedBy', 'IsContinuedBy', 'Continues', 'Describes', 'IsDescribedBy',
                'HasMetadata', 'IsMetadataFor', 'HasVersion', 'IsVersionOf', 'IsNewVersionOf', 'IsPreviousVersionOf', 'IsPartOf',
                'HasPart', 'IsPublishedIn', 'IsReferencedBy', 'References', 'IsDocumentedBy', 'Documents', 'IsCompiledBy', 'Complies',
                'IsVariantFormOf', 'IsOriginalFormOf', 'IsIdenticalTo', 'IsReviewedBy', 'Reviews', 'IsDerivedFrom', 'IsSourceOf',
                'IsRequiredBy', 'Requires', 'Obsoletes', 'IsObsoletedBy' )

class RelatedItem(ancestryModule.Ancestry2):

    moniker = 'relatedItem'
    keyName = 'name'

    def __init__(self, name, href, relationType):


        ancestryModule.Ancestry2.__init__(self)

        self.__name = textModule.raiseIfNotString(name, 'name')
        self.__href = textModule.raiseIfNotString(href, 'href')
        self.__relationType = abstractClassesModule.raiseIfNotInList(relationType, RelationType.allowed, 'relationType')

    @property
    def name(self):
        """."""

        return self.__name

    @property
    def href(self):
        """."""

        return self.__href

    @property
    def relationType(self):
        """."""

        return self.__relationType

    def toXMLList(self, **kwargs):

        indent = kwargs.get('indent', '')

        attributes = ' name="%s"' % self.__name
        if len(self.__href) > 0: attributes += ' href="%s"' % self.__href
        if len(self.__relationType) > 0: attributes += ' relationType="%s"' % self.__relationType

        return [ '%s<%s%s/>' % ( indent, self.moniker, attributes ) ]

    @staticmethod
    def parseConstructBareNodeInstance(node, xPath, linkData, **kwargs):

        name = node.get('name')
        href = node.get('href', '')
        relationType = node.get('relationType', '')

        return RelatedItem(name, href, relationType)

class RelatedItems(suiteModule.Suite):

    moniker = 'relatedItems'

    def __init__(self):

        suiteModule.Suite.__init__(self, [ RelatedItem ])

    def toXML(self, indent = '', **kwargs):

        return '\n'.join(self.toXMLList(**kwargs))
