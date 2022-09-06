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
from LUPY import ancestry as ancestryModule

from .. import suite as suiteModule
from .. import text as textModule
from . import abstractClasses as abstractClassesModule

class RelationType(enumsModule.Enum):

    none = enumsModule.auto()
    isCitedBy = 'IsCitedBy'
    cites = 'Cites'
    isSupplementTo = 'IsSupplementTo'
    isSupplementedBy = 'IsSupplementedBy'
    isContinuedBy = 'IsContinuedBy'
    continues = 'Continues'
    describes = 'Describes'
    isDescribedBy = 'IsDescribedBy'
    hasMetadata = 'HasMetadata'
    isMetadataFor = 'IsMetadataFor'
    hasVersion = 'HasVersion'
    isVersionOf = 'IsVersionOf'
    isNewVersionOf = 'IsNewVersionOf'
    isPreviousVersionOf = 'IsPreviousVersionOf'
    isPartOf = 'IsPartOf'
    hasPart = 'HasPart'
    isPublishedIn = 'IsPublishedIn'
    isReferencedBy = 'IsReferencedBy'
    references = 'References'
    isDocumentedBy = 'IsDocumentedBy'
    documents = 'Documents'
    isCompiledBy = 'IsCompiledBy'
    complies = 'Complies'
    isVariantFormOf = 'IsVariantFormOf'
    isOriginalFormOf = 'IsOriginalFormOf'
    isIdenticalTo = 'IsIdenticalTo'
    isReviewedBy = 'IsReviewedBy'
    reviews = 'Reviews'
    isDerivedFrom = 'IsDerivedFrom'
    isSourceOf = 'IsSourceOf'
    isRequiredBy = 'IsRequiredBy'
    requires = 'Requires'
    obsoletes = 'Obsoletes'
    isObsoletedBy = 'IsObsoletedBy'

class RelatedItem(ancestryModule.AncestryIO):

    moniker = 'relatedItem'
    keyName = 'name'

    def __init__(self, name, href, relationType):

        ancestryModule.AncestryIO.__init__(self)

        self.__name = textModule.raiseIfNotString(name, 'name')
        self.__href = textModule.raiseIfNotString(href, 'href')
        self.__relationType = RelationType.checkEnumOrString(relationType)

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

    def toXML_strList(self, **kwargs):

        indent = kwargs.get('indent', '')

        attributes = ' name="%s"' % self.__name
        if len(self.__href) > 0: attributes += ' href="%s"' % self.__href
        if len(self.__relationType) > 0: attributes += ' relationType="%s"' % self.__relationType

        return [ '%s<%s%s/>' % ( indent, self.moniker, attributes ) ]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        name = node.get('name')
        href = node.get('href', '')
        relationType = node.get('relationType', '')

        return cls(name, href, relationType)

class RelatedItems(suiteModule.Suite):

    moniker = 'relatedItems'
    suiteName = 'name'

    def __init__(self):

        suiteModule.Suite.__init__(self, [ RelatedItem ])

    def toXML(self, indent = '', **kwargs):

        return '\n'.join(self.toXML_strList(**kwargs))
