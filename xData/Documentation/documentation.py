# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the GNDS documentation class.
"""

from .. import ancestry as ancestryModule
from .. import date as dateModule
from . import author as authorModule
from . import contributor as contributorModule
from . import collaboration as collaborationsModule
from . import dates as datesModule
from . import copyright as copyrightModule
from . import acknowledgement as acknowledgementModule
from . import computerCode as computerCodeModule
from . import keyword as keywordModule
from . import relatedItem as relatedItemModule
from . import bibitem as bibitemModule
from . import exforDataSet as exforDataSetModule
from . import texts as textsModule
from .. import text as textParentModule

class Documentation(ancestryModule.Ancestry2):
    """A class representing a GNDS documentation node."""

    moniker = 'documentation'
    ancestryMembers = ( '[authors', '[contributors', '[collaborations', '[dates', 'copyright', '[acknowledgements', '[keywords', 
            '[relatedItems', 'title', 'abstract', 'body', 'computerCodes', 'experimentalDataSets', 'bibliography', 'endfCompatible' )

    def __init__(self, doi='', version='', publicationDate=None):

        ancestryModule.Ancestry2.__init__(self)

        self.doi = doi
        self.version = version
        self.publicationDate = publicationDate

        self.__authors = authorModule.Authors()
        self.__authors.setAncestor(self)

        self.__contributors = contributorModule.Contributors()
        self.__contributors.setAncestor(self)

        self.__collaborations = collaborationsModule.Collaborations()
        self.__collaborations.setAncestor(self)

        self.__dates = datesModule.Dates()
        self.__dates.setAncestor(self)

        self.__copyright = copyrightModule.Copyright()
        self.__copyright.setAncestor(self)

        self.__acknowledgements = acknowledgementModule.Acknowledgements()
        self.__acknowledgements.setAncestor(self)

        self.__keywords = keywordModule.Keywords()
        self.__keywords.setAncestor(self)

        self.__relatedItems = relatedItemModule.RelatedItems()
        self.__relatedItems.setAncestor(self)

        self.__title = textsModule.Title()
        self.__title.setAncestor(self)

        self.__abstract = textsModule.Abstract()
        self.__abstract.setAncestor(self)

        self.__body = textsModule.Body()
        self.__body.setAncestor(self)

        self.__computerCodes = computerCodeModule.ComputerCodes()
        self.__computerCodes.setAncestor(self)

        self.__experimentalDataSets = exforDataSetModule.ExperimentalDataSets()
        self.__experimentalDataSets.setAncestor(self)

        self.__bibliography = bibitemModule.Bibliography()
        self.__bibliography.setAncestor(self)

        self.__endfCompatible = textsModule.EndfCompatible()
        self.__endfCompatible.setAncestor(self)

    @property
    def doi(self):
        """Returns self's doi."""

        return self.__doi

    @doi.setter
    def doi(self, doi):
        """Sets self's doi to *doi*."""

        self.__doi = textParentModule.raiseIfNotString(doi, 'doi')

    @property
    def version(self):
        """Returns self's version."""

        return self.__version

    @version.setter
    def version(self, version):
        """Sets self's version to *version*."""

        self.__version = textParentModule.raiseIfNotString(version, 'version')

    @property
    def publicationDate(self):
        """Returns self's publicationDate."""

        return self.__publicationDate

    @publicationDate.setter
    def publicationDate(self, date):
        """Sets self's publicationDate to *date*."""

        if date is None:
            self.__publicationDate = dateModule.Date(None, dateModule.Resolution.undefined)
        else:
            if not isinstance(date, dateModule.Date): raise TypeError('Invalid date.')
            self.__publicationDate = date

    @property
    def authors(self):

        return self.__authors

    @property
    def contributors(self):

        return self.__contributors

    @property
    def collaborations(self):

        return self.__collaborations

    @property
    def dates(self):

        return self.__dates

    @property
    def copyright(self):

        return self.__copyright

    @property
    def acknowledgements(self):

        return self.__acknowledgements

    @property
    def keywords(self):

        return self.__keywords

    @property
    def relatedItems(self):

        return self.__relatedItems

    @property
    def title(self):

        return self.__title

    @property
    def abstract(self):

        return self.__abstract

    @property
    def body(self):

        return self.__body

    @property
    def computerCodes(self):

        return self.__computerCodes

    @property
    def experimentalDataSets(self):

        return self.__experimentalDataSets

    @property
    def bibliography(self):

        return self.__bibliography

    @property
    def endfCompatible(self):

        return self.__endfCompatible

    def convertUnits(self, unitMap):
        """This method does nothing as there are no units to convert in documentation."""

        return

    def toXMLList(self, indent = '', **kwargs):

        showEmptySuites = kwargs.get('showEmptySuites', False)
        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent

        XMLList  = self.__authors.toXMLList(indent2, **kwargs)
        XMLList += self.__contributors.toXMLList(indent2, **kwargs)
        XMLList += self.__collaborations.toXMLList(indent2, **kwargs)
        XMLList += self.__dates.toXMLList(indent2, **kwargs)
        XMLList += self.__copyright.toXMLList(indent2, **kwargs)
        XMLList += self.__acknowledgements.toXMLList(indent2, **kwargs)
        XMLList += self.__keywords.toXMLList(indent2, **kwargs)
        XMLList += self.__relatedItems.toXMLList(indent2, **kwargs)
        XMLList += self.__title.toXMLList(indent2, **kwargs)
        XMLList += self.__abstract.toXMLList(indent2, **kwargs)
        XMLList += self.__body.toXMLList(indent2, **kwargs)
        XMLList += self.__computerCodes.toXMLList(indent2, **kwargs)
        XMLList += self.__experimentalDataSets.toXMLList(indent2, **kwargs)
        XMLList += self.__bibliography.toXMLList(indent2, **kwargs)
        XMLList += self.__endfCompatible.toXMLList(indent2, **kwargs)

        attributes  = '' if self.doi == '' else ' doi="%s"' % self.doi
        attributes += '' if self.version == '' else ' version="%s"' % self.version
        attributes += self.__publicationDate.asXML_attribute(name = 'publicationDate')

        if not showEmptySuites and len(XMLList) == 0 and len(attributes) == 0: return []

        XMLList.insert(0, '%s<%s%s>' % (indent, self.moniker, attributes))
        XMLList[-1] += '</%s>' % self.moniker

        return XMLList

    def parseNode(self, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        self.parseAncestryMembers(node, xPath, linkData, **kwargs)

        xPath.pop()

    @staticmethod
    def parseConstructBareNodeInstance(node, xPath, linkData, **kwargs):

        doi = node.get('doi', '')
        version = node.get('version', '')
        publicationDate = dateModule.Date.parse(node.get('publicationDate', ''))

        return Documentation(doi, version, publicationDate)
