# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the Affiliation, Affiliations, Note and AuthorAbstract classes.

This module contains the following classes:

    +---------------------------+-------------------------------------------------------------------------------+
    | Class                     | Description                                                                   |
    +===========================+===============================================================================+
    | Affiliation               | This class represents a GNDS author/affiliations/affiliation node.            |
    +---------------------------+-------------------------------------------------------------------------------+
    | Affiliations              | This is the suite class for the GNDS/documentation/authors affiliations node. |
    +---------------------------+-------------------------------------------------------------------------------+
    | Note                      | This is a class representing a GNDS authors/author/note node.                 |
    +---------------------------+-------------------------------------------------------------------------------+
    | AuthorAbstract            | TThis is an abstract class for GNDS/documentation author nodes.               |
    +---------------------------+-------------------------------------------------------------------------------+
"""
from abc import ABC

from LUPY import ancestry as ancestryModule

from .. import suite as suiteModule
from .. import text as textModule


class Affiliation(textModule.Text):
    """
    This class represents a GNDS author/affiliations/affiliation node.

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | name          | The name of the institution.                                  |
    +---------------+---------------------------------------------------------------+
    | href          | The url for the institution.                                  |
    +---------------+---------------------------------------------------------------+
    """

    moniker = 'affiliation'
    keyName = 'label'

    def __init__(self, name='', href=''):
        """
        :param name:    The name of the institution.
        :param href:    The url for the institution.
        """

        textModule.Text.__init__(self)

        self.name = name
        self.href = href

    @property
    def name(self):
        """Returns self's name instance."""

        return self.__name

    @name.setter
    def name(self, value):

        self.__name = textModule.raiseIfNotString(value, 'name')

    @property
    def href(self):
        """Returns self's href instance."""

        return self.__href

    @href.setter
    def href(self, value):

        self.__href = textModule.raiseIfNotString(value, 'href')

    def XML_extraAttributes(self, **kwargs):
        """
        This method returns the XML attributes for *self* as a single python str.

        :kwargs:        This argument is not used.

        :returns:       A python str.
        """

        attributes = ''
        if len(self.__name) > 0: attributes += ' name="%s"' % self.__name
        if len(self.__href) > 0: attributes += ' href="%s"' % self.__href

        return attributes

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        name = node.get('name')
        href = node.get('href')

        return cls(name, href)


class Affiliations(suiteModule.Suite):
    """
    This is the suite class for the GNDS/documentation/authors affiliations node.
    """

    moniker = 'affiliations'
    suiteName = 'label'

    def __init__(self):
        suiteModule.Suite.__init__(self, [Affiliation])


class Note(textModule.Text):
    """
    This is a class representing a GNDS authors/author/note node.
    """

    moniker = 'note'


class AuthorAbstract(ancestryModule.AncestryIO, ABC):
    """
    This is an abstract class for GNDS/documentation author nodes.

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | name          | The author's name.                                            |
    +---------------+---------------------------------------------------------------+
    | orcid         | The orcid of the author.                                      |
    +---------------+---------------------------------------------------------------+
    | email         | The email address of the author.                              |
    +---------------+---------------------------------------------------------------+
    | affiliations  | A suite of auther affiliations.                               |
    +---------------+---------------------------------------------------------------+
    | note          | Notes about the author.                                       |
    +---------------+---------------------------------------------------------------+
    """

    ancestryMembers = ('affiliations', 'note')

    def __init__(self, name, orcid, email):
        """
        :param name:    The author's name.
        :param orcid:   The orcid of the author.
        :param email:   The email address of the author.
        """

        ancestryModule.AncestryIO.__init__(self)

        self.__name = textModule.raiseIfNotString(name, 'name')
        self.__orcid = textModule.raiseIfNotString(orcid, 'orcid')
        self.__email = textModule.raiseIfNotString(email, 'email')

        self.__affiliations = Affiliations()
        self.__affiliations.setAncestor(self)

        self.__note = Note()
        self.__note.setAncestor(self)

    @property
    def name(self):
        """
        This method returns that name of the author.

        :returns:   A python str.
        """

        return self.__name

    @property
    def orcid(self):
        """
        This method returns the orcid of the author.

        :returns:   A python str.
        """

        return self.__orcid

    @property
    def email(self):
        """
        This method returns that email of the author.

        :returns:   A python str.
        """

        return self.__email

    @property
    def affiliations(self):
        """
        This method returns a reference to the 'affiliations' suite.

        :returns:       An instance of :py:class:`Affiliations`.
        """

        return self.__affiliations

    @property
    def note(self):
        """
        This method returns a reference to the note member.

        :returns:       An instance of :py:class:`Note`.
        """

        return self.__note

    def XML_extraAttributes(self, **kwargs):
        """
        This method returns the XML attributes for *self* as a single python str.

        :kwargs:        This argument is not used.

        :returns:       A python str.
        """

        attributes = ' name="%s"' % self.__name
        if len(self.__orcid) > 0: attributes += ' orcid="%s"' % self.__orcid
        if len(self.__email) > 0: attributes += ' email="%s"' % self.__email

        return attributes

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:     The minimum amount of indentation.
        :param kwargs:     A dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :return:           List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XMLList = ['%s<%s%s>' % (indent, self.moniker, self.XML_extraAttributes(**kwargs))]
        XMLList += self.__affiliations.toXML_strList(indent=indent2, **kwargs)
        XMLList += self.__note.toXML_strList(indent=indent2, **kwargs)
        XMLList[-1] += '</%s>' % self.moniker

        return XMLList
