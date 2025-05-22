# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the GNDS documentation child node computerCodes and its child nodes classes.

This module contains the following classes:

    +---------------------------+-----------------------------------------------------------------------------------+
    | Class                     | Description                                                                       |
    +===========================+===================================================================================+
    | Version                   | This is the class for the GNDS computerCode version attribute.                    |
    +---------------------------+-----------------------------------------------------------------------------------+
    | CodeRepo                  | This is the class for the GNDS computerCode/codeRepo node.                        |
    +---------------------------+-----------------------------------------------------------------------------------+
    | ExecutionArguments        | This is the class for the GNDS computerCode/executionArguments node.              |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Note                      | This is the class for the GNDS computerCode/note node.                            |
    +---------------------------+-----------------------------------------------------------------------------------+
    | InputDeck                 | This is the class for the GNDS computerCode/inputDecks/inputDeck node.            |
    +---------------------------+-----------------------------------------------------------------------------------+
    | InputDecks                | This is the suite class for the GNDS computerCode/inputDecks node.                |
    +---------------------------+-----------------------------------------------------------------------------------+
    | OutputDeck                | This is the class for the GNDS computerCode/outputDecks/outputDeck node.          |
    +---------------------------+-----------------------------------------------------------------------------------+
    | OutputDecks               | This is the suite class for the GNDS computerCode/outputDecks node.               |
    +---------------------------+-----------------------------------------------------------------------------------+
    | ComputerCode              | This is the class for the GNDS documentation/computerCodes/computerCode node.     |
    +---------------------------+-----------------------------------------------------------------------------------+
    | ComputerCodes             | This is the suite class for the GNDS documentation/computerCodes node.            |
    +---------------------------+-----------------------------------------------------------------------------------+
"""

import datetime

from LUPY import ancestry as ancestryModule

from .. import suite as suiteModule
from .. import text as textModule
from .. import date as dateModule


class Version(textModule.Text):
    """
    This is the class for the GNDS computerCode version attribute.

    The following table list the primary members of this class added by this derived class:

    +-----------------------+---------------------------------------------------------------------------+
    | Member                | Description                                                               |
    +=======================+===========================================================================+
    | date                  | This is the date the code was checkout of the repo.                       |
    +-----------------------+---------------------------------------------------------------------------+
    """

    moniker = 'version'

    def __init__(self, date=None):
        """
        :param date:        An instance of :py:class:`dateModule.Date`.
        """

        textModule.Text.__init__(self, markup=textModule.Markup.none)

        self.__date = dateModule.Date(date)

    @property
    def date(self):
        """
        This method returns the *date* member of *self*.

        :returns:       An instance of :py:class:`dateModule.Date`.
        """

        return self.__date

    @date.setter
    def date(self, date):
        """
        This method sets the *date* member of *self* to *date*.

        :param date:   An instance of :py:class:`dateModule.Date`.
        """

        self.__date = dateModule.raiseIfNotDate(date)

    def XML_extraAttributes(self, **kwargs):
        """
        This method returns the XML attributes for *self* as a python str.

        :kwargs:        A python dict with key 'addExtraAttributes', used to determine if the *date* attribute is added.

        :returns:       A python str.
        """

        if not kwargs['addExtraAttributes']: return ''
        return ' date="%s"' % self.date

    def parseNode(self, node, xPath, linkData, **kwargs):
        """
        This method fills *self* by parsing the data in *node*.

        :param node:       Node to parse.
        :param xPath:      List containing xPath to current node, useful mostly for debugging.
        :param linkData:   dict that collects unresolved links.
        :param kwargs:     A dictionary of extra arguments controlling how *self* is converted to a list of XML strings.
        """

        textModule.Text.parseNode(self, node, xPath, linkData, **kwargs)
        xPath.append(node.tag)

        self.date = datetime.datetime.strptime(node.get('date'), '%Y-%m-%d')

        xPath.pop()


class CodeRepo(ancestryModule.AncestryIO_bare):
    """
    This is the class for the GNDS computerCode/codeRepo node.

    The following table list the primary members of this class:

    +-----------------------+---------------------------------------------------------------------------+
    | Member                | Description                                                               |
    +=======================+===========================================================================+
    | revisionSystem        | The name of the repo system used (e.g., git, subversion).                 |
    +-----------------------+---------------------------------------------------------------------------+
    | href                  | A url to the repo.                                                        |
    +-----------------------+---------------------------------------------------------------------------+
    | revisionID            | A unique identifier of the checkout version of the code.                  |
    +-----------------------+---------------------------------------------------------------------------+
    """

    moniker = 'codeRepo'

    def __init__(self, _revisionSystem, _href, _revisionID):
        """
        :param _revisionSystem:     The name of the repo system used (e.g., git, subversion). 
        :param _href:               URL to the repo.
        :param _revisionID:         A unique identifier of the checkout version of the code.
        """

        ancestryModule.AncestryIO_bare.__init__(self)

        self.__revisionSystem = textModule.raiseIfNotString(_revisionSystem, 'revisionSystem')
        self.__href = textModule.raiseIfNotString(_href, 'href')
        self.__revisionID = textModule.raiseIfNotString(_revisionID, 'revisionID')

    @property
    def revisionSystem(self):
        """
        This method returns the *revisionSystem* member of *self*.

        :returns:       A python str.
        """

        return self.__revisionSystem

    @property
    def href(self):
        """
        This method returns the *href* member of *self*.

        :returns:       A python str.
        """

        return self.__href

    @property
    def revisionID(self):
        """
        This method returns the *revisionID* member of *self*.

        :returns:       A python str.
        """

        return self.__revisionID

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:     The minimum amount of indentation.
        :param kwargs:     A dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :return:           List of str instances representing the XML lines of self.
        """

        if len(self.__revisionSystem + self.__href + self.__revisionID) == 0:
            if kwargs.get('showEmpty', False):
                return ['%s<!-- code repo -->' % indent]
            return []

        return ['%s<%s revisionSystem="%s" href="%s" revisionID="%s"/>' % (
            indent, self.moniker, self.__revisionSystem, self.__href, self.__revisionID)]

    def parseNode(self, node, xPath, linkData, **kwargs):
        """
        This method fills *self* by parsing the data in *node*.

        :param node:       Node to parse.
        :param xPath:      List containing xPath to current node, useful mostly for debugging.
        :param linkData:   dict that collects unresolved links.
        :param kwargs:     A dictionary of extra arguments controlling how *self* is converted to a list of XML strings.
        """

        xPath.append(node.tag)

        self.__revisionSystem = textModule.raiseIfNotString(node.get('revisionSystem'), 'revisionSystem')
        self.__href = textModule.raiseIfNotString(node.get('href'), 'href')
        self.__revisionID = textModule.raiseIfNotString(node.get('revisionID'), 'revisionID')

        xPath.pop()


class ExecutionArguments(textModule.Text):
    """
    This is the class for the GNDS computerCode/executionArguments node.
    """

    moniker = 'executionArguments'

    def __init__(self):
        textModule.Text.__init__(self, markup=textModule.Markup.none)


class Note(textModule.Text):
    """
    This is the class for the GNDS computerCode/note node.
    """

    moniker = 'note'


class InputDeck(textModule.Text):
    """
    This is the class for the GNDS computerCode/inputDecks/inputDeck node.

    The following table list the primary members of this class:

    +-----------------------+---------------------------------------------------------------------------+
    | Member                | Description                                                               |
    +=======================+===========================================================================+
    | label                 | A unique label for an instance within a inputDecks suite.                 |
    +-----------------------+---------------------------------------------------------------------------+
    | filename              | The path of the output file.                                              |
    +-----------------------+---------------------------------------------------------------------------+
    | text                  | Optional text to describe the contents of the output file.                |
    +-----------------------+---------------------------------------------------------------------------+
    """

    moniker = 'inputDeck'
    keyName = 'label'

    def __init__(self, label, filename, text=None):
        """
        :param label:       A unique label for an instance within a inputDecks suite.
        :param filename:    The path of the output file.
        :param text:        Optional text to describe the contents of the output file.
        """

        textModule.Text.__init__(self, text, markup=textModule.Markup.none, label=label)

        self.__filename = textModule.raiseIfNotString(filename, 'filename')

    @property
    def filename(self):
        """
        This method returns the *filename* member of *self*.

        :returns:       An instance of :py:class:`textModule.Tex`.
        """

        return self.__filename

    def XML_extraAttributes(self, **kwargs):
        """
        This method returns the XML attributes for *self* as a single python str.

        :kwargs:        This argument is not used.

        :returns:       A python str.
        """

        if self.filename == '': return ''

        return ' filename="%s"' % self.filename

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:        Form class to return.
        :param node:       Node to parse.
        :param xPath:      List containing xPath to current node, useful mostly for debugging.
        :param linkData:   dict that collects unresolved links.
        :param kwargs:     A dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :returns:          An instance of *cls* representing *node*.
        """

        _label = node.get('label')
        _filename = node.get('filename', '')

        inputDeck = cls(_label, _filename)
        inputDeck.parseNode(node, xPath, linkData, **kwargs)

        return inputDeck


class InputDecks(suiteModule.Suite):
    """
    This is the suite class for the GNDS computerCode/inputDecks node.
    """

    moniker = 'inputDecks'
    suiteName = 'label'

    def __init__(self):
        suiteModule.Suite.__init__(self, [InputDeck])


class OutputDeck(textModule.Text):
    """
    This is the class for the GNDS computerCode/outputDecks/outputDeck node.

    The following table list the primary members of this class:

    +-----------------------+---------------------------------------------------------------------------+
    | Member                | Description                                                               |
    +=======================+===========================================================================+
    | filename              | The path of the output file.                                              |
    +-----------------------+---------------------------------------------------------------------------+
    | text                  | Optional text to describe the contents of the output file.                |
    +-----------------------+---------------------------------------------------------------------------+
    """

    moniker = 'outputDeck'
    keyName = 'label'

    def __init__(self, filename, text=None):
        """
        :param filename:    The path of the output file.
        :param text:        Optional text to describe the contents of the output file.
        """

        textModule.Text.__init__(self, text, markup=textModule.Markup.none)

        self.__filename = textModule.raiseIfNotString(filename, 'filename')

    @property
    def filename(self):
        """
        This method returns the *filename* member of *self*.

        :returns:       An instance of :py:class:`textModule.Tex`.
        """

        return self.__filename

    def XML_extraAttributes(self, **kwargs):
        """
        This method returns the XML attributes for *self* as a single python str.

        :kwargs:        This argument is not used.

        :returns:       A python str.
        """

        if self.filename == '': return ''

        return ' filename="%s"' % self.filename

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:        Form class to return.
        :param node:       Node to parse.
        :param xPath:      List containing xPath to current node, useful mostly for debugging.
        :param linkData:   dict that collects unresolved links.
        :param kwargs:     A dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :returns:          An instance of *cls* representing *node*.
        """

        _filename = node.get('filename', '')

        return cls(_filename)


class OutputDecks(suiteModule.Suite):
    """
    This is the suite class for the GNDS computerCode/outputDecks node.
    """

    moniker = 'outputDecks'
    suiteName = 'label'

    def __init__(self):
        suiteModule.Suite.__init__(self, [OutputDeck])


class ComputerCode(ancestryModule.AncestryIO):
    """
    This is the class for the GNDS documentation/computerCodes/computerCode node.

    The following table list the primary members of this class:

    +-----------------------+---------------------------------------------------------------------------+
    | Member                | Description                                                               |
    +=======================+===========================================================================+
    | label                 | A unique label for an instance within a computerCodes suite.              |
    +-----------------------+---------------------------------------------------------------------------+
    | name                  | The name of the code.                                                     |
    +-----------------------+---------------------------------------------------------------------------+
    | version               | The version of the code.                                                  |
    +-----------------------+---------------------------------------------------------------------------+
    | codeRepo              | A url to the repo for the code.                                           |
    +-----------------------+---------------------------------------------------------------------------+
    | executionArguments    | A suite containing a list of :py:class:`ExecutionArgument` instances.     |
    +-----------------------+---------------------------------------------------------------------------+
    | note                  | A text node with note about the code.                                     |
    +-----------------------+---------------------------------------------------------------------------+
    | inputDecks            | A suite containing a list of :py:class:`InputDeck` instances.             |
    +-----------------------+---------------------------------------------------------------------------+
    | outputDecks           | A suite containing a list of :py:class:`OutputDeck` instances.            |
    +-----------------------+---------------------------------------------------------------------------+
    """

    # Also need buildParameters as a text node.

    moniker = 'computerCode'
    keyName = 'label'
    ancestryMembers = ('codeRepo', 'executionArguments', 'note', 'inputDecks', 'outputDecks')

    def __init__(self, label, name, version):
        """
        :param label:       A unique label for an instance within a computerCodes suite.
        :param name:        The name of the code.
        :param version:     The version of the code.
        """

        ancestryModule.AncestryIO.__init__(self)

        self.__label = textModule.raiseIfNotString(label, 'label')
        self.__name = textModule.raiseIfNotString(name, 'name')
        self.__version = textModule.raiseIfNotString(version, 'version')

        self.codeRepo = CodeRepo('', '', '')

        self.__executionArguments = ExecutionArguments()
        self.__executionArguments.setAncestor(self)

        self.__note = Note()
        self.__note.setAncestor(self)

        self.__inputDecks = InputDecks()
        self.__inputDecks.setAncestor(self)

        self.__outputDecks = OutputDecks()
        self.__outputDecks.setAncestor(self)

    @property
    def label(self):
        """
        This method returns a reference to the *label* member of *self*.

        :returns:       A python str.
        """

        return self.__label

    @property
    def name(self):
        """
        This method returns a reference to the *name* member of *self*.

        :returns:       A python str.
        """

        return self.__name

    @property
    def version(self):
        """
        This method returns a reference to the *version* member of *self*.

        :returns:       A python str.
        """

        return self.__version

    @property
    def codeRepo(self):
        """
        This method returns a reference to the *codeRepo* member of *self*.

        :returns:       An instance of :py:class:`CodeRepo`.
        """

        return self.__codeRepo

    @codeRepo.setter
    def codeRepo(self, _codeRepo):
        """
        This method sets the *codeRepo* member of *self* to *_codeRepo*.

        :param _codeRepo:   An instance of :py:class:`CodeRepo`.
        """

        if not (isinstance(_codeRepo, CodeRepo)): raise TypeError('Invalid codeRepo instance.')
        self.__codeRepo = _codeRepo
        self.__codeRepo.setAncestor(self)

    @property
    def executionArguments(self):
        """
        This method returns a reference to the *executionArguments* member of *self*.

        :returns:       An instance of :py:class:`ExecutionArguments`.
        """

        return self.__executionArguments

    @property
    def note(self):
        """
        This method returns a reference to the *note* member of *self*.

        :returns:       An instance of :py:class:`Note`.
        """

        return self.__note

    @property
    def inputDecks(self):
        """
        This method returns a reference to the *inputDecks* member of *self*.

        :returns:       An instance of :py:class:`InputDecks`.
        """

        return self.__inputDecks

    @property
    def outputDecks(self):
        """
        This method returns a reference to the *outputDecks* member of *self*.

        :returns:       An instance of :py:class:`OutputDecks`.
        """

        return self.__outputDecks

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:     The minimum amount of indentation.
        :param kwargs:     A dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :return:           List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XMLList = ['%s<%s label="%s" name="%s" version="%s">' % (
            indent, self.moniker, self.__label, self.__name, self.__version)]
        XMLList += self.__codeRepo.toXML_strList(indent2, **kwargs)

        XMLList += self.__executionArguments.toXML_strList(indent2, **kwargs)
        XMLList += self.__note.toXML_strList(indent2, **kwargs)
        XMLList += self.__inputDecks.toXML_strList(indent2, **kwargs)
        XMLList += self.__outputDecks.toXML_strList(indent2, **kwargs)

        XMLList[-1] += '</%s>' % self.moniker

        return XMLList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:        Form class to return.
        :param node:       Node to parse.
        :param xPath:      List containing xPath to current node, useful mostly for debugging.
        :param linkData:   dict that collects unresolved links.
        :param kwargs:     A dictionary of extra arguments controlling how *self* is converted to a list of XML strings.

        :returns:          An instance of *cls* representing *node*.
        """

        label = node.get('label')
        name = node.get('name')
        version = node.get('version', '')

        computerCode = cls(label, name, version)
        computerCode.parseAncestryMembers(node, xPath, linkData, **kwargs)

        return computerCode


class ComputerCodes(suiteModule.Suite):
    """
    This is the suite class for the GNDS documentation/computerCodes node.
    """

    moniker = 'computerCodes'
    suiteName = 'label'

    def __init__(self):
        suiteModule.Suite.__init__(self, [ComputerCode])
