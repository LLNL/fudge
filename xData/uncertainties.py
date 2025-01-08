# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the classes for representing an uncertainty instance that is a child node inside a GNDS function.

This module contains the following classes:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | Uncertainty                       | This class represents the uncertainty instance that is a child node   |
    |                                   | inside an xData function.                                             |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Covariance                        |                                                                       |
    +-----------------------------------+-----------------------------------------------------------------------+
    | ListOfCovariances                 |                                                                       |
    +-----------------------------------+-----------------------------------------------------------------------+
"""

from LUPY import ancestry as ancestryModule

from . import base as baseModule
from . import link as linkModule


class Uncertainty(ancestryModule.AncestryIO_bare):
    """
    This class represents the uncertainty instance that is a child node inside an xData function.

    The following table list the primary members of this class:

    +---------------+---------------------------------------------------------------+
    | Member        | Description                                                   |
    +===============+===============================================================+
    | functional    | This is an xData function defining the uncertainty.           |
    +---------------+---------------------------------------------------------------+
    """

    moniker = 'uncertainty'
    keyName = None

    defaultType='variance'
    defaultPdf='normal'
    defaultRelation='offset'

    def __init__(self, functional=None):
        """
        :param functional:  The xData function defining the uncertainty.
        """

        baseModule.XDataCoreMembers.__init__(self)

        self.__functional = functional

    @property
    def data(self):
        """
        This method returns a reference to *self*'s functional.
        """

        return self.__functional

    def convertUnits(self, unitMap):
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        if self.__functional is not None: self.__functional.convertUnits(unitMap)

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        if self.data is None:
            if kwargs.get('showEmpty', False):
                return ['%s<%s/>' % ( indent, self.moniker )]
            return []

        XML_strList = ['%s<%s>' % ( indent, self.moniker )]
        XML_strList += self.__functional.toXML_strList(indent=indent2, **kwargs)
        XML_strList[-1] += '</%s>' % self.moniker

        return XML_strList

    def parseNode(self, node, xPath, linkData, **kwargs):
        """
        This method sets data in *self* using the contents of *node*.

        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        """

        from . import XYs1d as XYs1dModule
        from . import series1d as series1dModule

        xPath.append( node.tag )

        if len(node) != 0:
            if len(node) != 1: raise TypeError("uncertainty node must contain exactly one functional.")

            kwargs = {}
            functionalClass = {
                linkModule.Link.moniker: linkModule.Link,
                XYs1dModule.XYs1d.moniker: XYs1dModule.XYs1d,
                series1dModule.Polynomial1d.moniker: series1dModule.Polynomial1d,
                Covariance.moniker: Covariance,
                ListOfCovariances.moniker: ListOfCovariances,
            }.get(node[0].tag)
            self.__functional = functionalClass.parseNodeUsingClass(node[0], xPath, linkData, **kwargs)
            if self.__functional:
                self.__functional.setAncestor(self)

        xPath.pop()

class Covariance(linkModule.Link):
    """
    This class represents a link to a covariance instance.
    """

    moniker = 'covariance'

class ListOfCovariances(ancestryModule.AncestryIO):

    moniker = 'listOfCovariances'

    def __init__(self):
        ancestryModule.AncestryIO.__init__(self)
        self.__items = []

    def __getitem__(self, item):
        return self.__items[item]

    def add(self, obj):
        if not isinstance(obj, Covariance):
            raise TypeError("Expected Covariance instance, got '%s'" % type(obj))
        obj.setAncestor(self)
        obj.label = 'cov%d' % len(self.__items)
        self.__items.append(obj)

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent','  ')
        XML_strList = ['%s<%s>' % (indent, self.moniker)]
        for item in self: XML_strList += item.toXML_strList(indent2, **kwargs)
        XML_strList[-1] += '</%s>' % self.moniker

        return XML_strList

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

        xPath.append(node.tag)

        listOfCovariances1 = cls()

        for child in node:
            covClass = { Covariance.moniker: Covariance }.get(child.tag)
            listOfCovariances1.add(covClass.parseNodeUsingClass(child, xPath, linkData, **kwargs))

        xPath.pop()

        return listOfCovariances1
