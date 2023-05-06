# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""This module defines the ProductYield and ProductYields classes, storing yields for multiple
fission products at various elapsed times after fission."""


from .. import misc as miscModule
from .. import suite as suiteModule
from . import nuclides as nuclidesModule
from . import elapsedTime as elapsedTimeModule

class ProductYield( miscModule.ClassWithLabelKey ) :
    """The ProductYield class stores a list of fission product nuclides and a list of elapsed times
    after fission.  The prompt fission product yield for each nuclide is stored with elapsedTime.time = 0,
    and delayed product yields are listed with elapsedTime.time > 0."""

    moniker = 'productYield'

    def __init__( self, label, nuclides ) :
        """Constructs a ProductYield with specified style label and list of nuclides.

        :param label: string corresponding to a data style, typically the evaluated style.
        :param nuclides: Nuclides instance or list of PoPs particle ids, e.g. ['Mn55', 'Mn56', ...]
        """

        miscModule.ClassWithLabelKey.__init__( self, label )

        self.__nuclides = nuclidesModule.Nuclides( nuclides )
        self.__nuclides.setAncestor( self )

        self.__elapsedTimes = elapsedTimeModule.Suite( )
        self.__elapsedTimes.setAncestor( self )

    @property
    def nuclides(self):
        """Returns the list of nuclides for which fission product yield data are given."""

        return self.__nuclides

    @property
    def elapsedTimes(self):
        """Returns the Suite of elapsedTimes."""

        return self.__elapsedTimes

    @property
    def durations(self):
        """Deprecated accessor using the GNDS-1.9 name. Returns the Suite of elapsedTimes."""

        return self.__elapsedTimes

    def convertUnits(self, unitMap):
        """This method currently does nothing."""
        # FIXME yields are unitless, but we need to support changing units of time.

        pass

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of self.

        :param indent:    The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:    A keyword list.

        :return:          List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XMLStringList = [ '%s<%s> ' % (indent, self.moniker) ]
        XMLStringList += self.__nuclides.toXML_strList(indent2, **kwargs)
        XMLStringList += self.__elapsedTimes.toXML_strList(indent2, **kwargs)
        XMLStringList[-1] += '</%s>' % self.moniker

        return XMLStringList

    def parseNode(self, node, xPath, linkData, **kwargs):
        """Parses contents of the given node into self (a ProductYield instance).

        :param node: pointer to a 'productYield' node in a GNDS file.
        :param xPath: list storing the xPath to the current node, useful for debugging.
        :param linkData: dictionary containing unresolved links and other parsing information."""

        xPath.append(node.tag)
        for child in node:
            if child.tag == nuclidesModule.Nuclides.moniker:
                self.nuclides.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == elapsedTimeModule.Suite.moniker:
                self.elapsedTimes.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, node.tag))

        xPath.pop()
        return self

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Create a new instance of class **cls** and parse contents of node into the instance.

        :param cls: FUDGE Python class to return.
        :param node: **FissionFragmentData** node to parse.
        :param xPath: List containing xPath to current node, useful mostly for debugging.
        :param linkData: dict that collects unresolved links.
        """

        xPath.append(node.tag)
        self = cls(label='', nuclides=[])
        xPath.pop()
        self.parseNode(node, xPath, linkData, **kwargs)

        return self

class Suite(suiteModule.Suite):
    """Contains one or more ProductYield instances, each with a label corresponding to a unique **style**."""

    moniker = 'productYields'

    def __init__(self):

        suiteModule.Suite.__init__(self, (ProductYield,))
