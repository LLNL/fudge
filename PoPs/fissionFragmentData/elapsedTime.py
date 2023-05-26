# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the ElapsedTime and ElapsedTimes classes,
used to organize independent and cumulative fission product yields by the time elapsed since fission.

Note: the GNDS-1.9 name 'duration' was renamed to 'elapsedTime' in GNDS-2.0.
GNDS-2.0 also removed the 'QMatrix', an unused method for transforming product yields to a later elapsedTime.
"""

from .. import misc as miscModule
from .. import suite as suiteModule
from . import time as timeModule
from . import yields as yieldsModule

#from . import QMatrix as QMatrixModule

class ElapsedTime(miscModule.ClassWithLabelKey):
    """Contains fission product yields at a specified time, measured since fission.

    Class contents:
     - time (elapsed time value and unit. Use time=0 for prompt yields)
     - yields (fission product yields at the given time)"""

    moniker = 'elapsedTime'

    def __init__(self, label):
        """Construct an ElapsedTime with the given label.
        The label is arbitrary but must be unique within the parent ElapsedTimes suite."""

        miscModule.ClassWithLabelKey.__init__(self, 'label')

        self.__label = label

        # CMM should this really be a Suite? Only one value should be allowed at this level,
        # since we're already inside a productYield corresponding to a specific style.
        self.__time = timeModule.Suite()
        self.__time.setAncestor(self)

        self.__yields = yieldsModule.Yields()
        self.__yields.setAncestor(self)

        #self.__QMatrix = QMatrixModule.Suite( )
        #self.__QMatrix.setAncestor( self )

    @property
    def label(self):
        """Returns string label corresponding to this elapsedTime."""

        return self.__label

    @property
    def time(self) :

        return self.__time

    @property
    def yields(self):

        return self.__yields

    """
    @property
    def QMatrix(self):

        return self.__QMatrix
    """

    # FIXME CMM is this used? Seems to be a stub.
    def sortCompare(self, other):

        return 1

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of self.

        :param indent:    The amount of indentation for each line. Child nodes and text may be indented more.
        :param kwargs:    A keyword list.

        :return:          List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XMLStringList = ['%s<%s label="%s">' % (indent, self.moniker, self.__label)]
        XMLStringList += self.time.toXML_strList(indent = indent2, **kwargs)
        XMLStringList += self.yields.toXML_strList(indent = indent2, **kwargs)
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseNode(self, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)
        for child in node:
            if child.tag == timeModule.Suite.moniker:
                self.time.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == yieldsModule.Yields.moniker:
                self.yields.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, node.tag))

        xPath.pop()
        return self

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """Create a new instance of class **cls** and parse contents of node into the instance.

        :param cls: FUDGE Python class to return.
        :param node: **elapsedTime** node to parse.
        :param xPath: List containing xPath to current node, useful mostly for debugging.
        :param linkData: dict that collects unresolved links.
        """

        xPath.append(node.tag)
        self = cls(node.get('label'))
        xPath.pop()
        self.parseNode(node, xPath, linkData, **kwargs)

        return self


class Suite(suiteModule.SortedSuite):
    """Contains a list of ElapsedTime instances."""

    moniker = 'elapsedTimes'

    def __init__(self):

        suiteModule.SortedSuite.__init__(self, allowedClasses = (ElapsedTime,))
