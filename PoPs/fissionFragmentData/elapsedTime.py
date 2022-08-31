# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the duration and durations classes,
used to organize independent and cumulative fission product yields by the time since fission.

The duration consists of a time constant (i.e. elapsed time since fission) and a set of product yields.
It may also contain one or more Q-matrices, which can be used to transform fission product yields to
a later time constant.
"""

from .. import misc as miscModule
from .. import suite as suiteModule
from . import time as timeModule
from . import yields as yieldsModule
from . import QMatrix as QMatrixModule

class ElapsedTime(miscModule.ClassWithLabelKey):

    moniker = 'elapsedTime'

    def __init__(self, label):

        miscModule.ClassWithLabelKey.__init__(self, 'label')

        self.__label = label

        self.__time = timeModule.Suite( )
        self.__time.setAncestor( self )

        self.__yields = yieldsModule.Yields( )
        self.__yields.setAncestor( self )

        self.__QMatrix = QMatrixModule.Suite( )
        self.__QMatrix.setAncestor( self )

    @property
    def label(self) :

        return self.__label

    @property
    def time(self) :

        return self.__time

    @property
    def yields(self):

        return self.__yields

    @property
    def QMatrix(self):

        return self.__QMatrix

    def sortCompare(self, other):

        return 1

    def toXML_strList(self, indent = '', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XMLStringList = ['%s<%s label="%s">' % (indent, self.moniker, self.__label)]
        XMLStringList += self.time.toXML_strList(indent = indent2, **kwargs)
        XMLStringList += self.yields.toXML_strList(indent = indent2, **kwargs)
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        for child in element:
            if child.tag == timeModule.Suite.moniker:
                self.time.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == yieldsModule.Yields.moniker:
                self.yields.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return self

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        self = cls(element.get('label'))
        xPath.pop()
        self.parseNode(element, xPath, linkData, **kwargs)

        return self


class Suite(suiteModule.SortedSuite):

    moniker = 'elapsedTimes'

    def __init__(self):

        suiteModule.SortedSuite.__init__(self, allowedClasses = (ElapsedTime,))
