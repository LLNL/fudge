# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from .. import misc as miscModule
from .. import suite as suiteModule
from . import nuclides as nuclidesModule
from . import elapsedTime as elapsedTimeModule

class ProductYield( miscModule.ClassWithLabelKey ) :

    moniker = 'productYield'

    def __init__( self, label, nuclides ) :

        miscModule.ClassWithLabelKey.__init__( self, label )

        self.__nuclides = nuclidesModule.Nuclides( nuclides )
        self.__nuclides.setAncestor( self )

        self.__elapsedTimes = elapsedTimeModule.Suite( )
        self.__elapsedTimes.setAncestor( self )

    @property
    def nuclides(self):

        return self.__nuclides

    @property
    def elapsedTimes(self):

        return self.__elapsedTimes

    @property
    def durations(self):
        # deprecated accessor using GNDS-1.9 name

        return self.__elapsedTimes

    def toXML_strList(self, indent = '', **kwargs):

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s> ' % ( indent, self.moniker ) ]
        XMLStringList += self.__nuclides.toXML_strList( indent2, **kwargs )
        XMLStringList += self.__elapsedTimes.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        for child in element:
            if child.tag == nuclidesModule.Nuclides.moniker:
                self.nuclides.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == elapsedTimeModule.Suite.moniker:
                self.elapsedTimes.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return self

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        self = cls(label='', nuclides=[])
        xPath.pop()
        self.parseNode(element, xPath, linkData, **kwargs)

        return self

class Suite( suiteModule.Suite ) :

    moniker = 'productYields'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, ( ProductYield, ) )
