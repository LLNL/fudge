# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from fudge import suites as suitesModule

from PoPs.fissionFragmentData import nuclides as nuclidesModule
from . import elapsedTime as elapsedTimeModule

class ProductYield( ancestryModule.AncestryIO ) :
    """
    Stores the fission product yields at one or more elapsed times
    (e.g. 'prompt' and 'delayed') and various incident energies.

    If the list of nuclides is common to all elapsedTime / incidentEnerg
    they may be stored at this level, otherwise they are found inside
    each elapsedTime / incidentEnergy.
    """
    # FIXME: lots of overlap with the version in PoPs/fissionFragmentData

    moniker = 'productYield'

    def __init__( self, label, nuclides=None ) :

        ancestryModule.AncestryIO.__init__( self )

        self.__label = label

        if nuclides is None:
            nuclides = []
        self.__nuclides = nuclidesModule.Nuclides(nuclides)
        self.__nuclides.setAncestor(self)

        self.__elapsedTimes = elapsedTimeModule.Suite( )
        self.__elapsedTimes.setAncestor( self )

    @property
    def label(self):

        return self.__label

    @property
    def nuclides(self):

        return self.__nuclides

    @property
    def elapsedTimes(self):

        return self.__elapsedTimes

    @property
    def durations(self):
        # deprecated property using old GNDS-1.9 name

        return self.__elapsedTimes

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s label="%s"> ' % (indent, self.moniker, self.__label) ]
        XMLStringList += self.__nuclides.toXML_strList(indent2, **kwargs)
        XMLStringList += self.__elapsedTimes.toXML_strList(indent2, **kwargs)
        XMLStringList[-1] += '</%s>' % self.moniker

        return XMLStringList

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        PY = cls(element.get('label'))
        for child in element:
            if child.tag == nuclidesModule.Nuclides.moniker:
                PY.nuclides.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == elapsedTimeModule.Suite.moniker:
                PY.elapsedTimes.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return PY

class Suite( suitesModule.Suite ) :

    moniker = 'productYields'

    def __init__( self ) :

        suitesModule.Suite.__init__( self, ( ProductYield, ) )
