# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the FissionFragmentData class, used to store fission product yields
and the delayed neutrons and gammas emitted by the decay of those fission products.
"""

from LUPY import ancestry as ancestryModule

from . import delayedNeutronData as delayedNeutronDataModule
from . import productYield as productYieldModule

class FissionFragmentData(ancestryModule.AncestryIO):

    moniker = 'fissionFragmentData'

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__(self)

        self.__delayedNeutronData = delayedNeutronDataModule.DelayedNeutronData( )
        self.__delayedNeutronData.setAncestor( self )

        self.__productYields = productYieldModule.Suite( )
        self.__productYields.setAncestor( self )
        
    @property
    def delayedNeutronData( self ) :

        return( self.__delayedNeutronData )

    @property
    def productYields( self ) :

        return( self.__productYields )

    def convertUnits( self, unitMap ) :

        self.__delayedNeutronData.convertUnits( unitMap )
        self.__productYields.convertUnits( unitMap )

    def replicate( self, other ) :

        self.__delayedNeutronData = other.delayedNeutronData    # FIXME not making a copy?
        self.__delayedNeutronData.setAncestor( self )

        self.__productYields = other.productYields
        self.__productYields.setAncestor( self )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s> ' % ( indent, self.moniker ) ]
        XMLStringList += self.__delayedNeutronData.toXML_strList( indent2, **kwargs )
        XMLStringList += self.__productYields.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        if( len( XMLStringList ) == 1 ) : XMLStringList = []

        return( XMLStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )
        FFD = cls()
        for child in element:
            if child.tag == delayedNeutronDataModule.DelayedNeutronData.moniker:
                FFD.delayedNeutronData.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == productYieldModule.Suite.moniker:
                FFD.productYields.parseNode(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Encountered unknown node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return FFD
