# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the fissionFragmentData class, used to store fission product yields
and the delayed neutrons and gammas emitted by the decay of those fission products.
"""
from xData import ancestry as ancestryModule

from . import delayedNeutronData as delayedNeutronDataModule
from . import productYield as productYieldModule

class fissionFragmentData( ancestryModule.ancestry ) :

    moniker = 'fissionFragmentData'

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

        self.__delayedNeutronData = delayedNeutronDataModule.delayedNeutronData( )
        self.__delayedNeutronData.setAncestor( self )

        self.__productYields = productYieldModule.suite( )
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

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s> ' % ( indent, self.moniker ) ]
        XMLStringList += self.__delayedNeutronData.toXMLList( indent2, **kwargs )
        XMLStringList += self.__productYields.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        if( len( XMLStringList ) == 1 ) : XMLStringList = []

        return( XMLStringList )

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):

        xPath.append( element.tag )
        FFD = cls()
        for child in element:
            if child.tag == delayedNeutronDataModule.delayedNeutronData.moniker:
                FFD.delayedNeutronData.parseXMLNode(child, xPath, linkData)
            elif child.tag == productYieldModule.suite.moniker:
                FFD.productYields.parseXMLNode(child, xPath, linkData)
            else:
                raise TypeError("Encountered unknown node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return FFD
