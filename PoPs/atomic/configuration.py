# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the atomic configuration and configurations classes.
These describe electron subshells inside atom, including the binding energy and decay info for each subshell.

Each subshell includes an electronNumber, i.e. the number of electrons occupying that shell for a neutral atom.
The electronNumber is not necessarily an integer. In Fluorine, for example, the final valence electron may be in either
the 2p1/2 or 2p3/2 subshells with probabilities 2/3 and 1/3 respectively.
"""

import abc

from .. import suite as suiteModule
from .. import misc as miscModule
from ..decays import decayData as decayDataModule
from ..quantities import bindingEnergy as bindingEnergyModule

class Configuration( miscModule.ClassWithSubshellKey ) :

    moniker = 'configuration'

    def __init__( self, subshell, electronNumber ) :
        """
        :param subshell: label for an electron subshell (e.g. 1s1/2, 2p1/2).
        :param electronNumber: number of electrons in the subshell when neutral
        """

        miscModule.ClassWithSubshellKey.__init__( self, subshell )

        if( not( isinstance( electronNumber, float ) ) ) : raise TypeError( 'electronNumber is not a float' )
        self.__electronNumber = electronNumber

        self.__bindingEnergy = bindingEnergyModule.Suite( )
        self.__bindingEnergy.setAncestor( self )

        self.__decayData = decayDataModule.DecayData( )
        self.__decayData.setAncestor( self )

    @property
    def electronNumber( self ) :

        return( self.__electronNumber )

    @property
    def bindingEnergy( self ) :
        """The total binding energy for the subshell"""

        return( self.__bindingEnergy )

    @property
    def decayData( self ) :

        return( self.__decayData )

    def convertUnits( self, unitMap ) :
        """See documentation in PoPs.database"""

        self.__bindingEnergy.convertUnits( unitMap )
        self.__decayData.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of self
        """

        _configuration = self.__class__( self.subshell, self.electronNumber )
        self.__decayData.copyItems( _configuration.decayData )
        return( _configuration )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s subshell="%s" electronNumber="%s">' % ( indent, self.moniker, self.subshell, self.electronNumber ) ]
        XMLStringList += self.__bindingEnergy.toXML_strList( indent = indent2, **kwargs )
        XMLStringList += self.__decayData.toXML_strList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        for child in element :
            if( child.tag == bindingEnergyModule.Suite.moniker ) :
                self.bindingEnergy.parseNode(child, xPath, linkData, **kwargs)
            elif( child.tag == decayDataModule.DecayData.moniker ) :
                self.decayData.parseNode(child, xPath, linkData, **kwargs)
            else :
                raise ValueError( 'Invalid tag = "%s"' % child.tag )

        xPath.pop( )
        return( self )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        self = cls( element.get('subshell'), float( element.get('electronNumber') ) )
        self.parseNode(element, xPath, linkData, **kwargs)

        xPath.pop( )
        return( self )

class Configurations( suiteModule.Suite ) :

    moniker = 'configurations'

    def __init__( self, replace = True ) :

        suiteModule.Suite.__init__( self, ( Configuration, ), replace = replace )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        for child in element :
            self.add( Configuration.parseNodeUsingClass(child, xPath, linkData, **kwargs))

        xPath.pop( )
        return( self )
