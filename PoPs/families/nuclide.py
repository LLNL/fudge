# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the nuclear level classes.
"""

from .. import misc as miscModule
from ..quantities import nuclearEnergyLevel as nuclearEnergyLevelModule
from ..groups import misc as chemicalElementMiscModule
from ..fissionFragmentData import fissionFragmentData as fissionFragmentDataModule

from . import particle as particleModule
from . import nucleus as nucleusModule

class alias( particleModule.alias ) :

    moniker = 'nuclideAlias'

    @property
    def chemicalElement( self ) :

        return( self.__particle.chemicalElement )

    @property
    def Z( self ) :

        return( self.__particle.Z )

    @property
    def A( self ) :

        return( self.__particle.A )

    @property
    def index( self ) :

        return( self.__particle.index )

    @property
    def energy( self ) :

        return( self.__particle.energy )

class particle( particleModule.particle ) :

    moniker = 'nuclide'
    familyOrder = 4
    alias = alias

    def __init__( self, id ) :

        particleModule.particle.__init__( self, id )

        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( id )
        self.__nucleus = nucleusModule.particle( chemicalElementMiscModule.nucleusIDFromNuclideID( id ), levelID )
        self.__nucleus.setAncestor( self, attribute = self.keyName )

        self.__fissionFragmentData = fissionFragmentDataModule.fissionFragmentData( )
        self.__fissionFragmentData.setAncestor( self )

    def __lt__( self, other ) :

        if( self.familyOrderLessThan( other ) ) : return( True )
        if( self.familyOrder != other.familyOrder ) : return( False )
        if( self.Z < other.Z ) : return( True )
        if( self.Z != other.Z ) : return( False )
        return( self.A < other.A )

    @property
    def fissionFragmentData( self ) :

        return( self.__fissionFragmentData )

    @property
    def nucleus( self ) :

        return( self.__nucleus )

    @property
    def A( self ) :

        return( self.__nucleus.__A )

    @property
    def chemicalElementSymbol( self ) :

        return( self.__nucleus.chemicalElementSymbol )

    @property
    def index( self ) :

        return( self.__nucleus.index )

    @property
    def isotope( self ) :

        return( self.ancestor.ancestor )

    @property
    def energy( self ) :

        return( self.__nucleus.energy )

    @property
    def Z( self ) :

        return( self.__nucleus.Z )

    def check( self, info ):

        from .. import warning as warningModule
        warnings = []

        subWarnings = self.__nucleus.check(info)
        if subWarnings:
            warnings.append( warningModule.context('nucleus', subWarnings) )
        # FIXME other checks to perform on the nuclide? Will decay info ever live in the nuclide?

        return warnings

    def convertUnits( self, unitMap ) :

        particleModule.particle.convertUnits( self, unitMap )
        self.__nucleus.convertUnits( unitMap )
        self.__fissionFragmentData.convertUnits( unitMap )

    def copy( self ) :

        _particle = particle( self.id )
        self.__copyStandardQuantities( _particle )
        _particle.__nucleus.replicate( self.__nucleus )
        _particle.__fissionFragmentData.replicate( self.__fissionFragmentData )

        return( _particle )

    def extraXMLElements( self, indent, **kwargs ) :

        XMLStringList = self.__nucleus.toXMLList( indent, **kwargs )
        XMLStringList += self.__fissionFragmentData.toXMLList( indent, **kwargs )

        return( XMLStringList )

    def getMass( self, unit ) :

        if( len( self.mass ) > 0 ) : return( self.mass[0].float( unit ) )
        if( self.index == 0 ) : raise Exception( 'Recursion detected as group-state does not have a mass: ID = %s.' % self.id )
        return( self.ancestor[0].mass[0].float( unit ) + self.energy[0].float( unit + ' * c**2' ) )

    def parseExtraXMLElement( self, element, xPath, linkData ) :

        if( element.tag == nucleusModule.particle.moniker ) :
            nucleus = nucleusModule.particle.parseXMLNodeAsClass( element, xPath, linkData )
            self.__nucleus.replicate( nucleus )
            return( True )
        elif( element.tag == fissionFragmentDataModule.fissionFragmentData.moniker ) :
            fissionFragmentData = fissionFragmentDataModule.fissionFragmentData.parseXMLNode( element, xPath, linkData )
            self.__fissionFragmentData.replicate( fissionFragmentData )
            return( True )

        return( False )

    def sortCompare( self, other ) :

        if( not( isinstance( other, particle ) ) ) : raise TypeError( 'Invalid other.' )
        return( self.index - other.index )

class suite( particleModule.suite ) :

    moniker = 'nuclides'
    particle = particle
