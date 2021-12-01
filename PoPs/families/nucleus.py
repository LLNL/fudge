# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the nucleus class definitions.
The nucleus consists only of baryons, so quantities like spin, charge, etc. do not include atomic electrons.
For the compound particle formed of a nucleus + electrons, see the `nuclide` module.
"""

from .. import misc as miscModule

from ..groups import misc as chemicalElementMiscModule

from ..quantities import nuclearEnergyLevel as nuclearEnergyLevelModule

from . import particle as particleModule

class alias( particleModule.alias ) :

    moniker = 'nucleusAlias'

    @property
    def chemicalElementSymbol( self ) :

        return( self.__particle.chemicalElementSymbol )

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

    moniker = 'nucleus'
    familyOrder = 3
    alias = alias

    def __init__( self, id, index ) :
        """
        :param id: Unique id for this nucleus (string). Per naming convention, should be similar to 'he4'
        :param index: Nuclear excited level index (int). 0 = ground state, 1 = 1st excited, etc.
        """

        particleModule.particle.__init__( self, id )

        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( id )
        if( not( isNucleus ) ) : raise ValueError( 'Invalid nucleus id = "%s"' % id )

        index = int( index )
        if( index != levelID ) : raise ValueError( 'id = "%s" does not agree with index = "%s"' % 
                ( miscModule.toLimitedString( id ), miscModule.toLimitedString( index ) ) )

        self.__chemicalElementSymbol = chemicalElementSymbol
        self.__Z = chemicalElementMiscModule.ZFromSymbol[chemicalElementSymbol]
        self.__A = chemicalElementMiscModule.checkA( A )
        self.__index = chemicalElementMiscModule.checkIndex( index )
        self.__energy = self.addSuite( nuclearEnergyLevelModule )

    def __lt__( self, other ) :

        if( self.familyOrderLessThan( other ) ) : return( True )        # Also checks isinstance of other.
        if( self.familyOrder != other.familyOrder ) : return( False )
        if( self.Z < other.Z ) : return( True )
        if( self.Z != other.Z ) : return( False )
        return( self.A < other.A )

    @property
    def A( self ) :

        return( self.__A )

    @property
    def chemicalElementSymbol( self ) :

        return( self.__chemicalElementSymbol )

    @property
    def index( self ) :

        return( self.__index )

    @property
    def energy( self ) :

        return( self.__energy )

    @property
    def nuclide( self ) :

        return( self.ancestor )

    @property
    def Z( self ) :

        return( self.__Z )

    def convertUnits( self, unitMap ) :
        """ See convertUnits documentation in PoPs.database """

        particleModule.particle.convertUnits( self, unitMap )
        self.__energy.convertUnits( unitMap )

    def copy( self ) :
        """
        :return: deep copy of self
        """

        _particle = particle( self.id, self.index )
        self.__copyStandardQuantities( _particle )
        for item in self.__energy : _particle.energy.add( item.copy( ) )
        return( _particle )

    def getMass( self, unit ) :
        """
        Evaluate the nucleus mass in the desired unit. Mass is adjusted for excitation energy,
        but not for electron masses / binding energy.
        :param unit: desired unit (string)
        :return: mass (float)
        """

# FIXME Still need to correct for electron masses and binding energy.
        if( len( self.mass ) > 0 ) : return( self.mass[0].float( unit ) )
        return( self.nuclide.getMass( unit ) )

    def replicate( self, other ) :
        """
        Copy data from other into self
        :param other: another nucleus.particle instance
        """

        particleModule.particle.replicate( self, other )
        self.__energy.replicate( other.energy )

    def extraXMLAttributes( self ) :

        return( ' index="%s"' % self.index )

    def extraXMLElements( self, indent, **kwargs ) :

        return( self.energy.toXMLList( indent, **kwargs ) )

    def parseExtraXMLElement( self, element, xPath, linkData ) :

        if( element.tag == nuclearEnergyLevelModule.suite.moniker ) :
            nuclearEnergyLevelModule.suite.parseXMLNode( self.energy, element, xPath, linkData )
            return( True )

        return( False )

    def sortCompare( self, other ) :

        if( not( isinstance( other, particle ) ) ) : raise TypeError( 'Invalid other.' )
        return( self.index - other.index )
