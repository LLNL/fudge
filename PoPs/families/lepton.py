# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the lepton class definitions.
Leptons include electrons, neutrinos, etc.
"""

from LUPY import enums as enumsModule

from .. import misc as miscModule
from . import particle as particleModule

class Generation(enumsModule.Enum):

    electronic = enumsModule.auto()
    muonic = enumsModule.auto()
    tauonic = enumsModule.auto()

class Alias( particleModule.Alias ) :

    moniker = 'leptonAlias'

class Particle( particleModule.Particle ) :

    moniker = 'lepton'
    familyOrder = 1
    alias = Alias
    generations = { Generation.electronic : 0, Generation.muonic : 1,  Generation.tauonic: 2 }  # Only used for sorting particles.
    chargeOrder = { -1 : 0, 1 : 1, 0 : 2 }                                                      # Only used for sorting particles.

    def __init__( self, id, **kwargs ) :

        particleModule.Particle.__init__( self, id )

        if( 'generation' not in kwargs ) : raise ValueError( 'generation attribute not present' )
        generation = kwargs['generation']
        if isinstance(generation, str):
            generation = Generation.fromString(generation)
        if generation not in Generation:
            raise Exception( 'Invalid generation "%s" for particle "%s" of family "%s"' % ( generation, id, self.family ) )
        self.__generation = generation

    def __lt__( self, other ) :

        if( self.familyOrderLessThan( other ) ) : return( True )        # Also checks isinstance of other.
        if( self.familyOrder != other.familyOrder ) : return( False )
        if( self.generation < other.generation ) : return( True )
        if( self.chargeOrder[self.charge[0].value] < self.chargeOrder[other.charge[0].value] ) : return( True )
        return( False )

    @property
    def generation( self ) :

        return( self.__generation )

    def copy( self ) :

        _particle = Particle( self.id, generation = self.generation )
        self.__copyStandardQuantities( _particle )
        return( _particle )

    def extraXMLAttributes( self ) :

        return( ' generation="%s"' % self.generation )

    def intid(self, intidDB={}):

        sign = -1 if self.isAnti else 1

        base, anti, qualifier = miscModule.baseAntiQualifierFromID(self.id)
        leptonIndex = {'e-': 0}.get(base)
        if leptonIndex is None:
            ValueError('Baryon "%s" does not have a defined intid.' % (self.id))

        return sign * (2**30 + 2**26 + leptonIndex)

class Suite( particleModule.Suite ) :

    moniker = 'leptons'
    particle = Particle
