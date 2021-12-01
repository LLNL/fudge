# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the lepton class definitions.
Leptons include electrons, neutrinos, etc.
"""

class Generation :

    electronic = 'electronic'
    muonic = 'muonic'
    tauonic = 'tauonic'

    allowed = ( electronic, muonic, tauonic )

from . import particle as particleModule

class alias( particleModule.alias ) :

    moniker = 'leptonAlias'

class particle( particleModule.particle ) :

    moniker = 'lepton'
    familyOrder = 1
    alias = alias
    generations = { Generation.electronic : 0, Generation.muonic : 1,  Generation.tauonic: 2 }  # Only used for sorting particles.
    chargeOrder = { -1 : 0, 1 : 1, 0 : 2 }                                                      # Only used for sorting particles.

    def __init__( self, id, **kwargs ) :

        particleModule.particle.__init__( self, id )

        if( 'generation' not in kwargs ) : raise ValueError( 'generation attribute not present' )
        generation = kwargs['generation']
        if( generation not in Generation.allowed ) :
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

        _particle = particle( self.id, generation = self.generation )
        self.__copyStandardQuantities( _particle )
        return( _particle )

    def extraXMLAttributes( self ) :

        return( ' generation="%s"' % self.generation )

class suite( particleModule.suite ) :

    moniker = 'leptons'
    particle = particle
