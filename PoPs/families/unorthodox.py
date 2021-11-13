# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the 'unorthodox' class definitions.
Unorthodox particles are used when a particle-like object needs to be described in PoPs even though
it does not fit in the normal definition of a particle.
Currently used to store average fission product particles
(these are computed from a weighted average over actual fission product yields).
In the future, may also be used for 'thermal scattering law' particles.
"""

from . import particle as particleModule

class alias( particleModule.alias ) :

    moniker = 'unorthodoxAlias'

class particle( particleModule.particle ) :

    moniker = 'unorthodox'
    familyOrder = 5
    alias = alias

    def __lt__( self, other ) :

        if( self.familyOrderLessThan( other ) ) : return( True )        # Also checks isinstance of other.
        if( self.familyOrder != other.familyOrder ) : return( False )
        return( self.id < other.id )

class suite( particleModule.suite ) :

    moniker = 'unorthodoxes'
    particle = particle