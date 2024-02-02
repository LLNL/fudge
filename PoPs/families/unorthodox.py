# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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

from .. import IDs as IDsModule
from .. import intId as intIdModule

from . import particle as particleModule

class Particle( particleModule.Particle ) :

    moniker = 'unorthodox'
    familyOrder = 5

    def __lt__( self, other ) :

        if( self.familyOrderLessThan( other ) ) : return( True )        # Also checks isinstance of other.
        if( self.familyOrder != other.familyOrder ) : return( False )
        return( self.id < other.id )

    def intid(self, intidDB={}):

        if self.id == IDsModule.FissionProductENDL99120:
            unorthodoxIndex = 99120
        elif self.id == IDsModule.FissionProductENDL99125:
            unorthodoxIndex = 99125
        else:
            unorthodoxIndex = intidDB.get(self.id)
            if unorthodoxIndex is None:
                unorthodoxIndex = -1
            return unorthodoxIndex

        return intIdModule.intidHelper(False, intIdModule.Family.ENDL_fissionProduct, unorthodoxIndex)            # Only for FissionProductENDL9912[05].

class Suite( particleModule.Suite ) :

    moniker = 'unorthodoxes'
    particle = Particle
