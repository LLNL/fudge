# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the gaugeBoson class definitions.
Gauge bosons include the photon, gluon, W and Z particles
"""

from .. import misc as miscModule
from . import particle as particleModule

class Alias( particleModule.Alias ) :

    moniker = 'gaugeBosonAlias'

class Particle( particleModule.Particle ) :

    moniker = 'gaugeBoson'
    familyOrder = 0
    alias = Alias

    def __lt__( self, other ) :

        if( self.familyOrderLessThan( other ) ) : return( True )        # Also checks isinstance of other.
        if( self.familyOrder != other.familyOrder ) : return( False )
        return( False )                                                 # Currently, only supports photon.

    def intid(self, intidDB={}):

        base, anti, qualifier = miscModule.baseAntiQualifierFromID(self.id)
        baryonIndex = {'photon': 0}.get(base)
        if baryonIndex is None:
            ValueError('Gauge boson "%s" does not have a defined intid.' % (self.id))

        return 2**30

class Suite( particleModule.Suite ) :

    moniker = 'gaugeBosons'
    particle = Particle
