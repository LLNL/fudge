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
from .. import IDs as IDsModule
from .. import intId as intIdModule
from . import particle as particleModule

class Particle( particleModule.Particle ) :

    moniker = 'gaugeBoson'
    familyOrder = 0

    def __lt__( self, other ) :

        if( self.familyOrderLessThan( other ) ) : return( True )        # Also checks isinstance of other.
        if( self.familyOrder != other.familyOrder ) : return( False )
        return( False )                                                 # Currently, only supports photon.

    def intid(self, intidDB={}):

        base, anti, qualifier = miscModule.baseAntiQualifierFromID(self.id)
        gaugeBosonIndex = {IDsModule.photon: 0}.get(base)
        if gaugeBosonIndex is None:
            raise ValueError('Gauge boson "%s" does not have a defined intid.' % (self.id))

        return intIdModule.intidHelper(anti, intIdModule.Family.gaugeBoson, gaugeBosonIndex)

class Suite( particleModule.Suite ) :

    moniker = 'gaugeBosons'
    particle = Particle
