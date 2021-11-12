# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the gaugeBoson class definitions.
Gauge bosons include the photon, gluon, W and Z particles
"""

from . import particle as particleModule

class alias( particleModule.alias ) :

    moniker = 'gaugeBosonAlias'

class particle( particleModule.particle ) :

    moniker = 'gaugeBoson'
    familyOrder = 0
    alias = alias

    def __lt__( self, other ) :

        if( self.familyOrderLessThan( other ) ) : return( True )        # Also checks isinstance of other.
        if( self.familyOrder != other.familyOrder ) : return( False )
        return( False )                                                 # Currently, only supports photon.

class suite( particleModule.suite ) :

    moniker = 'gaugeBosons'
    particle = particle
