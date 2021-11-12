# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the baryon class definitions.
Baryons include the neutron and proton.
"""

from . import particle as particleModule

class alias( particleModule.alias ) :

    moniker = 'baryonAlias'

class particle( particleModule.particle ) :

    moniker = 'baryon'
    familyOrder = 2
    alias = alias

    def __lt__( self, other ) :

        if( self.familyOrderLessThan( other ) ) : return( True )        # Also checks isinstance of other.
        if( self.familyOrder != other.familyOrder ) : return( False )
        return( self.getMass( 'amu' ) < other.getMass( 'amu' ) )

class suite( particleModule.suite ) :

    moniker = 'baryons'
    particle = particle
