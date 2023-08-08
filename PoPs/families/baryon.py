# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the baryon class definitions.
Baryons include the neutron and proton.
"""

from .. import misc as miscModule
from . import particle as particleModule

class Alias( particleModule.Alias ) :

    moniker = 'baryonAlias'

class Particle( particleModule.Particle ) :

    moniker = 'baryon'
    familyOrder = 2
    alias = Alias

    def __lt__( self, other ) :

        if( self.familyOrderLessThan( other ) ) : return( True )        # Also checks isinstance of other.
        if( self.familyOrder != other.familyOrder ) : return( False )
        return( self.getMass( 'amu' ) < other.getMass( 'amu' ) )

    def intid(self, intidDB={}):

        sign = -1 if self.isAnti else 1
        base, anti, qualifier = miscModule.baseAntiQualifierFromID(self.id)
        baryonIndex = {'n': 0, 'p': 1}.get(base)
        if baryonIndex is None:
            ValueError('Baryon "%s" does not have a defined intid.' % (self.id))

        return sign * (2**30 + 2 * 2**26 + baryonIndex)

class Suite( particleModule.Suite ) :

    moniker = 'baryons'
    particle = Particle
