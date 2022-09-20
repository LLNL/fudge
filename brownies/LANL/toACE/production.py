# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module adds the method toACE to the production class.
"""

from fudge.reactions import production as productionModule

def toACE( self, temperature, EMin, data, verbose ) :

    MT = self.ENDF_MT
    if( verbose > 1 ) : print( '   %s: MT = %d: Skipping production reaction' % ( str( self ), MT ) )

productionModule.Production.toACE = toACE
