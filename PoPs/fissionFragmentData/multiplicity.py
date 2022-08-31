# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from .. import suite as suiteModule
from ..quantities import quantity as quantityModule

class Double( quantityModule.Double ) :

    pass

class Suite( suiteModule.SortedSuite ) :

    moniker = 'multiplicity'

    def __init__( self ) :

        suiteModule.SortedSuite.__init__( self, allowedClasses = ( Double, ) )
