# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from .. import suite as suiteModule
from ..quantities import quantity as quantityModule

class double( quantityModule.double ) :

    pass

class suite( suiteModule.sortedSuite ) :

    moniker = 'rate'

    def __init__( self ) :

        suiteModule.sortedSuite.__init__( self, allowedClasses = ( double, ) )
