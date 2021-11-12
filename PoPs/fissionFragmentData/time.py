# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines classes for storing the elapsed time since fission.
Each 'duration' instance has an associated time, where time = 0 for prompt fission products and time > 0 for delayed
"""

from .. import suite as suiteModule
from ..quantities import quantity as quantityModule

class double( quantityModule.double ) :

    pass

class string( quantityModule.string ) :

    pass

class suite( suiteModule.sortedSuite ) :

    moniker = 'time'

    def __init__( self ) :

        suiteModule.sortedSuite.__init__( self, allowedClasses = ( double, string ) )
