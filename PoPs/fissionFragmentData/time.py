# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines classes for storing the elapsed time since fission.
Each 'duration' instance has an associated time, where time = 0 for prompt fission products and time > 0 for delayed. Time may also be given as a string, e.g. 'Unspecified'.
"""

from .. import suite as suiteModule
from ..quantities import quantity as quantityModule

class Double( quantityModule.Double ) :

    pass

class String( quantityModule.String ) :

    pass

class Suite( suiteModule.SortedSuite ) :

    moniker = 'time'

    def __init__( self ) :

        suiteModule.SortedSuite.__init__( self, allowedClasses = ( Double, String ) )
