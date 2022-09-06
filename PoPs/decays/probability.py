# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the decay probability classes.
"""

from ..quantities import quantity as quantityModule

baseUnit = quantityModule.stringToPhysicalUnit( '' )

class Double( quantityModule.Double ) :

    __baseUnit = baseUnit

    def __init__( self, label, value, unit = None, documentation = '' ) :

        if( unit is None ) : unit = quantityModule.stringToPhysicalUnit( '' )
        quantityModule.Double.__init__( self, label, value, unit, documentation = documentation )

class Suite( quantityModule.NumberSuite ) :

    moniker = 'probability'
    _allowedClasses = [ Double ]
