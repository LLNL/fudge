# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the decay probability classes.
"""

from ..quantities import quantity as quantityModule

baseUnit = quantityModule.stringToPhysicalUnit( '' )

class double( quantityModule.double ) :

    __baseUnit = baseUnit

    def __init__( self, label, value, unit = None, documentation = '' ) :

        if( unit is None ) : unit = quantityModule.stringToPhysicalUnit( '' )
        quantityModule.double.__init__( self, label, value, unit, documentation = documentation )

class suite( quantityModule.numberSuite ) :

    moniker = 'probability'
    _allowedClasses = [ double ]
