# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the parity classes.
Parity is unitless, and is stored as an integer: +1 for positive parity, -1 for negative.
"""

from . import quantity as quantityModule

baseUnit = quantityModule.stringToPhysicalUnit( '' )

class Integer( quantityModule.Integer ) :

    baseUnit = baseUnit

    def __init__( self, label, value, unit, documentation = '' ) :

        quantityModule.Integer.__init__( self, label, value, unit, documentation = documentation )
        if( self.value not in [ -1, 1 ] ) : raise ValueError( 'parity value %d not allowed' % self.value )

class Suite( quantityModule.Suite ) :

    moniker = 'parity'
    _allowedClasses = [ Integer ]
