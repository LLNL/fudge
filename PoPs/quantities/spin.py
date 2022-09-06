# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the spin classes.
Spins are typically stored in units of hbar.
"""

from . import quantity as quantityModule

baseUnit = quantityModule.stringToPhysicalUnit( 'hbar' )

class Fraction( quantityModule.Fraction ) :

    baseUnit = baseUnit

class Suite( quantityModule.NumberSuite ) :

    moniker = 'spin'
    _allowedClasses = [ Fraction ]
