# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
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

class fraction( quantityModule.fraction ) :

    baseUnit = baseUnit

class suite( quantityModule.numberSuite ) :

    moniker = 'spin'
    _allowedClasses = [ fraction ]
