# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the decay Q classes.
"""

from ..quantities import quantity as quantityModule

baseUnit = quantityModule.stringToPhysicalUnit( 'eV' )

class Double( quantityModule.Double ) :

    __baseUnit = baseUnit

class Suite( quantityModule.NumberSuite ) :

    moniker = 'Q'
    _allowedClasses = [ Double ]
