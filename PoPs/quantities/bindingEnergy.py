# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the binding energy classes.
Binding energy is typically stored in units of eV.
"""

from . import quantity as quantityModule

baseUnit = quantityModule.stringToPhysicalUnit( 'eV' )

class Double( quantityModule.Double ) :

    baseUnit = baseUnit

class Suite( quantityModule.NumberSuite ) :

    moniker = 'bindingEnergy'
    _allowedClasses = [ Double ]
