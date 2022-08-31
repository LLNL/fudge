# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the mass classes.
Particle masses are typically stored in units like 'amu', 'MeV/c**2' or possibly 'kg'
"""

from . import quantity as quantityModule

baseUnit = quantityModule.stringToPhysicalUnit( 'kg' )

class Double( quantityModule.Double ) :

    baseUnit = baseUnit

class Suite( quantityModule.NumberSuite ) :

    moniker = 'mass'
    _allowedClasses = [ Double ]
