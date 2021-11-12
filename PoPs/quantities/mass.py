# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
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

class double( quantityModule.double ) :

    baseUnit = baseUnit

class suite( quantityModule.numberSuite ) :

    moniker = 'mass'
    _allowedClasses = [ double ]
