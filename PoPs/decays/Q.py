# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the decay Q classes.
"""

from ..quantities import quantity as quantityModule

baseUnit = quantityModule.stringToPhysicalUnit( 'eV' )

class double( quantityModule.double ) :

    __baseUnit = baseUnit

class suite( quantityModule.numberSuite ) :

    moniker = 'Q'
    _allowedClasses = [ double ]
