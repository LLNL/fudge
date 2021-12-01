# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the charge classes.
Charge is typically given in units of the electron charge 'e'.
"""

from . import quantity as quantityModule

baseUnit = quantityModule.stringToPhysicalUnit( 'e' )

class integer( quantityModule.integer ) :

    baseUnit = baseUnit

class suite( quantityModule.numberSuite ) :

    moniker = 'charge'
    _allowedClasses = [ integer ]
