# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
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

class double( quantityModule.double ) :

    baseUnit = baseUnit

class suite( quantityModule.numberSuite ) :

    moniker = 'bindingEnergy'
    _allowedClasses = [ double ]
