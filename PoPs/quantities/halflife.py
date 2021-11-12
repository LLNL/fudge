# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the halflife classes.
Halflife may be given as a value with units of time (i.e. '4.523 ms'),
or as one of the strings 'stable' or 'unstable'
"""

STABLE = 'stable'
UNSTABLE = 'unstable'

from . import quantity as quantityModule

baseUnit = quantityModule.stringToPhysicalUnit( 's' )

class double( quantityModule.double ) :

    baseUnit = baseUnit

class string( quantityModule.string ) :

    baseUnit = baseUnit

class suite( quantityModule.numberSuite ) :

    moniker = 'halflife'
    _allowedClasses = [ double, string ]
