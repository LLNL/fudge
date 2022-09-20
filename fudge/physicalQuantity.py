# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import physicalQuantity as physicalQuantityModule

class Temperature( physicalQuantityModule.PhysicalQuantity ) :

    moniker = 'temperature'

class U( physicalQuantityModule.PhysicalQuantity ) :

    moniker = 'U'

class EFL( physicalQuantityModule.PhysicalQuantity ) :

    moniker = 'EFL'

class EFH( physicalQuantityModule.PhysicalQuantity ) :

    moniker = 'EFH'

class Mass( physicalQuantityModule.PhysicalQuantity ) :

    moniker = 'mass'
