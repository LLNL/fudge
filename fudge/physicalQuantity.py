# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import physicalQuantity as physicalQuantityModule

class temperature( physicalQuantityModule.physicalQuantity ) :

    moniker = 'temperature'

class  U( physicalQuantityModule.physicalQuantity ) :

    moniker = 'U'

class  EFL( physicalQuantityModule.physicalQuantity ) :

    moniker = 'EFL'

class  EFH( physicalQuantityModule.physicalQuantity ) :

    moniker = 'EFH'

class mass( physicalQuantityModule.physicalQuantity ) :

    moniker = 'mass'
