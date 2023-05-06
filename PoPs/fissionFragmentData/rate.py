# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""This module defines classes for storing the characteristic beta-decay rate
for beta-delayed neutron products."""

from .. import suite as suiteModule
from ..quantities import quantity as quantityModule

class Double( quantityModule.Double ) :
    """Stores a characteristic decay rate."""

    pass

# FIXME CMM: specifications say only one 'double' can appear inside 'rate'.
# Should we define suiteModule.OneChildOnly class?
class Suite( suiteModule.SortedSuite ) :
    """Suite containing one or more rates, each corresponding to a unique data **style**."""

    moniker = 'rate'

    def __init__( self ) :

        suiteModule.SortedSuite.__init__( self, allowedClasses = ( Double, ) )
