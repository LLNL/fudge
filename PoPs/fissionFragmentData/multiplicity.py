# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""This module defines the Suite and Double classes, used to store
spontaneous fission delayed neutron multiplicities."""

from .. import suite as suiteModule
from ..quantities import quantity as quantityModule

class Double( quantityModule.Double ) :
    """Average delayed neutron multiplicity, stored as a double."""

    pass

class Suite( suiteModule.SortedSuite ) :
    """This class stores one or more multiplicity values, e.g. corresponding
    to the 'evaluated' data style."""

    moniker = 'multiplicity'

    def __init__( self ) :

        suiteModule.SortedSuite.__init__( self, allowedClasses = ( Double, ) )
