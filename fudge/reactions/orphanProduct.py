# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the orphanProduct class, for storing orphanProduct reaction data.
"""

from . import base as baseModule

__metaclass__ = type

class orphanProduct( baseModule.base_reaction ):
    """
    An orphanProduct reaction is another special case of the <reaction> class, used to store product information
    for a product not in the exclusive reactions. Currently only used for photons (i.e., gammas) from ENDL's C=55 or the ENDF
    equivalent.
    """

    moniker = 'orphanProduct'
