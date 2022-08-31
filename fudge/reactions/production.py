# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the production class, for storing production cross section of radioactive products.
"""

from . import base as baseModule

class Production(baseModule.Base_reaction):
    """
    Production reactions are another special case of the <reaction> class, used to store the cross section for
    producing specific radioactive products. Requires a cross section and product id, but no distributions.
    As with the summedReaction, 'outputChannel' should be a channels.simpleOutputChannel instance.
    """

    moniker = 'production'
