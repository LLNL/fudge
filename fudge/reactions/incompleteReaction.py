# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the incompleteReaction class, for storing partially evaluated reactions
that do not have enough information to be used for transport but may still contain useful info.

Originally designed to store sub-actinide fission and standards cross-sections
"""

from . import base as baseModule

__metaclass__ = type

class incompleteReaction( baseModule.base_reaction ):
    """
    Similar to the base <reaction> class, except it generally does not contain
    a complete description of products, and cannot be used for transport. It must contain a crossSection
    and outputChannel, but the list of products may be empty.
    """

    moniker = 'reaction'
