# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the fissionComponent class
"""

from . import base as baseModule


class FissionComponent(baseModule.Base_reaction1):
    """
    Special case of reaction, for storing 1st, 2nd, 3rd, etc fission chances when we also have total.
    """

    moniker = 'fissionComponent'
