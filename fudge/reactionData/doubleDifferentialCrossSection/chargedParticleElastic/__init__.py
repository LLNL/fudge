# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains a form for storing the differential cross section for charged-particle elastic scattering.
Internally, data can be represented three ways:

    - pure Rutherford scattering
    - Rutherford scattering along with Legendre expansions for nuclear scattering
        and for real and imaginary nuclear/Coulomb interference
    - Rutherford scattering along with effective cross sections and distributions,
        obtained by summing the nuclear and interference terms
"""

from . import CoulombPlusNuclearElastic
from . import nuclearAmplitudeExpansion
from . import nuclearPlusInterference
