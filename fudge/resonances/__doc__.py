# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the resonances, resolved and unresolved classes.

class Hierarchy:

    * resonances contains one or more of scatteringRadius, resolved and unresolved

    * resolved and unresolved are components, containing one or more resonance forms. Supported forms:

        - resolved resonances may be stored in the BreitWigner or RMatrix containers
        - unresolved average widths and level spacings are stored in a 'tabulatedWidths' container
        - They may also use a (deprecated) 'energyIntervals' container with two or more resonance sections inside

      each of these subsections has a list of resonances and attributes
"""
