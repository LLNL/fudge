# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__doc__ = """The PoPs.fissionFragmentData module defines classes for storing information about prompt
and delayed products from spontaneous fission. These containers are similar to classes for induced fission,
except with no dependence on projectile energy.

Data hierarchy:
    ```
    - fissionFragmentData
        - delayedNeutrons
            - delayedNeutron (1 or more)
                - rate
                - product
        - fissionEnergyReleased
            (defines components of fission energy release, currently only used for induced fission)
        - productYields
            - productYield  (for spontaneous fission)
                - nuclides
                - elapsedTimes
                    - elapsedTime
                        - time
                        - yields
                            - values
                            - uncertainty

            - productYield (for induced fission)
                - nuclides (may appear at this level or below)
                - elapsedTimes
                    - elapsedTime
                        - time
                        - incidentEnergies
                            - incidentEnergy
                                - energy
                                - yields
                                    - nuclides (overrides definition above)
                                    - values
                                    - uncertainty```
"""
