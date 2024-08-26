# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module represents classes that are for internal use.

This module contains the following classes:

    +-------------------+-------------------------------------------------------------------------------------------------------+
    | Class             | Description                                                                                           |
    +===================+=======================================================================================================+
    | Floats            | This class represents several parameters to limit a relative smallness float value.                   |
    +-------------------+-------------------------------------------------------------------------------------------------------+
"""

import sys

class Floats :
    """
    This class represent some basic float parameters used by some xData classes and functions.

    This module contains the following members:

    +-------------------+-------------------------------------------------------------------------------------------------------+
    | Class             | Description                                                                                           |
    +===================+=======================================================================================================+
    | epsilonMultiplier | This is the smallest relative distance between two x-value should that xData should store/generate.   |
    |                   | specified in terms of the system parameters in :py:class:`sys.float_info`.                            |
    +-------------------+-------------------------------------------------------------------------------------------------------+
    | epsilon           | This is the smallest relative distance between two x-value should that xData should store/generate.   |
    +-------------------+-------------------------------------------------------------------------------------------------------+

    """

    epsilonMultiplier = 4
    epsilon = sys.float_info.epsilon * epsilonMultiplier
