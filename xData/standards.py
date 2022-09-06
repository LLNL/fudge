# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains standard xData types.
"""

import sys

class Floats :

    epsilonMultiplier = 4
    epsilon = sys.float_info.epsilon * epsilonMultiplier
