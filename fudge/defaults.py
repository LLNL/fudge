# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Define user-configurable defaults.
"""

# numTasks: how many parallel jobs to spawn when doing CPU-heavy tasks like resonance reconstruction,
# URR probability table generation, etc.
# Default behavior is to use all available CPUs
import multiprocessing
numTasks = multiprocessing.cpu_count()