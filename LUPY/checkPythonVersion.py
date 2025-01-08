# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
if sys.version_info.major != 3:
    raise ValueError('FUDGE only tested with Python 3.')
elif sys.version_info.minor < 7:
    raise ValueError('FUDGE only tested with Python 3 with minor version 7 or higher.')

try:
    import numericalFunctions
except ImportError:
    raise ImportError("Cannot import required modules! Have extensions been compiled? You are using python %s" % sys.version)
