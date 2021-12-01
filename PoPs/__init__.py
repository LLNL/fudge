# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Top-level __init__ for the PoPs package.
"""

POPS_MAJORVERSION = 0
POPS_MINORVERSION = 1
POPS_PATCHLEVEL = 0
__version__ = '%s.%s.%s' % ( POPS_MAJORVERSION, POPS_MINORVERSION, POPS_PATCHLEVEL )

from .__doc__ import __doc__
