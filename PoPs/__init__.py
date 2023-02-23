# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import checkPythonVerions

POPS_MAJORVERSION = 0
POPS_MINORVERSION = 1
POPS_PATCHLEVEL = 0
__version__ = '%s.%s.%s' % ( POPS_MAJORVERSION, POPS_MINORVERSION, POPS_PATCHLEVEL )

from .__doc__ import __doc__
