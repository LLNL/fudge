# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Top-level __init__ for the fudge project.
"""

from . import fudgeVersion as versionModule

FUDGE_MAJORVERSION = versionModule.FUDGE_MAJORVERSION
FUDGE_MINORVERSION = versionModule.FUDGE_MINORVERSION

__version__ = f'{versionModule.FUDGE_MAJORVERSION}.{versionModule.FUDGE_MINORVERSION}'
if versionModule.FUDGE_RELEASECANDIDATE != '':
    __version__ += f'{versionModule.FUDGE_RELEASECANDIDATE}post{versionModule.FUDGE_POSTRELEASE}'

from .__doc__ import __doc__

try:
    import numericalFunctions
except ImportError:
    raise ImportError("Cannot import required modules! Have extensions been compiled?")

from .core import *
from . import vis

# if we want to export a smaller set of files with 'from fudge import *':
#__all__ = ['core',...]
