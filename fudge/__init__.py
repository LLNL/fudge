# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from . import fudgeVersion as versionModule

FUDGE_MAJORVERSION = versionModule.FUDGE_MAJORVERSION
FUDGE_MINORVERSION = versionModule.FUDGE_MINORVERSION

__version__ = f'{versionModule.FUDGE_MAJORVERSION}.{versionModule.FUDGE_MINORVERSION}'
if versionModule.FUDGE_RELEASECANDIDATE != '':
    __version__ += f'{versionModule.FUDGE_RELEASECANDIDATE}post{versionModule.FUDGE_POSTRELEASE}'

from .__doc__ import __doc__
