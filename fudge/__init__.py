# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

"""
top-level __init__ for the fudge project
cmattoon, 3/18/2011
"""
from __future__ import division, nested_scopes

FUDGE_MAJORVERSION = 4
FUDGE_MINORVERSION = 1
FUDGE_PATCHLEVEL = 0
__version__ = '%s.%s.%s' % ( FUDGE_MAJORVERSION, FUDGE_MINORVERSION, FUDGE_PATCHLEVEL )

__docformat__ = 'epytext en'
from __doc__ import __doc__

import fudgeDefaults
import fudgeParameters
from core import *
import vis
import gnd
import particles

# if we want to export a smaller set of files with 'from fudge import *':
#__all__ = ['core','gnd',...]
