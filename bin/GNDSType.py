#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module import the LUPY/GNDSType.py and calls it for each files in the argument list with 
a *True* value for the *show* argument. That is, information about each files will be printed.
"""

import sys

from LUPY import GNDSType as GNDSTypeModule

if( __name__ == '__main__' ) :
    """
    Loop over each input file and print information about them.
    """

    for file in sys.argv[1:] : GNDSTypeModule.type( file, show = True )
