#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module imports fudge/GNDS_file.py and calls it for each files in the argument list with
a *True* value for the *show* argument. That is, information about each file will be printed.
"""

import sys

from fudge import GNDS_file as GNDS_fileModule

if( __name__ == '__main__' ) :
    """
    Loop over each input file and print information about them.
    """

    for file in sys.argv[1:] : GNDS_fileModule.type(file, show = True)
