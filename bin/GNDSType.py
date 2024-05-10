#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from fudge import GNDS_file as GNDS_fileModule

summaryDocStringFUDGE = """This script prints the GNDS type of each file listed."""

description = """
This script imports fudge/GNDS_file.py and calls it for each file in the argument list with
a *True* value for the *show* argument. That is, information about the GNDS type of each file is printed.
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('paths', nargs='*',                     help='The list of file whose GNDS type is printed.')
args = parser.parse_args()

if( __name__ == '__main__' ) :
    """
    Loop over each input file and print information about them.
    """

    for path in args.paths:
        GNDS_fileModule.type(path, show=True)
