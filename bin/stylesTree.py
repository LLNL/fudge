#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
import pathlib

from fudge import GNDS_file as GNDS_fileModule
from fudge import styles as stylesModule
from fudge import reactionSuite as reactionSuiteModule

summaryDocStringFUDGE = """Prints a tree representation of the styles in each specified GNDS reactionSuite."""

description = """
Prints a tree representation of the styles in each specified GNDS reactionSuite.
"""

parser = argparse.ArgumentParser(description=description)

parser.add_argument('files', nargs="*",                         help='The list of GNDS reactionSuite files to process.')
args = parser.parse_args()

def printLimb(limb, prefix):
    """
    This function print a limb of the tree.

    :param limb:    The current limb of the tree to print.
    :param prefix:  The prefix for the limb.
    """

    if len(limb) != 0:
        prefix2 = prefix + '  |   '
        if len(prefix) > 4:
            prefix = prefix[:-3] + '-->'
        for label in limb:
            print('%s%s' % (prefix, label))
            printLimb(limb[label], prefix2)

for file in args.files:
    GNDS = GNDS_fileModule.preview(file)
    if isinstance(GNDS, reactionSuiteModule.ReactionSuite):
        print(file)
        chains = GNDS.styles.chains()
        limb= {}
        for chain in chains:
            style = chain[0]
            subLimb = limb
            for style in reversed(chain):
                if style.label not in subLimb:
                    subLimb[style.label] = {}
                subLimb = subLimb[style.label]
        printLimb(limb, '    ')
    else:
        continue
