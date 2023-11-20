#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

summaryDocStringPoPs = '''Given a particle id, prints its intid.'''

description = '''Given a particle id, this script prints its intid. This script requires a GNDS PoPs file as its first argument and the 
PoPs file must contain the requested id.
'''

import argparse
import pathlib

from PoPs import database as databaseModule

parser = argparse.ArgumentParser(description=description)
parser.add_argument('pops', type=pathlib.Path,      help='PoPs file that must contain the requested id.')
parser.add_argument('ids', nargs='*',               help='The list of ids whose intid is printed.')

if __name__ == '__main__':
    args = parser.parse_args()

    pops = databaseModule.read(args.pops)

    for id in args.ids:
        particle = pops[id]
        print('%-16s: %s' % (id, particle.intid()))
