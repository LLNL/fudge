#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

summaryDocString__PoPs = """Translates a intid (integer id) into a GNDS PoPs id."""

description = """
This scripts translates an intid (integer id) into a GNDS PoPs id. An intid is a unique mapping 
of a GNDS PoPs id, which is a string, into an integer id which is used by MCGIDI.
"""

import argparse

from PoPs import IDs as IDsModule
from PoPs import intId as intIdModule
from PoPs import specialNuclearParticleID as specialNuclearParticleIDModule

parser = argparse.ArgumentParser(description=description)
parser.add_argument('intids', type=int, nargs='*',          help="""The list particle intids to parse.""")

if __name__ == '__main__':
    args = parser.parse_args()

    intidDB = intIdModule.IntidDB()

    for intid in args.intids:
        print()
        print('intid:', intid)

        id1 = intIdModule.idFromIntid(intid, intidDB, mode=specialNuclearParticleIDModule.Mode.familiar)
        id2 = intIdModule.idFromIntid(intid, intidDB, mode=specialNuclearParticleIDModule.Mode.nuclide)
        id3 = intIdModule.idFromIntid(intid, intidDB, mode=specialNuclearParticleIDModule.Mode.nucleus)
        if id1 == id2:
            print('  id:', id1)
        elif id2 == IDsModule.photon:
            print('             id:', id2)
            print('    familiar id:', id1)
        else:
            print('    familiar id:', id1)
            print('    nuclide  id:', id2)
            print('    nucleus  id:', id3)

        family, isAnti, *others = intIdModule.parseIntid(intid)
        print()
        print('    Parts:')
        print('       family:', family)
        print('       isAnti:', isAnti)
        print('       others:', *others)
