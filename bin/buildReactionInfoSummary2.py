#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import pathlib
import argparse

from fudge import enums as enumsModule
from fudge import map as mapModule
from fudge import reactionSuite as reactionSuiteModule
from fudge import protareProductInfo as protareProductInfoModule

description = '''Builds **ris** (reaction infomation summary) data for each protare in the specified map file and writes the results to a **ris** file.'''

def processMap(map, recursive=False):
    '''
    Creates a **ris** file for the protares in map. The output file is in the same location as the map file but with the extension ".ris".

    :param map:         A Map instance.
    :param recursive:   If **True**, **ris** files for imported map files are also built.
    '''

    input = pathlib.Path(map.path)
    output = input.with_suffix('.ris')
    with output.open('w') as fOut:
        fOut.write('#ris: 1.0\n')
        for entry in map:
            if isinstance(entry, mapModule.Import):
                fOut.write('#import : %s\n' % pathlib.Path(entry.path).with_suffix('.ris'))
                if recursive:
                    processMap(entry.map, recursive=recursive)
            else:
                if isinstance(entry, mapModule.TNSL):
                    continue
                if entry.interaction == enumsModule.Interaction.atomic:
                    continue
                protare = entry.read(lazyParsing=True)
                fOut.write('#protare : %s : %s : %s : %s\n' % (entry.projectile, entry.target, entry.evaluation, protare.domainUnit))

                aliases = [[alias.id, protare.PoPs.final(alias.id).id] for alias in protare.PoPs.aliases if alias.isMetaStable()]
                if len(aliases) > 0:
                    fOut.write('#aliases : %s\n' % len(aliases))
                    for pid, alias in aliases:
                        fOut.write('    %s : %s\n' % (pid, alias))

                fOut.write('#reactions : %s\n' % len(protare.reactions))
                lines = protareProductInfoModule.protareProductInfo2(protare, '    ') + ['']
                fOut.write('\n'.join(lines))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('input',                                                help='The map file whose protares are analyzed.')
    parser.add_argument('--recursive', action='store_true',                     help='If present, imported map files a also processed; otherwise, import map files are not processed.')

    args = parser.parse_args()

    processMap(mapModule.read(args.input, lazyParsing=True), recursive=args.recursive)
