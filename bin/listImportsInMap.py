#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from fudge import map as mapModule

description = '''Reads the specified map file and recursively prints the path to all map files imported.'''

parser = argparse.ArgumentParser(description=description)

parser.add_argument('map',                                      help='''The path to a map file.''')

args = parser.parse_args( )

def listImports(mapFile, indent=''):

    print('%s%s' % (indent, mapFile.path))
    for entry in mapFile:
        if isinstance(entry, mapModule.Import):
            subMap = entry.readMap()
            listImports(subMap, indent + '    ')

listImports(mapModule.Map.readXML(args.map))
