#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import glob
import argparse
import xml.sax

from fudge import map as mapModule, GNDS_file as GNDS_fileModule

summaryDocString__FUDGE = """This script examines all map files in the specified directory and prints the minimun list of map files that covers all map files in the specified directory."""

description = """
This script examines all map files in the specified directory and prints the minimun list of map files that
covers all map files in the specified directory. No map files in sub-directories are examined. This script
can be used with the script checkMap.py to check that all the map files in the specified directory are valid.
"""

parser = argparse.ArgumentParser( description = description, formatter_class = argparse.RawTextHelpFormatter )
parser.add_argument('directory', type = str,                 help = 'The directory whose map files are examined.')

args = parser.parse_args( )

directoryRealPath = os.path.realpath(args.directory)

def checkMap(mapFile, uniqueMapFiles):

    map = mapModule.Map.readXML_file(mapFile)
    for entry in map :
        if( isinstance( entry, mapModule.Import ) ) :
            realpath = os.path.realpath(entry.fileName)
            if directoryRealPath == os.path.dirname(realpath):
                try:
                    index = uniqueMapFiles.index(os.path.basename(realpath))
                    uniqueMapFiles.pop(index)
                except:
                    pass
            checkMap(entry.fileName, uniqueMapFiles)

files = glob.glob(os.path.join(args.directory, '*'))
mapFiles = []
for file in files:
    if os.path.isfile(file) :
        try:
            moniker, data = GNDS_fileModule.type(file)
        except xml.sax._exceptions.SAXParseException:
            continue
        if moniker == mapModule.Map.moniker:
            mapFiles.append(file)

uniqueMapFiles = [os.path.basename(mapFile) for mapFile in mapFiles]
for mapFile in mapFiles: checkMap(mapFile, uniqueMapFiles)

for uniqueMapFile in uniqueMapFiles: print(' %s' % uniqueMapFile, end='')
print()
