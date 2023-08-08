#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from fudge import map as mapModule

summaryDocStringFUDGE = '''Prints the protare differences between two map files.'''

description = '''For each projectile in the two map files, prints targets that are in one map file and not in the other.'''

parser = argparse.ArgumentParser(description=description)
parser.add_argument('map1', type=str,                           help='A map file.')
parser.add_argument('map2', type=str,                           help='Another map file.')
parser.add_argument('--dump', action='store', default='',       help='If present, projectile/target information for each map file is dumped to a file in the current directory with value prefixing each files name.')

args = parser.parse_args()

map1 = mapModule.read(args.map1)
map2 = mapModule.read(args.map2)

def checkMaps(map1):
    '''
    Returns a dictionary with the keys being the list of projectiles in the *map1* and the values being
    a dictionary of tarets in for that projectile. The values in the second dictionary are the list
    of protares in *map2* for that projectile/taret pair.
    '''

    projectiles = {}
    for protare in map1.iterate():
        if protare.projectile not in projectiles:
            projectiles[protare.projectile] = {}
        if protare.target not in projectiles[protare.projectile]:
            projectiles[protare.projectile][protare.target] = []
        projectiles[protare.projectile][protare.target].append(protare)

    return projectiles

def printMissing(projectiles1, projectiles2, map, index):
    '''
    For each projectile, prints the list of targets whose values are empty.
    '''

    print()
    print(map.fileName)

    fOut = None
    if args.dump != '':
        fOut = open('%s_%s.dat' % (args.dump, index), 'w')

    for projectile in sorted(projectiles1):
        targets1 = projectiles1[projectile]

        if fOut is not None:
            print('%s - %s' % (projectile, index), file=fOut)
            for target in sorted(targets1):
                print('    %s' % target, file=fOut)

        if projectile not in projectiles2:
            print('    Projectile "%s" missing from map file.' % projectile)
        else:
            targets2 = projectiles2[projectile]
            missingTargets = []
            for target in sorted(targets1):
                if target not in targets2:
                    missingTargets.append(target)

            if len(missingTargets) > 0:
                print('    %s as projectile is missing the following targets:' % projectile)
                for target in missingTargets:
                    print('        %s' % target)

projectiles1 = checkMaps(map1)
projectiles2 = checkMaps(map2)

printMissing(projectiles2, projectiles1, map1, 2)
printMissing(projectiles1, projectiles2, map2, 1)
