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

args = parser.parse_args()

map1 = mapModule.read(args.map1)
map2 = mapModule.read(args.map2)

def checkMaps(map1, map2):
    '''
    Returns a dictionary with the keys being the list of projectiles in the *map1* and the values being
    a dictionary of tarets in for that projectile. The values in the second dictionary are the list
    of protares in *map2* for that projectile/taret pair.
    '''

    projectiles = {}
    for protare in map1.iterate():
        if protare.projectile not in projectiles:
            projectiles[protare.projectile] = {}
        projectiles[protare.projectile][protare.target] = map2.findAllOf(projectile=protare.projectile, target=protare.target)

    return projectiles

def printMissing(projectiles, map):
    '''
    For each projectile, prints the list of targets whose values are empty.
    '''

    print()
    print(map.fileName)
    for projectile in projectiles:
        targets = projectiles[projectile]
        missingTargets = []
        for target in targets:
            if len(targets[target]) == 0:
                missingTargets.append(target)

        if len(missingTargets) > 0:
            print('    %s as projectile is missing the following targets:' % projectile)
            for target in missingTargets:
                print('        %s' % target)

projectiles1 = checkMaps(map1, map2)
projectiles2 = checkMaps(map2, map1)

printMissing(projectiles1, map2)
printMissing(projectiles2, map1)
