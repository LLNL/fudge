#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from PoPs import database as databaseModule
from PoPs.families import nucleus as nucleusModule

summaryDocString__PoPs = '''Simple diff of two GNDS PoPs files.'''

description = '''
Diffs two GNDS PoPs files. Currently, only prints:

    * particles missing from one file that is in the other,
    * if option "--mass" present, then for particles in each database, their mass difference if not 0.

Note, mass unit cannot be an energy unit. For example use '--massUnit "MeV/c**2"' and not '--massUnit MeV'.
'''

massUnit = 'amu'
missingPerLine = 10
massRelDiff = 1e-8
maxIdLength = 6             # Set to correct value later.
idFormat = ''               # Set to correct value later.

parser = argparse.ArgumentParser(description=description)
parser.add_argument('pops1', type=str,                                      help='First GNDS PoPs file for diffing.')
parser.add_argument('pops2', type=str,                                      help='Second GNDS PoPs file for diffing.')
parser.add_argument('--missingPerLine', type=int, default=missingPerLine,   help='Number of missing ids to print per line. Default is %s.' % missingPerLine)
parser.add_argument('--mass', action='store_true',                          help='If present, masses of common particles are diffed.')
parser.add_argument('--massRelDiff', type=float, default=massRelDiff,       help='If option "--mass" present, only print mass difference if relative difference greater than this value. Default is "%s".' % massRelDiff)
parser.add_argument('--massUnit', default=massUnit,                         help='Unit for masses. Default is "%s".' % massUnit)
parser.add_argument('--ignoreNuclei', action='store_true',                  help='If present, nuclei diffing will not be done.')

def missingIds(pops, args, set1, set2, first):
    '''Prints the list of particles in *set1* but not in *set2*.'''

    s1, s2 = 'first', 'second'
    if not first:
        s1, s2, = s2, s1
    print('The following particles are in the %s PoPs file but not the %s PoPs file:' % (s1, s2))
    diffSet = set1.difference(set2)
    if args.ignoreNuclei:
        diffList = []
        for pid in diffSet:
            if not isinstance(pops[pid], nucleusModule.Particle):
                diffList.append(pid)
    else:
        diffList = list(diffSet)
    if len(diffList) > 0:
        idFormat = '%%-%ds' % max([len(pid) for pid in diffList])
        diff = sorted([[pid.lower(), pid] for pid in diffList])
        pids = [idFormat % pid for _, pid in diff]
        for index in range(0, len(pids), args.missingPerLine):
            print('    %s' % ' '.join(pids[index:index+args.missingPerLine]))

if __name__ == '__main__':
    args = parser.parse_args()

    massUnit = args.massUnit
    massRelDiff = args.massRelDiff

    pops1 = databaseModule.read(args.pops1)
    pops2 = databaseModule.read(args.pops2)

    set1 = set(particle.id for particle in pops1)
    set2 = set(particle.id for particle in pops2)

    maxIdLength = max(max([len(pid) for pid in set1]), max([len(pid) for pid in set2]))
    idFormat = '%%-%ds' % maxIdLength

    missingIds(pops1, args, set1, set2, True)
    missingIds(pops2, args, set2, set1, False)

    intersection = sorted([[pid.lower(), pid] for pid in set1.intersection(set2)])
    intersection = list(pid for _, pid in intersection)
    if args.mass:
        print('Mass difference (%s):' % massUnit)
        header = '    %-12s %-20s %-20s  %-10s %-10s' % ('id', 'mass-1', 'mass-2', 'diff', 'rel-diff')
        print(header)
        print(len(header) * '=')
        for pid in intersection:
            particle1 = pops1[pid]
            if args.ignoreNuclei and isinstance(particle1, nucleusModule.Particle):
                continue
            particle2 = pops2[pid]
            mass1 = particle1.getMass(massUnit)
            mass2 = particle2.getMass(massUnit)
            norm = max(mass1, mass2)
            diff = mass1 - mass2
            if abs(diff) > norm * massRelDiff:
                print('    %-12s %-20s %-20s % 10.2e %10.2e' % (pid, mass1, mass2, diff, diff / norm))
