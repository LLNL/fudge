#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from PoPs import database as databaseModule
from PoPs import alias as aliasModule
from PoPs.families import nuclide as nuclideModule
from PoPs.families import nucleus as nucleusModule
from PoPs.families import unorthodox as unorthodoxModule

summaryDocString__PoPs = '''This script does some rudimentary checking of a GNDS PoPs file.'''

description = '''Checks a GNDS PoPs file for the following:

    * Any undefined masses,
    * Any undefined charges.
    * Any nucleus with an undefined energy.
'''

args = None                     # Set properly later.
missingPerLine = 10

parser = argparse.ArgumentParser(description=description)
parser.add_argument('pops', type=str,                                       help='GNDS PoPs file to check.')
parser.add_argument('--missingPerLine', type=int, default=missingPerLine,   
                        help='Number of missing ids to print per line. Default is %s.' % missingPerLine)

def printMissing(message, missing):
    '''Prints *message* and then the list of particles with undefined data per the message.'''

    print(message)
    if len(missing) > 0:
        idFormat = '%%-%ds' % max([len(pid) for pid in missing])
        diff = sorted([[pid.lower(), pid] for pid in missing])
        pids = [idFormat % pid for dummy, pid in diff]
        for index in range(0, len(pids), args.missingPerLine):
            print('    %s' % ' '.join(pids[index:index+args.missingPerLine]))

if __name__ == '__main__':
    args = parser.parse_args()

    pops = databaseModule.read(args.pops)

    missingChargeNodes = []
    for particle in pops:
        if isinstance(particle, (aliasModule.BaseAlias, unorthodoxModule.Particle)):
            continue
        if len(particle.charge) == 0:
            missingChargeNodes.append(particle.id)
    printMissing('Particles with undefined charges:', missingChargeNodes)

    missingMassNodes = []
    for particle in pops:
        if isinstance(particle, (aliasModule.BaseAlias, nucleusModule.Particle)):
            continue
        if isinstance(particle, nuclideModule.Particle):
            if particle.index > 0:
                continue
        if len(particle.mass) == 0:
            missingMassNodes.append(particle.id)
    printMissing('Particles with undefined masses:', missingMassNodes)

    missingEnergyNodes = []
    for particle in pops:
        if not isinstance(particle, nucleusModule.Particle):
            continue
        if len(particle.energy) == 0:
            missingEnergyNodes.append(particle.id)
    printMissing('Nuclei with undefined energies:', missingEnergyNodes)
