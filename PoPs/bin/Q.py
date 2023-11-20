#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import userFUDGE_defaults

summaryDocStringPoPs = '''Calculates the Q-value for the list of ingoing and outgoing particles.'''

description = '''
    Calculates the Q-value for the list of ingoing and outgoing particles using data in a GNDS PoPs database.
    If only ingoing particles are listed, the output is the mass-energy of those particles. If only one ingoing
    particle is specified, it is treated as decaying to the outgoing particles. To get unit in "amu", enter 
    "--unit 'amu*c**2'" as unit must be that of energy. Note, all ingoing particles 
    must be entered before the "-o" option is entered and only one "-o" option should be entered as all 
    particles entered after the "-o" option are added to the "-o" list.
    '''

import pathlib
import argparse

from PoPs import database as databaseModule
from PoPs import misc as miscModule
from PoPs.chemicalElements import misc as chemicalElementsMiscModule

parser = argparse.ArgumentParser(description=description)
userFUDGE_defaults.add_argument(parser, 'PoPs_path', True, '--pops',                help='GNDS PoPs file to use.')
parser.add_argument('ingoingParticles', nargs='+',                                  help='List of ingoing particles.')
parser.add_argument('-o', '--outgoingParticles', nargs='*', default=[],             help='List of outgoing particles.')
parser.add_argument('-u', '--unit', default='MeV',                                  help='Energy unit.')
parser.add_argument('--balanceZA', action='store_true',                             help='If present, this script will calculate the heavy residual (outgoing particle) to balance Z and A.')

def compoundZandA(pops, particleList):
    '''
    Calculates the totol Z and A for the list particles in *particleList*.
    '''

    Z = 0
    A = 0

    particleIds = miscModule.parseParticleList(particleList)

    for particleId in particleIds:
        particle = pops[particleId]
        Zp, Ap, dummy, dummy = chemicalElementsMiscModule.ZAInfo(particle)
        Z += Zp
        A += Ap

    return Z, A

if __name__ == '__main__':
    args = parser.parse_args()

    ingoingParticles = args.ingoingParticles
    outgoingParticles = args.outgoingParticles

    pops = databaseModule.read(args.pops)
    ingoingParticles  = args.ingoingParticles
    outgoingParticles = args.outgoingParticles

    if args.balanceZA:
        Zi, Ai = compoundZandA(pops, ingoingParticles)
        Zo, Ao = compoundZandA(pops, outgoingParticles)
        Z = Zi - Zo
        A = Ai - Ao
        if Z > 0:
            outgoingParticles.append(chemicalElementsMiscModule.idFromZAndA(Z, A))

    if not outgoingParticles:
        comment = '# mass energy for particle(s) %s.' % ' + '.join(ingoingParticles)
    elif len(ingoingParticles) == 1:
        comment = '# Q-value for decay "%s -> %s".' % (ingoingParticles[0], ' + '.join(outgoingParticles))
    else:
        comment = '# Q-value for reaction "%s -> %s".' % (' + '.join(ingoingParticles), ' + '.join(outgoingParticles))

    print(pops.QValue(ingoingParticles, outgoingParticles, unit=args.unit), args.unit, comment)
