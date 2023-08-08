#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

description = '''
    Calculates the Q-value for the list of input and output particles using data in a GNDS PoPs database.
    Note, the energy mass of a particle can be obtained as "Q.py pops.xml -n 1 id" where id is the particle's PoPs id.
    '''

import pathlib
import argparse

from PoPs import database as databaseModule

parser = argparse.ArgumentParser(description=description)
parser.add_argument('pops', type=pathlib.Path,                      help='GNDS PoPs file to use.')
parser.add_argument('particles', nargs='+',                         help='List of input and output particles.')
parser.add_argument('-n', '--numberOfInputParticles', action='store', default=2, type=int,
                                                                    help='Number of input particles. Default is 2.')
if __name__ == '__main__':
    args = parser.parse_args()

    ingoingParticles  = args.particles[:args.numberOfInputParticles]
    outgoingParticles = args.particles[args.numberOfInputParticles:]
    pops = databaseModule.read(args.pops)
    print(pops.QValue(ingoingParticles, outgoingParticles, unit='MeV'), 'MeV      # Q-value for reacton "%s -> %s".' % (' + '.join(ingoingParticles), ' + '.join(outgoingParticles)))
