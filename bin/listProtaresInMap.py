#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from fudge import map as mapModule

summaryDocStringFUDGE = """This script reads the specified map file, and prints all protares matching the specified projectile, target and evaluation."""

description = '''
This script reads the specified map file, and prints all protares matching the specified projectile, 
target and evaluation. This script recursively scans all imported map files. If a target 
for any projectile is desired, specify the projectile as a blank string (e.g., ''). A 
specified projectile, target or evaluation can be a Python re expression.

Some examples are:

listProtaresInMap.py ENDF-VIII.0/all.map '' '.*s1.*4.*'         # Blank string for projectile selects all projectiles.
n   Cs134   ENDF/B-8.0
n   Os184   ENDF/B-8.0

listProtaresInMap.py all.map '' '.*Os.*4.*'
n     Os184   endl2009.4-rc1-1.0
H1    Os184   endl2009.4-rc1-1.0
H2    Os184   endl2009.4-rc1-1.0
H3    Os184   endl2009.4-rc1-1.0
He3   Os184   endl2009.4-rc1-1.0
He4   Os184   endl2009.4-rc1-1.0

listProtaresInMap.py all.map 'He.*' '.*Os.*4.*'
He3   Os184   endl2009.4-rc1-1.0
He4   Os184   endl2009.4-rc1-1.0
'''

parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('map',                                      help='''The path to a map file.''')
parser.add_argument('projectile',   default='', nargs='?',      help='''The projectile's PoPs ID.''')
parser.add_argument('target',       default='', nargs='?',      help='''The target's PoPs ID.''')
parser.add_argument('evaluation',   default='', nargs='?',      help='''The evalaution label.''')

args = parser.parse_args( )

mapFile = mapModule.Map.readXML_file(args.map)

projectile = args.projectile if args.projectile != '' else None
target = args.target if args.target != '' else None
evaluation = args.evaluation if args.evaluation != '' else None

protares = mapFile.findAllOf(projectile=projectile, target=target, evaluation=evaluation)

projectileMaximumLength = 0
targetMaximumLength = 0
evaluationMaximumLength = 0
for protare in protares:
    projectileMaximumLength = max(projectileMaximumLength, len(protare.projectile))
    targetMaximumLength = max(targetMaximumLength, len(protare.target))
    evaluationMaximumLength = max(evaluationMaximumLength, len(protare.evaluation))

format = '%%-%ds   %%-%ds   %%s' % (projectileMaximumLength, targetMaximumLength)
for protare in protares:
    print(format % (protare.projectile, protare.target, protare.evaluation))
