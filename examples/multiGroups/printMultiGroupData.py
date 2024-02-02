#! /usr/bin/env python3
  
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import argparse
import pathlib

from LUPY import argumentsForScripts as argumentsForScriptsModule
from xData import vector as vectorModule
from fudge import reactionSuite as reactionSuiteModule
from fudge.processing import transporting as transportingModule
from fudge.processing.deterministic import multiGroupSettings as multiGroupSettingsModule

summaryDocStringFUDGE_example = '''This script shows how to use GIDI to access multi-group data in a GNDS file.'''

description = '''
This script shows how to use GIDI to access multi-group data in a GNDS file.
'''

doubleFormat = '%14.7e'
file=sys.stdout

parser = argparse.ArgumentParser(description=description)
singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments(parser)
parser.add_argument('outputPath', type=pathlib.Path, default=None, nargs='?',   help='File path where output is written.')
parser.add_argument('--histogram', action='store_true',                         help='If present, all vector data are printed as one (boundary, value) pair per line and will plot as a histogram.')

args = parser.parse_args()

if args.outputPath is not None:
    file = outputPath.open(mode='w')

#
# Read in the protare.
#
protare = singleProtareArguments.protare(args)
settings = multiGroupSettingsModule.MG(protare.projectile, transportingModule.Mode.multiGroup, transportingModule.DelayedNeutrons.off, raiseOnError=True)

#
# This section prints a list of temperatures and their various styles.
#
temperatures = protare.styles.temperatures()
header = '| %-5s | %-12s | %-20s | %-20s | %-20s | %-20s | %-20s |' % ('index', ' temperature', 'multi-group', 'upScatter', 'heated', 'Monte Carlo', 'URR tables')
print('# %s' % (len(header) * '-'))
print('# %s' % header)
print('# %s' % (len(header) * '='))

multiGroups = {}
for index, temperature in enumerate(temperatures):
    print('# | %5d | %12.4e | %-20s | %-20s | %-20s | %-20s | %-20s |' % (index, temperature.temperature, temperature.heatedMultiGroup, 
            temperature.SnElasticUpScatter, temperature.heated, temperature.griddedCrossSection, temperature.URR_probabilityTables))
    if temperature.heatedMultiGroup != '':
        multiGroups[temperature.heatedMultiGroup] = {}
        multiGroupStyle = protare.styles[temperature.heatedMultiGroup]
        for transportable in multiGroupStyle.transportables:
            multiGroups[temperature.heatedMultiGroup][transportable.label] = [boundary for boundary in transportable.group.boundaries.values]
print('# %s' % (len(header) * '-'))
print()

#
# This section prints multi-group boundaries.
#
for label in multiGroups:
    print(label)
    multiGroups2 = multiGroups[label]
    for pid in multiGroups2:                            # Keys are particle GNDS ids.
        print('    %s %s' % (pid, multiGroups2[pid]))
    break

def printVector(prefix, vector, temperature):
    '''
    Prints the vector as a histogram (group boundary, value) pairs if option '--histogram' present. Otherwise, calls printVector2.
    '''

    if not args.histogram:
        printVector2(prefix, vector)
    else:
        print('# %s::' % prefix)
        boundaries = multiGroups[temperature.heatedMultiGroup][protare.projectile]
        for index, value in enumerate(vector):
            print('%s %s' % (doubleFormat % boundaries[index], doubleFormat % value))
            print('%s %s' % (doubleFormat % boundaries[index+1], doubleFormat % value))

def printVector2(prefix, vector):
    '''Prints the prefix and vector data on one line.'''

    print('# %s::' % prefix, end='')
    for value in vector:
        print(doubleFormat % value, end=' ')
    print()

def printMatrix(prefix, matrix):
    '''Prints a matrix by calling printVector2 for each row of the matrix.'''

    prefix = '%s (shape = %s)' % (prefix, matrix.shape)
    rowPrefix = len(prefix) * ' '
    for row in matrix.matrix:
        printVector2(prefix, row)
        prefix = rowPrefix

temperature = temperatures[0]
print()
printVector('Inverse speeds', protare.multiGroupInverseSpeed(temperature), temperature)
printVector('Cross section', protare.multiGroupCrossSection(settings, temperature), temperature)

#
# This section prints the multiplicity for each particle in the for loop.
#
print()
for pid in ['n', 'p', 'H1', 'H2', 'H3', 'He3', 'He4', 'photon']:
    try:
        printVector('Multiplicity for %s' % pid, protare.multiGroupMultiplicity(settings, temperature, pid), temperature)
    except:
        print('# Particle %s has incomplete data.' % pid)

#
# This section prints deposition energy assuming that the particles in the variable particles are being transported.
#
print()
particles = ['n', 'photon']
printVector('Deposition energy', protare.multiGroupDepositionEnergy(settings, temperature, particles), temperature)

depositionEnergySum = vectorModule.Vector()
for reaction in protare.reactions:          # Deposition energy by reaction.
    reactionDepositionEnergy = reaction.multiGroupDepositionEnergy(settings, temperature, particles)
    depositionEnergySum += reactionDepositionEnergy
    printVector('Deposition energy for %s' % reaction.label, reactionDepositionEnergy, temperature)
#
# Need to substract average energy for each tracked particle to get deposition energy sum.
#
for orphanProduct in protare.orphanProducts:
    for particle in particles:
        averageEnergy = orphanProduct.multiGroupAverageEnergy(settings, temperature, particle)
        depositionEnergySum -= averageEnergy
        printVector('Average energy for %s:%s' % (orphanProduct.label, particle), averageEnergy, temperature)
printVector('Deposition energy sum', depositionEnergySum, temperature)

#
# This section prints the Legendre order = 0 transfer matrix for each particle in the for loop.
#
print()
for pid in ['n']:
    printMatrix('Product maxtrix for %s' % pid, protare.multiGroupProductMatrix(settings, temperature, particles, pid, 0))
