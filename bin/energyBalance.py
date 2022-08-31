#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import argparse
import pathlib

from LUPY import argumentsForScripts as argumentsForScriptsModule
from xData import XYs1d as XYs1dModule
from fudge import styles as stylesModule
from fudge.productData import averageProductEnergy as averageProductEnergyModule

description = '''
This script calculates, for each reaction, the available energy, the energy to each product and the energy balance. The results are
written into files.  Information (e.g., index, MT and label) about each reaction is written to the "energyBalance.log" file
in the output directory.  All files for a reaction are written into its own sub-directory under the output directory. The name 
of a reaction's sub-directory is its index.
'''

apdLabel = 'apd'
lowerUpperEpsilon = 1e-6

parser = argparse.ArgumentParser(description=description)

singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments(parser)

parser.add_argument('outputDir', type=pathlib.Path,                         help='Director where all output files are written.')
parser.add_argument('-o', '--output', type=pathlib.Path, default = '',      help='If present, file name to write the protare which contains the added energy data.')
parser.add_argument('--plot', action='store_true',                          help='If present, the energy balance by reaction curves are plotted.')
parser.add_argument('-v', '--verbose', action='count', default=0,           help='If present, prints infomation.')
args = parser.parse_args()

protare = singleProtareArguments.protare(args)
averageProductEnergyAxes = averageProductEnergyModule.defaultAxes(protare.domainUnit)

def addCurves1d(curve1, curve2):

    if not isinstance(curve1, XYs1dModule.XYs1d):
        curve1 = curve1.toPointwise_withLinearXYs(accuracy=1e-5, lowerEps=lowerUpperEpsilon, upperEps=lowerUpperEpsilon)
    if not isinstance(curve2, XYs1dModule.XYs1d):
        curve2 = curve2.toPointwise_withLinearXYs(accuracy=1e-5, lowerEps=lowerUpperEpsilon, upperEps=lowerUpperEpsilon)

    if len(curve1) > 0:
        if curve1.interpolation != curve2.interpolation:
            if curve1.interpolation != xDataEnumsModule.Interpolation.linlin:
                curve1 = curve1.changeInterpolation(xDataEnumsModule.Interpolation.linlin, accuracy, lowerUpperEpsilon, lowerUpperEpsilon)
            if curve2.interpolation != xDataEnumsModule.Interpolation.linlin:
                curve2 = curve2.changeInterpolation(xDataEnumsModule.Interpolation.linlin, accuracy, lowerUpperEpsilon, lowerUpperEpsilon)
    if not curve1.areDomainsMutual(curve2):
        curve1, curve2 = curve1.mutualify(lowerUpperEpsilon, lowerUpperEpsilon, True, curve2, lowerUpperEpsilon, lowerUpperEpsilon, True)

    return curve1 + curve2

def output(reactionOutputDir, fileName, curve):

    outputName = reactionOutputDir / fileName
    with outputName.open('w') as fOut:
        print(curve.toString(), file=fOut, end='')

def outputProductData(reactionOutputDir, outputChannel, productSums, level):

    for productIndex, product in enumerate(outputChannel):
        if args.verbose > 1:
            print((level + 1) * '    ', '    %s' % product)
        if product.outputChannel is None:
            if apdLabel in product.averageProductEnergy:
                averageProductEnergy = product.averageProductEnergy[apdLabel]
            else:
                averageProductEnergy = XYs1dModule.XYs1d(axes=averageProductEnergyAxes)
            fileName = 'averageProductEnergy_%.2d_%.2d_%s.dat' % (level, productIndex, product.pid)
            output(reactionOutputDir, fileName, averageProductEnergy)
            if product.pid not in productSums:
                productSums[product.pid] = XYs1dModule.XYs1d(axes=averageProductEnergyAxes)
            productSums[product.pid] = addCurves1d( productSums[product.pid], averageProductEnergy )
        else:
            outputProductData(reactionOutputDir, product.outputChannel, productSums, level + 1)

rootStyle = protare.styles.preProcessingChainHead()

AEPStyle = protare.styles.getStyleOfClass(stylesModule.AverageProductData)
if AEPStyle is None:
    AEPStyle = stylesModule.AverageProductData(apdLabel, rootStyle.label)
    protare.styles.add(AEPStyle)
    protare.calculateAverageProductData(AEPStyle, indent='  ')

if args.outputDir.exists():
    os.system('rm -rf %s' % args.outputDir)
args.outputDir.mkdir()

fLog = open(args.outputDir / 'energyBalance.log', 'w')
print('# file name: %s' % protare.sourcePath, file=fLog)
print('# absolute path: %s' % pathlib.Path(protare.sourcePath).resolve(), file=fLog)
print('# energy unit = %s' % protare.domainUnit, file=fLog)
print(file=fLog)

print(' index |  MT | reaction label', file=fLog)
print( 50 * '=', file=fLog)

reactionIndex = 0
energyBalanceCurves = []

def process(reaction, isReaction):

    global reactionIndex

    if args.verbose > 0:
        print('%5d: %s' % (reactionIndex, reaction))
    print(' %5d | %3d | %s' % (reactionIndex, reaction.ENDF_MT, reaction), file=fLog)

    reactionOutputDir = args.outputDir / ('%.4d_%.3d' % (reactionIndex, reaction.ENDF_MT))
    reactionOutputDir.mkdir()

    if isReaction:
        availableEnergy = reaction.availableEnergy[apdLabel]
        output(reactionOutputDir, 'availableEnergy.dat', availableEnergy)

    productSums = {}
    outputProductData(reactionOutputDir, reaction.outputChannel, productSums, 0)

    totalProductEnergy = XYs1dModule.XYs1d(axes=averageProductEnergyAxes)
    for pid in productSums:
        totalProductEnergy = addCurves1d(totalProductEnergy, productSums[pid])
        output(reactionOutputDir, '%s_averageProductEnergy.dat' % pid, productSums[pid])
    output(reactionOutputDir, 'totalProductEnergy.dat', totalProductEnergy)

    if isReaction:
        energyBalance = addCurves1d(availableEnergy, -totalProductEnergy)
        output(reactionOutputDir, 'energyBalance.dat', energyBalance)
    else:
        energyBalance = totalProductEnergy
    energyBalance.plotLabel = str(reaction)
    energyBalanceCurves.append(energyBalance)

    reactionIndex += 1

for reaction in protare.reactions:
    process(reaction, True)
for reaction in protare.orphanProducts:
    process(reaction, False)

if args.plot:
    XYs1dModule.XYs1d.multiPlot(energyBalanceCurves, title='Energy balance per reaction')

if str(args.output) != '.':
    protare.saveToFile(args.output)

fLog.close()
