#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import shutil
import argparse

from pqu import PQU as PQUModule
from LUPY import argumentsForScripts as argumentsForScriptsModule
from xData import XYs1d as XYs1dModule

from PoPs import IDs as  PoPsIDsModule

from fudge import enums as enumsModule
from fudge.reactionData import crossSection as crossSectionModule
from fudge.reactionData.doubleDifferentialCrossSection.photonScattering import coherent as coherentModule
from fudge.reactionData.doubleDifferentialCrossSection.photonScattering import incoherent as incoherentModule

description = '''
Print cross section data about the requested protare. If input options '--outputDir' and '--plot' are not specified, then
the total cross section is printed to the screen. If the input options '--outputDir', cross section data for each reaction
and total are written to files in the directory specified by '--outputDir' and the protares input channel name.
'''

energyUnitDefault = 'MeV'

parser = argparse.ArgumentParser(description=description)

singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments(parser)

parser.add_argument('--energyUnit', choices=('eV', 'MeV'), default=energyUnitDefault,   help='Energy unit to use. Default is "%s".' % energyUnitDefault)
parser.add_argument('-o', '--outputDir',                                                help='Write cross section for each reaction and total to the specified directory.')
parser.add_argument('--plot', action = 'store_true',                                    help='Draw plot of all reaction cross sections (requires PyQt5).')
parser.add_argument('-v', '--verbose', action='count', default = 0,                     help='Increase verbosity (-v, -vv, -vvv).')

args = parser.parse_args()

protare = singleProtareArguments.protare(args, verbosity=args.verbose, lazyParsing=True)
if args.energyUnit != protare.domainUnit: protare.convertUnits({protare.domainUnit: args.energyUnit})
protareLabel = '%s+%s' % (protare.projectile, protare.target)

outputDir = args.outputDir
if outputDir is not None:
    outputDir = os.path.join(outputDir, protareLabel)

def output(MT, reactionStr, prefix, crossSection, yLabel, subDir=None):

    outputDir2 = outputDir
    if subDir is not None:
        outputDir2 = os.path.join(outputDir2, subDir)
        if not os.path.exists(outputDir2):
            os.makedirs(outputDir2)
    if len(crossSection) == 0: return

    with open(os.path.join(outputDir2, prefix + '_%.3d.dat' % MT), 'w') as fOut:
        fOut.write('# reaction = "%s"\n' % reactionStr)
        fOut.write('# x axes label = "Projectile energy [%s]"\n' % crossSection.axes[1].unit)
        fOut.write('# y axes label = "%s [%s]"\n' % (yLabel, crossSection.axes[0].unit))
        fOut.write(crossSection.toString(format='  %16.8e %16.8e'))

if outputDir is not None:
    if os.path.exists(outputDir): shutil.rmtree(outputDir)
    os.makedirs(outputDir)
    
    outputLog = open(os.path.join(outputDir, 'index'), 'w')
    outputLog.write('index   MT  label                          \n')
    outputLog.write('-------------------------------------------\n')

crossSections = []
if args.verbose > 0: print(protare.sourcePath)
total = crossSectionModule.XYs1d(axes=protare.reactions[0].crossSection[-1].axes)
for reactionCounter, reaction in enumerate(protare.reactions):
    if args.verbose > 1: print('    %s' % reaction)
    crossSection = reaction.crossSection.toPointwise_withLinearXYs(lowerEps=1e-6, upperEps=1e-6)
    if not crossSection.areDomainsMutual(total):
        if args.verbose > 2: print('        Mutualifing domains.')
        total, crossSection = total.mutualify(1e-6, -1e-6, True, crossSection, 1e-6, -1e-6, True)
    crossSection.plotLabel = str(reaction)
    crossSections.append(crossSection)
    total = total + crossSection
    if outputDir is not None:
        outputLog.write('%5d  %3d  %s\n' % (reactionCounter, reaction.ENDF_MT, str(reaction)))
        output(reaction.ENDF_MT, str(reaction), '%3.3d' % reactionCounter, crossSection, 'Cross section')

if outputDir is not None:
    if len(protare.productions) > 0:
        outputLog.write('\nProduction reactions:\n')
    for reactionCounter, reaction in enumerate(protare.productions):
        if args.verbose > 1: print('    Production: %s' % reaction)
        crossSection = reaction.crossSection.toPointwise_withLinearXYs(lowerEps=1e-6, upperEps=1e-6)
        crossSection.plotLabel = str(reaction)
        outputLog.write('%5d  %3d  %s\n' % (reactionCounter, reaction.ENDF_MT, str(reaction)))
        output(reaction.ENDF_MT, str(reaction), '%3.3d' % reactionCounter, crossSection, 'Production cross section', subDir='Productions')

    if protare.interaction == enumsModule.Interaction.atomic and protare.projectile == PoPsIDsModule.photon:
        for reaction in protare.reactions:
            if reaction.ENDF_MT == 502:
                prefix = 'coherentFormFactor'
                formClass = coherentModule.Form
                dataName = 'formFactor'
            elif reaction.ENDF_MT == 504:
                prefix = 'incoherentScatteringFactor'
                formClass = incoherentModule.Form
                dataName = 'scatteringFactor'
            else:
                continue
            for form in reaction.doubleDifferentialCrossSection:
                if isinstance(form, formClass):
                    data = getattr(form, dataName)
                    if hasattr(data, 'data'):
                        data = data.data
                    invLengthToEnergy = (PQUModule.PQU(1, data.axes[1].unit) * PQUModule.PQU('1 hplanck * c')).getValueAs(args.energyUnit)
                    data = data.toPointwise_withLinearXYs(lowerEps=1e-6, upperEps=1e-6)
                    data.setData(data.nf_pointwiseXY.scaleOffsetXAndY(xScale=invLengthToEnergy))
                    data.axes[1].unit = args.energyUnit
                    output(reaction.ENDF_MT, str(reaction), prefix, data, prefix)
                    break

total.plotLabel = 'total'
crossSections.insert(0, total)

if outputDir is None and not args.plot:
    print(total.toString(format='  %16.8e %16.8e'))
elif outputDir is not None:
    output(1, 'total', 'total', total, 'Cross section')
if args.plot: crossSections[0].multiPlot(crossSections, title=protareLabel, xLabel='Energy [%s]' % total.axes.axes[1].unit,
        yLabel='Cross section [%s]' % total.axes[0].unit)
