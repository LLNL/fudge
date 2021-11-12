# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import numpy

from xData.interactivePlot import multiplot as interactivePlotModule

from fudge import reactionSuite as reactionSuiteModule

# assume input GNDS file is in this script's folder
filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'n-094_Pu_239_reactionSuite.gnds.xml')
reactionSuite = reactionSuiteModule.readXML(filename)

fissionReaction = reactionSuite.reactions['m(E)*n + photon [total fission]']
fissionDistributionXYs = fissionReaction.outputChannel.products['n'].distribution.evaluated.energySubform.data
fissionDistributionNested = fissionDistributionXYs.copyDataToNestedLists()
variableLabels = [x.label for x in fissionDistributionXYs.axes]
variableLabels.reverse()

xValues = fissionDistributionXYs.domainGrid
xLabel = '%s (%s)' % (' '.join([x.capitalize() for x in variableLabels[0].split('_')]),
                      fissionDistributionXYs.domainUnit)

potentialYValues = [x.domainGrid for x in fissionDistributionXYs]
boolIdentical2ndGrid = all([set(x).difference(potentialYValues[0]) == set() for x in potentialYValues])
assert boolIdentical2ndGrid, 'Non-equal 2nd domain grid values'
yValues = potentialYValues[0]
yLabel = '%s (%s)' % (' '.join([y.capitalize() for y in variableLabels[1].split('_')]),
                      fissionDistributionXYs[0].domainUnit)

zValues = numpy.array([x.copyDataToXsAndYs()[1] for x in fissionDistributionXYs])
zLabel = '(%s)' % fissionDistributionXYs.axes[0].unit

yValues, xValues = numpy.meshgrid(yValues, xValues)
plotOptions = {'title': 'Total Fission Distribution', 'xLabel': xLabel, 'yLabel': yLabel, 'zLabel': zLabel,
               'xMin': str(xValues.min()), 'xMax': str(xValues.max()), 'yMin': str(yValues.min()),
               'yMax': str(yValues.max()), 'zMin': str(zValues.min()), 'zMax': str(zValues.max())}
plotData = {'totalFission': [xValues, yValues, zValues]}
interactivePlotModule.MultiPlotWithPyQt5(plotOptions, plotData, plotType='3d')
