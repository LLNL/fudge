# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse

from fudge import reactionSuite as reactionSuiteModule


parser = argparse.ArgumentParser(description='Demonstration of the interactive plot method for a single XYs1d instance')
parser.add_argument('filename', metavar='filename', type=str, help='GNDS XML file with a reactionSuite')
parser.add_argument('--reactionLabel', metavar='reactionLabel', type=str, default=None, help='Reaction label')
parser.add_argument('--allReactions', default=False, action='store_true',
    help='Option to plot the evaluated cross sections for either a single or all reactions in the reactionSuite')
args = parser.parse_args()

reactionSuite = reactionSuiteModule.ReactionSuite.readXML_file(args.filename)

if args.allReactions:
    curve1ds = []
    for reaction in reactionSuite.reactions:
        crossSection = reaction.crossSection.toPointwise_withLinearXYs(accuracy=1e-3, lowerEps=1e-8)
        crossSection.plotLegendKey = str(reaction)
        curve1ds.append(crossSection)

    if len(curve1ds)>0:
        curve1ds[0].multiPlot(curve1ds)

else:
    reactionLabel = '%s + %s' % (reactionSuite.projectile, reactionSuite.target) if args.reactionLabel is None \
        else args.reactionLabel
    crossSectionXYs = reactionSuite.reactions[reactionLabel].crossSection.\
        toPointwise_withLinearXYs(accuracy=1e-3, lowerEps=1e-8)
    crossSectionXYs.plot()
