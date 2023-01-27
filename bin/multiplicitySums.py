#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import argparse
import pathlib

from LUPY import argumentsForScripts as argumentsForScriptsModule
from fudge import sums as sumsModule
from PoPs import IDs as PoPsIDsModule

description = '''Creates multiplicitySum nodes for delayed and total fission neutron multiplicity if needed.'''

parser = argparse.ArgumentParser(description=description)
singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments(parser)
parser.add_argument('--output', nargs='?', default=None,                                    help='Path for the output file.')

args = parser.parse_args()

protare = singleProtareArguments.protare(args)

fissionReaction = protare.getReaction('fission')
if fissionReaction is not None:
    MTs = {}
    for multiplicitySum in protare.sums.multiplicitySums:
        MTs[multiplicitySum.ENDF_MT] = multiplicitySum
    if 455 in MTs:
        MultiplicitySum = MTs[455]
    else:
        delayedNeutronMultiplicity = None
        MultiplicitySum = sumsModule.MultiplicitySum('delayed fission neutron multiplicity', 455)
        for delayedNeutron in fissionReaction.outputChannel.fissionFragmentData.delayedNeutrons:
            MultiplicitySum.summands.summands.append(sumsModule.Add(delayedNeutron.product.multiplicity))
            if delayedNeutronMultiplicity is None:
                delayedNeutronMultiplicity = delayedNeutron.product.multiplicity[0].copy()
                delayedNeutronMultiplicity.label = delayedNeutron.product.multiplicity[0].label
            else:
                delayedNeutronMultiplicity += delayedNeutron.product.multiplicity[0]
        if delayedNeutronMultiplicity is not None:
            MultiplicitySum.multiplicity.add(delayedNeutronMultiplicity)
            protare.sums.multiplicitySums.add(MultiplicitySum)

    if 452 not in MTs:
        neutrons = [product for product in fissionReaction.outputChannel.products if product.pid == PoPsIDsModule.neutron]
        if len(neutrons) != 1:
            raise Exception('Number of prompt neutrons not 1 but %s' % len(neutrons))

        MultiplicitySum2 = sumsModule.MultiplicitySum('total fission neutron multiplicity', 452)
        MultiplicitySum2.summands.summands.append(sumsModule.Add(neutrons[0].multiplicity))
        MultiplicitySum2.summands.summands.append(sumsModule.Add(MultiplicitySum))

        multiplicity = neutrons[0].multiplicity[0] + MultiplicitySum.multiplicity[0]
        multiplicity.label = neutrons[0].multiplicity[0].label
        MultiplicitySum2.multiplicity.add(multiplicity)
        protare.sums.multiplicitySums.add(MultiplicitySum2)

output = args.output
if output is None:
    output = pathlib.Path(protare.sourcePath).with_suffix('.msums.xml')

protare.saveToFile(output)
