#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import argparse

from pqu import PQU as PQUModule
from fudge import reactionSuite as reactionSuiteModule
from fudge.reactionData import crossSection as crossSectionModule

summaryDocString__FUDGE = '''Prints summary information about resonances for a GNDS reactionSuite file.'''

description = '''
Prints summary information about the resolved and unresolved resonance regions, and each reaction that has resonance data.
'''

floatToShortestString = PQUModule.floatToShortestString

parser = argparse.ArgumentParser(description=description)
parser.add_argument('file', type=str,                                       help='Input reactionSuite file.')
parser.add_argument('--energyUnit', type=str, default=None,                 help='Convert all energies in the gnds file to this unit.')

args = parser.parse_args()

protare = reactionSuiteModule.read(args.file)
if args.energyUnit is not None and args.energyUnit != protare.domainUnit:
    protare.convertUnits({protare.domainUnit: args.energyUnit})

def function1dInfo(label, resonanceRegion):

    if resonanceRegion is None:
        return
    function1d = resonanceRegion.data
    print('            %-10s: %-10s (%s to %s)' % \
            (label, function1d.moniker, floatToShortestString(function1d.domainMin), floatToShortestString(function1d.domainMax)))
    if isinstance(function1d, crossSectionModule.Regions1d):
        for region in function1d:
            print('                %s (%s to %s)' % (region.moniker, floatToShortestString(region.domainMin), floatToShortestString(region.domainMax)))

def reactionInfo(self):

    crossSection = self.crossSection
    printReactionLabel = True
    for form in crossSection:
        if isinstance(form, crossSectionModule.ResonancesWithBackground):
            if printReactionLabel:
                print('        %s (%s)' % (self, self.moniker))
            background = form.background
            function1dInfo('resolved',   background.resolvedRegion)
            function1dInfo('unresolved', background.unresolvedRegion)
            function1dInfo('fast',       background.fastRegion)
            return True

    return False

resolved = None
unresolved = None

if protare.resonances is not None:
    resolved = protare.resonances.resolved
    unresolved = protare.resonances.unresolved
if resolved is None and unresolved is None:
    sys.exit()

print('    Energy unit is "%s".' % protare.domainUnit)

resonanceDomainMax = -1
if resolved is not None:
    resonanceDomainMax = resolved.domainMax
    print('    Resolved domain: %s %s' % (floatToShortestString(resolved.domainMin), floatToShortestString(resolved.domainMax)))

if unresolved is not None:
    print('    Unresolved domain: %s %s' % (floatToShortestString(unresolved.domainMin), floatToShortestString(unresolved.domainMax)))

reactionsInResonanceRegion = []
print('    Cross sections:')
for reaction in protare.reactions:
    if not reactionInfo(reaction):
        if reaction.effectiveThreshold() < resonanceDomainMax:
            reactionsInResonanceRegion.append(reaction)
for crossSectionSum in protare.sums.crossSectionSums:
    reactionInfo(crossSectionSum)

if len(reactionsInResonanceRegion) > 0:
    print('    Other reaction in resonance domain:')
    for reaction in reactionsInResonanceRegion:
        print('        %s' % reaction)
