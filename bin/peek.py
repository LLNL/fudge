#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from argparse import ArgumentParser

from LUPY import argumentsForScripts as argumentsForScriptsModule

from PoPs import IDs as PoPsIDsModule

from fudge.reactions import base as reactionBaseModule
from fudge import outputChannel as outputChannelModule
from fudge import product as productModule
from fudge import sums as sumsModule
from fudge.reactionData import crossSection as crossSectionModule
from fudge.productData import multiplicity as multiplicityModule
from fudge.productData.distributions import unspecified as unspecifiedModule
from fudge.productData.distributions import uncorrelated as uncorrelatedModule
from fudge.productData.distributions import energy as energyModule
from fudge.outputChannelData.fissionFragmentData import fissionFragmentData as fissionFragmentDataModule
from fudge.outputChannelData.fissionFragmentData import delayedNeutron as delayedNeutronModule

indentIncrement = '  '

description = '''
Prints each reaction and brief information about each reaction's products. Product information include its id, label, and distribution type and frame.
If option '--summaryGammas' is present, only a summary of each photon (gamma) for each output channel is printed.
'''


parser = ArgumentParser(description = description)

singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments(parser)

parser.add_argument('--summaryGammas', action = 'store_true',              help = 'If present, photons (e.g., gammas) for each output channel are summerized.')
parser.add_argument('--productPath', action = 'store_true',                help = 'If present, the product path needed by other scripts (e.g., spectrum.py) are included with each product.')
parser.add_argument('--unspecified', action = 'store_true',                help = 'If present, only information for products with unspecified multiplicities or distributions are printed.')
parser.add_argument('--doNotShowProducts', action = 'store_true',          help = 'If present, no product data are displayed.')

args = parser.parse_args()

protare = singleProtareArguments.protare(args, verbosity=0, lazyParsing=True)

class PhotonDummy:

    pid = PoPsIDsModule.photon
    label = 'many'

def PID_label(indent, product):

    return '%-50s' % ('%sid = %-10s label = %-16s' % (indent, product.pid, product.label))

def delayedNeutronsPeek(self, indent):

    if len(self) > 0: print('%s-- delayedNeutrons --' % indent)
    indent2 = indentIncrement + indent
    for delayedNeutron in self:
        delayedNeutron.product.__peek(indent2, [], delayedNeutronRate = 'label = %-6s rate = %s' % (delayedNeutron.label, delayedNeutron.rate[0]))

delayedNeutronModule.DelayedNeutrons.__peek = delayedNeutronsPeek

def fissionFragmentDataPeek(self, indent):

    self.delayedNeutrons.__peek(indent)

fissionFragmentDataModule.FissionFragmentData.__peek = fissionFragmentDataPeek

def productPeek(self, indent, productPath, doPrint = True, delayedNeutronRate = None):

    productPath = productPath + [ self.label ]

    unspecifiedPresent = False
    unspecifiedMultiplicity = ''
    distributionInfo = 'unknown'

    multiplicity = self.multiplicity
    if len(multiplicity) == 0:
        print('%s missing multiplicity data present' % PID_label(indent, self))
        return True, distributionInfo
    multiplicityFrom = multiplicity[0]
    if isinstance(multiplicityFrom, multiplicityModule.Unspecified):
        unspecifiedPresent = True
        unspecifiedMultiplicity = ' (multiplicity-unspecified)'

    distribution = self.distribution
    if len(distribution) == 0:
        print('%s missing distribution data present' % PID_label(indent, self))
        return True, distributionInfo
    distributionForm = distribution[0]
    if isinstance(distributionForm, unspecifiedModule.Form): unspecifiedPresent = True
    distributionInfo = '%s (%s)' % (distribution[0].moniker, distribution[0].productFrame)
    if isinstance(distributionForm, uncorrelatedModule.Form) and not args.summaryGammas:
        energy = distributionForm.energySubform.data
        if isinstance(energy, energyModule.DiscreteGamma):
            distributionInfo += ': discrete gamma %g %s' % (energy.value, energy.axes[1].unit)
        elif isinstance(energy, energyModule.PrimaryGamma):
            distributionInfo += ': primary gamma %g %s' % (energy.value, energy.axes[1].unit)

    if args.unspecified:
        if len(unspecifiedMultiplicity) == 0 and distributionForm.moniker != unspecifiedModule.Form.moniker:
            if self.outputChannel is not None: self.outputChannel.__peek(indent + indentIncrement, productPath)
            return unspecifiedPresent, distributionInfo

    productPathStr = ''
    if delayedNeutronRate is None:
        if args.productPath: productPathStr = ' (' + ':'.join(productPath) + ')'
        if doPrint: print('%s distribution[0] = %s%s%s' % (PID_label(indent, self), distributionInfo, unspecifiedMultiplicity, productPathStr))
        if self.outputChannel is not None: self.outputChannel.__peek(indent + indentIncrement, productPath)
    else:
        print('%-50s distribution[0] = %s%s' % (indent + delayedNeutronRate, distributionInfo, unspecifiedMultiplicity))

    return unspecifiedPresent, distributionInfo

productModule.Product.__peek = productPeek

def outputChannelPeek(self, indent, productPath):

    photonsToSummarize = []
    if args.summaryGammas:
        photonsToSummarize = [ product.pid for product in self.products if product.pid == PoPsIDsModule.photon ]
        if len(photonsToSummarize) == 1: photonsToSummarize = []

    photonsInfo = {}
    photonMonikers = []
    unspecifiedPresent = False
    for product in self.products:
        if product.pid in photonsToSummarize:
            unspecifiedPresent2, distributionInfo = product.__peek(indent, productPath, False)
            unspecifiedPresent = unspecifiedPresent or unspecifiedPresent2
            if distributionInfo not in photonsInfo: photonsInfo[distributionInfo] = 0
            photonsInfo[distributionInfo] += 1
        else:
            product.__peek(indent, productPath)

    if len(photonsInfo) > 0:
        if not args.unspecified or args.unspecified and unspecifiedPresent:
            photonsInfoStr = ', '.join([ '%d %s' % (photonsInfo[form], form) for form in photonsInfo ])
            print('%s distribution[0] = %s' % (PID_label(indent + indentIncrement, PhotonDummy), photonsInfoStr))

    self.fissionFragmentData.__peek(indent)

outputChannelModule.OutputChannel.__peek = outputChannelPeek

def reactionPeek(self, prefix, index, indent):

    crossSectionStr = ''
    crossSection = self.crossSection[0]
    if isinstance(crossSection, crossSectionModule.Reference):
        try:
            crossSectionSum = crossSection.link.findClassInAncestry(sumsModule.CrossSectionSum)
            crossSectionStr = ': %s' % crossSectionSum.label
        except:
            pass
    print('%s%-32s (%4d): domainMin = %s, domainMax = %s %s%s' % (indent, prefix % str(self), index, self.domainMin, self.domainMax, self.domainUnit, crossSectionStr))
    productPath = [str(index)]
    if not args.doNotShowProducts: self.outputChannel.__peek(indent + indentIncrement, [str(reactionIndex)])

reactionBaseModule.Base_reaction.__peek = reactionPeek

indent = indentIncrement
print(protare.sourcePath, '(GNDS-%s)' % protare.format)
print('%sreactions:' % indent)
for reactionIndex, reaction in enumerate(protare.reactions): reaction.__peek('%s', reactionIndex, 2 * indentIncrement)
for reactionIndex, reaction in enumerate(protare.orphanProducts): reaction.__peek('orphanProduct %s', reactionIndex, 2 * indentIncrement)
