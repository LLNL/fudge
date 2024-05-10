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
from PoPs import specialNuclearParticleID as specialNuclearParticleIDModule

from fudge.productData import multiplicity as multiplicityModule
from fudge.productData.distributions import unspecified as unspecifiedModule

indentIncrement = '  '

summaryDocStringFUDGE = '''Reads in all protares in a map file and prints all products which have missing/unspecified multiplicity/distribution data.'''

description = '''
This script iterates over all protares in the specified map file. For each protare, it loops over each product (i.e., outgoing particle)
and checks if the product has good multiplicity and distribution data. If a product does not have good data, its xpath is printed.
Good data means that the data exist and it is not unspecified data.
The inputted path can also reference a protare (i.e., GNDS reactionSuite) file in which case just it is checked.

If a list of products is specified after the map file, only those products will be checked. For example, for the reaction 'n + O16 -> n + p + N15',
if n and p are specified after the map file, then the product N15 in this reaction is not checked but n and p are checked.
'''

parser = argparse.ArgumentParser(description=description)

parser.add_argument('mapPath', type=pathlib.Path,                           help='Path to map file whose protares are checked or a protare file.')
parser.add_argument('products', nargs='*',                                  help='The list of products to check. If none present, all products are checked.')
parser.add_argument('--skipProductions', action='store_true',               help='If present, skips the "productions" node.')
parser.add_argument('--skipIncompleteReactions', action='store_true',       help='If present, skips the "incompleteReactions" node.')

productList = []

def checkDelayedNeutrons(self, indent):
    """
    This function checks all the delayed neutron products.
    """

    for delayedNeutron in self:
        checkProduct(delayedNeutron.product, indent)

def checkfissionFragmentData(self, indent):
    """
    This function calls **checkDelayedNeutrons** for self.delayedNeutrons.
    """

    checkDelayedNeutrons(self.delayedNeutrons, indent)

def checkProduct(self, indent):
    """
    This function check a product.
    """

    if len(productList) > 0:
        if self.pid not in productList:
            return

    info = []

    multiplicity = self.multiplicity
    if len(multiplicity) == 0:
        info.append('multiplicity missing')
    else:
        multiplicityFrom = multiplicity[0]
        if isinstance(multiplicityFrom, multiplicityModule.Unspecified):
            info.append('multiplicity unspecified')

    distribution = self.distribution
    if len(distribution) == 0:
        info.append(' distribution missing')
    else:
        distributionForm = distribution[0]
        if isinstance(distributionForm, unspecifiedModule.Form):
            info.append(' distribution unspecified')

    if len(info) > 0:
        print('    %s: %s' % (self.toXLink(), ','.join(info)))

def checkOutputChannel(self, indent):
    """
    This function checks all the products in an outputChannel (i.e., *self*) and its fission fragment data.
    """

    for product in self.products:
        checkProduct(product, indent)

    checkfissionFragmentData(self.fissionFragmentData, indent)

def checkReaction(self, indent):
    """
    This function calls *8checkOutputChannel** on a reactions outputChannel.
    """

    checkOutputChannel(self.outputChannel, indent + indentIncrement)

def showChildReaction(reactions, skip):
    """
    This function calls **checkReaction** on all reaction in the reaction suite *reactions*.
    """

    if not skip:
        for reaction in reactions:
            checkReaction(reaction, indentIncrement)

if __name__ == '__main__':
    args = parser.parse_args()

    productList = []
    for product in args.products:
        ids = specialNuclearParticleIDModule.findID(product)
        if ids is None:
            productList.append(product)
        else:
            for id in ids:
                productList.append(id)

    map = argumentsForScriptsModule.mapFromMapOrProtarePath(args.mapPath)

    for mapProtare in map.iterate( ):
        protare = mapProtare.read(lazyParsing=True)

        print(protare.sourcePath)
        showChildReaction(protare.reactions, False)
        showChildReaction(protare.orphanProducts, False)
        showChildReaction(protare.productions, args.skipProductions)
        showChildReaction(protare.incompleteReactions, args.skipIncompleteReactions)
