#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from argparse import ArgumentParser

from LUPY import argumentsForScripts as argumentsForScriptsModule
from fudge.reactions import base as reactionBaseModule
from fudge import outputChannel as outputChannelModule
from fudge import product as productModule

indentIncrement = '  '

summaryDocStringFUDGE = '''Prints disc usage of various nodes in a GNDS reactionSuite file.'''

description = '''
Prints disc usage of various nodes in a GNDS file. Currently, only looks at disc usage for a file written to XML.
If no '-v' options are present, empty nodes are not printed.
'''

parser = ArgumentParser(description=description)
singleProtareArguments = argumentsForScriptsModule.SingleProtareArguments(parser)
parser.add_argument('-b', '--brief', action='store_true',                   help='If present, limits depth of calls.')
parser.add_argument('-v', '--verbose', action='count', default=0,           help='The more, the more verbose.')

def printDiscSpaceUsage(indent, name, self):

    discSpaceUsage = len(self.toXML())
    if discSpaceUsage > 0 or args.verbose > 0:
        header = '%s %s' % (indent, name)
        print('%-64s: %12d' % (header, discSpaceUsage))

def productDiscSpaceUsage(self, indent):

    printDiscSpaceUsage(indent, self, self)

productModule.Product.__discSpaceUsage = productDiscSpaceUsage

def outputChannelDiscSpaceUsage(self, indent):

    indent2 = indent + indentIncrement

    printDiscSpaceUsage(indent, self, self)
    for product in self.products:
        product.__discSpaceUsage(indent2)

    printDiscSpaceUsage(indent2, 'fissionFragmentData', self.fissionFragmentData)

outputChannelModule.OutputChannel.__discSpaceUsage = outputChannelDiscSpaceUsage

def reactionDiscSpaceUsage(self, indent):

    indent2 = indent + indentIncrement

    printDiscSpaceUsage(indent, self, self)
    if args.brief:
        return
    printDiscSpaceUsage(indent2, 'double differential cross section', self.doubleDifferentialCrossSection)
    printDiscSpaceUsage(indent2, 'cross section', self.crossSection)
    self.outputChannel.__discSpaceUsage(indent2)
    printDiscSpaceUsage(indent2, 'available energy', self.availableEnergy)
    printDiscSpaceUsage(indent2, 'available momentum', self.availableMomentum)

reactionBaseModule.Base_reaction.__discSpaceUsage = reactionDiscSpaceUsage

def showReaction(reactions, indent):

    indent2 = indent + indentIncrement

    printDiscSpaceUsage(indent, reactions.moniker, reactions)
    for reaction in reactions:
        reaction.__discSpaceUsage(indent2)

if __name__ == '__main__':
    args = parser.parse_args()
    protare = singleProtareArguments.protare(args, verbosity=0, lazyParsing=False)

    indent = indentIncrement
    indent2 = indent + indentIncrement

    printDiscSpaceUsage(indent, protare, protare)
    showReaction(protare.reactions, indent2)
    showReaction(protare.orphanProducts, indent2)
    showReaction(protare.productions, indent2)
    showReaction(protare.incompleteReactions, indent2)
