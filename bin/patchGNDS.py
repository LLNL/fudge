#! /usr/bin/env python3

# What is the "source" options?

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

# Still to look for links in distributions/photonScattering.py and distributions/reference.py and resonances.

description = """
Applies a patch file to a GNDS file. First the two GNDS files specified in the patch file are read in. 
A new "evaluated" style is create in the first GNDS file with the label specified by the "label" argument 
and a version specified by the "version" argument.  If the option "--library" is specified, the new "evaluated" 
style's library attribute will be set to it; otherwise, the library from the prior "evaluated" style will be used.
Next, the list of diffs specfied in the patch file are scanned and applied to the first file (e.g., data from 
the second GNDS file may be copied to the first file).  All new data inserted into the patched GNDS file will 
have the label specified by the "label" argument.

If the "--forENDF" option is present, for any reaction that products (outgoing particles) are replaced, all other
prouducts for the reaction are check to replaced if need to insure that the data can be written to an ENDF file.
Currently, this means that the products' frame fit one of the ENDF MF=6 LCT values.
"""

import sys
import os
import argparse
import textwrap

from PoPs import IDs as PoPsIdsModule
from PoPs.chemicalElements import misc as PoPsChemicalElementsMiscModule
from xData import date as dateModule
from xData.Documentation import computerCode as computerCodeModule

import fudge
from fudge import enums as enumsModule
from fudge import styles as stylesModule
from fudge import reactionSuite as reactionSuiteModule
from fudge import outputChannel as outputChannelModule
from fudge import product as productModule
from fudge.reactions import base as reactionBaseModule
from fudge.productData.distributions import angular as angularModule
from fudge.productData.distributions import unspecified as unspecifiedModule

reconstructResonancesLabelDefault = 'recon'

parser = argparse.ArgumentParser( description = description )
parser.add_argument( "label", type = str,                                               help = """The label for the new "evaluated" style.""" )
parser.add_argument( "version", type = str,                                             help = """The library version for the new "evaluated" style.""" )
parser.add_argument( "patchFile",                                                       help = """The name of the path file.""" )
parser.add_argument( "outputFile",                                                      help = """The name of the file to write the patched GNDS data to.""" )
parser.add_argument( "--writeAlways", action = 'store_true',                            help = """If there are no changes, the default action is to not write the output file, with this option, the output file is always written.""" )
parser.add_argument( "--energyUnit", action = 'store', default = None,                  help = """Energy units for the new patched file.""" )
parser.add_argument( "--library", type = str, action = 'store', default = "",           help = """The library for the new "evaluated" style. If not present, use the latest in the file to patch.""" )
parser.add_argument( "--forENDF", action = 'store_true',                                help = """If present, attempts to make the final file ENDF compatible.""" )
parser.add_argument( "--reconstructResonances", action = 'store_true', default = False, help = """If present resonances are reconstructed.""" )
parser.add_argument( "--reconstructResonancesLabel", action = 'store', default = reconstructResonancesLabelDefault,
                                                                                        help = '''Label to use if reconstructing resonances. Default = "%s".''' % reconstructResonancesLabelDefault )
parser.add_argument( '--documentation', type = str, default = '',                       help = '''A simple text put into the new evaluated style's documentation body.''' )
parser.add_argument( '--documentationFile', default = None,                             help = '''A file whose content will be put into the new evaluated style's documentation body. If present, the option '--documentation' will be ignored.''' )
parser.add_argument( '--endfCompatible', type = str, default = '',                      help = '''A simple text put into the new evaluated style's endfCompatible documentation.''' )
parser.add_argument( '--endfCompatibleFile', default = None,                            help = '''A file whose content will be put into the new evaluated style's endfCompatible documentation. If present, the option '--endfCompatible' will be ignored.''' )
parser.add_argument('--verbose', '-v', action = 'count', default = 0,                   help = '''Amount of verbosity''')

args = parser.parse_args( )

if not os.path.exists(args.patchFile):
    print('WARNING: patch file "%s" does not exist.' % args.patchFile)
    sys.exit(1)

fIn = open( args.patchFile )
lines = fIn.readlines( )
fIn.close( )

if( ( len( lines ) - 2 ) % 3 != 0 ) : raise Exception( "Input file does not have the correct number of lines to be a valid patch file." )

class ChangeLog:

    def __init__(self):

        self.lines = []

    def addText(self, text):

        if args.verbose > 1: print(text)
        self.lines.append(text)

    def addLines(self, lines):

        for line in lines: self.addText(line)

    def toString(self, endfCompatible):

        lines = []
        for line in self.lines:
            if endfCompatible: line = textwrap.fill(line, width = 66, subsequent_indent = '    ')
            lines.append(line)

        return '\n'.join(lines)

def getReactionLabel(instance):

    return instance.findClassInAncestry(reactionBaseModule.Base_reaction).label

changeLog = ChangeLog()
changeLog.addLines([ '', 'Changes by program %s:' % os.path.basename(__file__)])

def checkFileLine( index ) :

    line = lines[index]
    if( ":" not in line ) : raise Exception( "Invalid file line at line %d is not valid." % ( index + 1 ) )
    tag, fileName = line.split( ":" )
    if( "FILE%s" % ( index + 1 ) != tag ) : raise Exception( "Invalid file line at line %d is not valid." % ( index + 1 ) )

    return( fileName.strip( ) )

def replaceProduct(protare1, protare2, fromComponent, toComponent, forENDF):

    reaction = toComponent.findClassInAncestry(reactionBaseModule.Base_reaction)
    product = toComponent.findClassInAncestry(productModule.Product)
    changeLog.addText('  %s distribution for MT=%s: %s.' % ( product.pid, reaction.ENDF_MT, reaction.label ))
    if forENDF: changeLog.addText('      per --forENDF option as needed by ENDF LCT flag.')
    keys = fromComponent.keys()
    for key in keys: toComponent.remove(key)
    toComponent.add( fromComponent.evaluated )
    toComponent[fromComponent.evaluated.label].label = args.label

protare1 = reactionSuiteModule.ReactionSuite.readXML_file( checkFileLine( 0 ) )
productsInitially = protare1.listOfProducts()
protare2 = reactionSuiteModule.ReactionSuite.readXML_file( checkFileLine( 1 ) )

if args.energyUnit == None: args.energyUnit = protare1.domainUnit
if protare1.domainUnit != args.energyUnit: protare1.convertUnits({ protare1.domainUnit : args.energyUnit })
if protare2.domainUnit != args.energyUnit: protare2.convertUnits({ protare2.domainUnit : args.energyUnit })

chains = protare1.styles.chains( True )
if( len( chains ) != 1 ) : raise Exception( "Only one styles chain is supported." )

chain = chains[0]
for style in chain :
    if( isinstance( style, stylesModule.Evaluated ) ) : break
if( not( isinstance( style, stylesModule.Evaluated ) ) ) : raise Exception( "At least one Evaluated style must be present in file to patch." )

if args.verbose > 0: print('Patch file = "%s".' % args.patchFile)

if( len( lines ) > 2 ) :

    args.writeAlways = True
    if( args.label in protare1.styles ) : raise ValueError( "Label already in file to patch." )
    if( args.library == "" ) : args.library = style.library
    date = dateModule.Date(resolution=dateModule.Resolution.time)
    evaluateStyle = stylesModule.Evaluated(args.label, style.label, None, None, args.library, args.version, date = date)
    if( style.library == evaluateStyle.library ) :
        if( style.version == evaluateStyle.version ) : raise ValueError( "Library and version must be different than one in patch file." )
    protare1.styles.add( evaluateStyle )

    if( args.documentationFile is None ):
        documentation = args.documentation
    else :
        with open( args.documentationFile ) as file:
            documentation = ''.join( file.readlines( ) )

    if( args.endfCompatibleFile is None ):
        endfCompatible = args.endfCompatible
    else :
        with open( args.endfCompatibleFile ) as file:
            endfCompatible = ''.join( file.readlines( ) )
    if 'CURRENT_DATE' in endfCompatible: endfCompatible = endfCompatible.replace('$CURRENT_DATE', str(date))
    if 'CURRENT_DATE' in documentation: documentation = documentation.replace('$CURRENT_DATE', str(date))

    links = protare1.findLinks( )

    sumFixes = []
    outputChannelsToCheck = set()
    for index in range( 2, len( lines ), 3 ) :
        line = lines[index]
        if( ":" not in line ) : raise Exception( "Input tag at line %d is not valid." % ( 3 * index ) )
        tag = line.split( ":" )[0]

        toXLink = lines[index+1].strip( )
        fromXLink = lines[index+2].strip( )

        if( tag == "Distribution unspecified - 1" ) :
            fromComponent = protare2.followXPath( fromXLink )
            toComponent = protare1.followXPath( toXLink )
            replaceProduct(protare1, protare2, fromComponent, toComponent, False)
            toOutputChannel = toComponent.findClassInAncestry(outputChannelModule.OutputChannel)
            if args.forENDF:            # Check that other outgoing products with a non unspecified distribution have the same frame as the one replaced.
                toOutputChannel = toComponent.findClassInAncestry(outputChannelModule.OutputChannel)
                fromOutputChannel = fromComponent.findClassInAncestry(outputChannelModule.OutputChannel)
                for toProduct in toOutputChannel.products:   
                    distribution1 = toProduct.distribution[0]
                    if not isinstance(distribution1, unspecifiedModule.Form):
                        if toComponent[0].productFrame != distribution1.productFrame:
                            if toProduct.pid == PoPsIdsModule.photon:
                                distribution1.productFrame = toComponent[0].productFrame
                            else:
                                fromProducts = [ product2 for product2 in fromOutputChannel.products if product2.pid == toProduct.pid ]
                                if len(fromProducts) == 1:
                                    print('INFO: per "--forENDF" option, replacing distribution for product:\n    %s' % toProduct.toXLink())
                                    replaceProduct(protare1, protare2, fromProducts[0].distribution, toProduct.distribution, True)
                                else:
                                    print('WARNING: could not figure out reguired "--forENDF" option for outputChannel. Found %s want 1 in:\n    %s.' %
                                            ( len(fromProducts), toComponent.toXLink() ))
            outputChannelsToCheck.add(toOutputChannel)
        elif( tag == "reactions missing - 1" ) :
            fromReaction = protare2.followXPath( fromXLink )
            fromReaction.amendForPatch( style.label, evaluateStyle.label )
            protare1.reactions.add( fromReaction )
            changeLog.addText('  Adding reaction "%s"' % ( fromReaction.label ))
        elif( tag == "reactions missing - 2" ) :
            sumFixes.append( [ 'removed', toXLink ] )
            toReaction = protare1.followXPath( toXLink )
            changeLog.addText('  Removing reaction "%s"' % ( toReaction.label ))
            protare1.reactions.pop( toReaction.label )
        elif( 'non-constant Q - 2' in tag ) :
            toComponent = protare1.followXPath( toXLink ).ancestor
            fromComponent = protare2.followXPath( fromXLink )
            toComponent.replace( fromComponent )
            fromComponent.label = evaluateStyle.label
            changeLog.addText('  Replacing constant Q with non-constant in reaction')
            changeLog.addText('    %s' % getReactionLabel(toComponent))
        elif 'More than one product with same id and unspecified - 1' in tag:
            pid = line.split(":")[1].split('=')[1].strip().strip('"')
            toChannel = protare1.followXPath(toXLink)
            fromChannel = protare2.followXPath(fromXLink)
            productsToRemove = [product for product in toChannel.products if product.pid == pid]
            for product in productsToRemove:
                toChannel.products.pop(product.label)

            for product in fromChannel.products:
                if product.pid == pid:
                    toChannel.products.add(product)
                    product.multiplicity[0].label = evaluateStyle.label
                    product.distribution[0].label = evaluateStyle.label
            reaction = protare1.getReaction(getReactionLabel(toChannel))
            changeLog.addText('  %s distribution for MT=%s: %s.' % ( product.pid, reaction.ENDF_MT, reaction.label ))
#           changeLog.addText('  Replacing product "%s" in reaction' % pid)
#           changeLog.addText('    %s' % getReactionLabel(toChannel))
        elif( 'More than one product with same id' in tag ) :
            pid = line.split( ":" )[1].split( '=' )[1].strip( ).strip( '"' )
            toChannel = protare1.followXPath( toXLink )
            fromChannel = protare2.followXPath( fromXLink )
            productsToRemove = [ product for product in toChannel.products if product.pid == pid ]
            for product in productsToRemove : toChannel.products.pop( product.label )

            for product in fromChannel.products :
                if( product.pid == pid ) :
                    toChannel.products.add( product )
                    product.multiplicity[0].label = evaluateStyle.label
                    product.distribution[0].label = evaluateStyle.label
            changeLog.addText('  Replacing product "%s" in reaction' % pid)
            changeLog.addText('    %s' % getReactionLabel(toChannel))
        elif( 'Product missing from channel - 2' == tag ) :
            pid = line.split(":")[1].split('=')[1].strip().strip('"')
            toChannel = protare1.followXPath(toXLink)
            productsToRemove = [product for product in toChannel.products if product.pid == pid]
            for product in productsToRemove : toChannel.products.pop(product.label)
        elif( 'Product missing from channel - 1' == tag ) :
            pid = line.split(":")[1].split('=')[1].strip().strip('"')
            toChannel = protare1.followXPath(toXLink)
            fromChannel = protare2.followXPath(fromXLink)
            for product in fromChannel.products:
                if product.pid == pid: toChannel.products.add(product)
        else :
            print("""WARNING: unsupported diff tag "%s".""" % tag)

    crossSectionReconstructedStyles = protare1.styles.getStylesOfClass( stylesModule.CrossSectionReconstructed )
    if( len( crossSectionReconstructedStyles ) > 0 ) :
        if( len( crossSectionReconstructedStyles ) > 1 ) : raise Exception( "Multiple cross section reconstructed styles present which is not supported." )
        crossSectionReconstructedStyle = crossSectionReconstructedStyles[0]
        protare1.removeStyle( crossSectionReconstructedStyle.label )

    if( args.reconstructResonances ) :
        crossSectionReconstructedStyle = stylesModule.CrossSectionReconstructed( args.reconstructResonancesLabel, args.label )
        protare1.reconstructResonances( crossSectionReconstructedStyle, 1e-3 )

    changeLog.addText('')
    changeLog.addText('')
    if args.forENDF:
        endfCompatible += changeLog.toString(True)
    else:
        documentation += changeLog.toString(False)
    evaluateStyle.documentation.body.body = documentation
    evaluateStyle.documentation.endfCompatible.body = endfCompatible

    computerCode = computerCodeModule.ComputerCode('code', os.path.basename(__file__), fudge.__version__)
    computerCode.executionArguments.body = ' '.join(sys.argv[1:])

    inputDeck = computerCodeModule.InputDeck('input', args.patchFile)
    inputDeck.body = '\n' + ''.join(lines)
    computerCode.inputDecks.add(inputDeck)

    evaluateStyle.documentation.computerCodes.add(computerCode)

    protare1.fixLinks( )

    for outputChannel in outputChannelsToCheck:
        twoBodyCounter = 0
        for product in outputChannel.products:
            if isinstance(product.distribution[0], angularModule.TwoBody): twoBodyCounter += 1
        if twoBodyCounter > 0:
            if twoBodyCounter == 2 and len(outputChannel) == 2:
               outputChannel.genre = enumsModule.Genre.twoBody
            else:
                print('ERROR: outputChannel genre cannot be determined: %s' % outputChannel.toXLink())
        else:
            outputChannel.genre = enumsModule.Genre.NBody

    crossSectionSums = protare1.sums.crossSectionSums
    for label in [ crossSectionSum.label for crossSectionSum in crossSectionSums if len(crossSectionSum.summands) == 0]: crossSectionSums.pop(label)

    multiplicitySums = protare1.sums.multiplicitySums
    for label in [ multiplicitySum.label for multiplicitySum in multiplicitySums if len(multiplicitySum.summands) == 0]: multiplicitySums.pop(label)

    for particle in protare1.listOfProducts().difference(productsInitially):        # The following lines check for any particles that may need to be added.
        protare1.PoPs.add(protare2.PoPs[particle].copy())

    for pid in protare1.listOfProducts():                                           # The following lines check for any masses that may need to be added.
        particle1 = protare1.PoPs[pid]
        try:
            particle1.getMass('amu')
        except:
            pidBase = pid
            if PoPsChemicalElementsMiscModule.hasNucleus(particle1, True): pidBase = particle1.ancestor[0].id
            if pidBase in protare2.PoPs:
                particle2 = protare2.PoPs[pidBase]
                if len(particle2.mass) > 0:
                    particle1.mass.add(particle2.mass[0].copy())
            try:                                                                    # Did it get fixed.
                particle1.getMass('amu')
            except:
                print('WARNING: no mass for particle "%s" in PoPs.' % pid)

if args.writeAlways: protare1.saveToFile( args.outputFile, formatVersion = protare1.format )
