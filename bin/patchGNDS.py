#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
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
"""

from argparse import ArgumentParser

from fudge import styles as stylesModule
from fudge import reactionSuite as reactionSuiteModule

reconstructResonancesLabelDefault = 'recon'

parser = ArgumentParser( description = description )
parser.add_argument( "label", type = str,                                               help = """The label for the new "evaluated" style.""" )
parser.add_argument( "version", type = str,                                             help = """The library version for the new "evaluated" style.""" )
parser.add_argument( "patchFile",                                                       help = """The name of the path file.""" )
parser.add_argument( "outputFile",                                                      help = """The name of the file to write the patched GNDS data to.""" )
parser.add_argument( "--library", type = str, default = "",                             help = """The library for the new "evaluated" style. If not present, use the latest in the file to patch.""" )
parser.add_argument( "--reconstructResonances", action = 'store_true', default = False, help = """If present resonances are reconstructed.""" )
parser.add_argument( "--reconstructResonancesLabel", action = 'store', default = reconstructResonancesLabelDefault,
                                                                                        help = """Label to use if reconstructing resonances. Default = "%s".""" % reconstructResonancesLabelDefault )

args = parser.parse_args( )

fIn = open( args.patchFile )
lines = fIn.readlines( )
fIn.close( )

if( ( len( lines ) - 2 ) % 3 != 0 ) : raise Exception( "Input file does not have the correct number of lines to be a valid patch file." )

def checkFileLine( index ) :

    line = lines[index]
    if( ":" not in line ) : raise Exception( "Invalid file line at line %d is not valid." % ( index + 1 ) )
    tag, fileName = line.split( ":" )
    if( "FILE%s" % ( index + 1 ) != tag ) : raise Exception( "Invalid file line at line %d is not valid." % ( index + 1 ) )

    return( fileName.strip( ) )

protare1 = reactionSuiteModule.readXML( checkFileLine( 0 ) )
protare2 = reactionSuiteModule.readXML( checkFileLine( 1 ) )

if( protare1.domainUnit != protare2.domainUnit ) : protare1.convertUnits( { protare1.domainUnit : protare2.domainUnit } )

chains = protare1.styles.chains( True )
if( len( chains ) != 1 ) : raise Exception( "Only one styles chain is supported." )

chain = chains[0]
for style in chain :
    if( isinstance( style, stylesModule.evaluated ) ) : break
if( not( isinstance( style, stylesModule.evaluated ) ) ) : raise Exception( "At least one evaluated style must be present in file to patch." )

if( len( lines ) > 2 ) :

    if( args.label in protare1.styles ) : raise ValueError( "Label already in file to patch." )
    if( args.library == "" ) : args.library = style.library
    evaluateStyle = stylesModule.evaluated( args.label, style.label, None, None, args.library, args.version )
    if( style.library == evaluateStyle.library ) :
        if( style.version == evaluateStyle.version ) : raise ValueError( "Library and version must be different than one in patch file." )
    protare1.styles.add( evaluateStyle )

    links = protare1.findLinks( )

    sumFixes = []
    for index in range( 2, len( lines ), 3 ) :
        line = lines[index]
        if( ":" not in line ) : raise Exception( "Input tag at line %d is not valid." % ( 3 * index ) )
        tag = line.split( ":" )[0]

        toXLink = lines[index+1].strip( )
        fromXLink = lines[index+2].strip( )

        if( tag == "Distribution unspecified - 1" ) :
            fromComponent = protare2.followXPath( fromXLink )
            toComponent = protare1.followXPath( toXLink )
            for key in fromComponent.keys( ) : toComponent.remove( key )
            toComponent.add( fromComponent.evaluated )
        elif( tag == "reactions missing - 1" ) :
            toReaction = protare2.followXPath( fromXLink )
            toReaction.amendForPatch( style.label, evaluateStyle.label )
            protare1.reactions.add( toReaction )
        elif( tag == "reactions missing - 2" ) :
            sumFixes.append( [ 'removed', toXLink ] )
            protare1.reactions.pop( protare1.followXPath( toXLink ).label )
        elif( 'non-constant Q - 2' in tag ) :
            toComponent = protare1.followXPath( toXLink ).ancestor
            fromComponent = protare2.followXPath( fromXLink )
            toComponent.replace( fromComponent )
            fromComponent.label = evaluateStyle.label
        elif( 'More than one product with same id' ) :
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
        else :
            print( """Unsupported diff tag "%s".""" % tag )

crossSectionReconstructedStyles = protare1.styles.getStylesOfClass( stylesModule.crossSectionReconstructed )
if( len( crossSectionReconstructedStyles ) > 0 ) :
    if( len( crossSectionReconstructedStyles ) > 1 ) : raise Exception( "Multiple cross section reconstructed styles present which is not supported." )
    crossSectionReconstructedStyle = crossSectionReconstructedStyles[0]
    protare1.removeStyle( crossSectionReconstructedStyle.label )

if( args.reconstructResonances ) :
    crossSectionReconstructedStyle = stylesModule.crossSectionReconstructed( args.reconstructResonancesLabel, args.label )
    protare1.reconstructResonances( crossSectionReconstructedStyle )

protare1.fixLinks( )

protare1.saveToFile( args.outputFile, formatVersion = protare1.format )
