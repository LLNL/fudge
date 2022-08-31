#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys, os

binDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert( 0, os.path.dirname( binDir ) )

from pqu import PQU as PQUModule
from xData import enums as xDataEnumsModule

import brownies.legacy.toENDF6.toENDF6     # this import adds 'toENDF6' methods to many GNDS classes
from brownies.legacy.converting import endfFileToGNDS

from fudge.reactionData import crossSection as crossSectionModule

from argparse import ArgumentParser

FIRSTLINE=" $Rev:: 700      $  $Date:: 2016-01-13#$                             1 0  0    0"

usage = """
A FUDGE script to replicate PREPRO. In particular the PREPRO modules recent, linear and sigma1.
Each of these modules does the follow::

    +--------+-------------------------------------------------------------------+
    | recent | Reconstructs cross sections from resonance parameters.            |
    +--------+-------------------------------------------------------------------+
    | linear | converts all cross sections to 'lin-lin' interpolation.           |
    +--------+-------------------------------------------------------------------+
    | sigma1 | Heat cross section to a specified temperature.                    |
    +--------+-------------------------------------------------------------------+

This module automatically performs the equivalent for recent and linear.
Heating is performs if the '-t' option is used to specify a temperature.
"""

__doc__ = usage

parser = ArgumentParser( description = usage )
parser.add_argument( 'inputFile', type = str,
        help = 'ENDF-6 file to translate to GNDS and write back to ENDF with PREPRO-like changes.')
parser.add_argument( '-s', '--skipBadData', default = False, action = 'store_true',
        help = 'Recover from format errors if possible when converting ENDF to GNDS.' )
parser.add_argument( '-t', '--temperature', default = 0, action = 'store', type = float,
        help = '''Set the temperature to heat cross section to. "k" times temperature (i.e., "k * T") 
        unit is assumed to be the same as the incident energy unit (e.g., "eV" for an ENDF file).''' )
parser.add_argument( '--printBadNK14', default = False, action = 'store_true',
        help = 'If true print warning for NK mismatch between MF 12 and 14.' )
parser.add_argument( '--continuumSpectraFix', default = False, action = 'store_true',
        help = 'Skip unnormalizeable continuum gamma distributions.' )
parser.add_argument( '-v', '--verbose', default = False, action = 'store_true',
        help = 'enable verbose output.' )
parser.add_argument( '-r', '--useRedsFloatFormat', default = True, action = 'store_false',
        help = '''If True (default) write floats to the ENDF file using Red Cullen's high precision format.''' )
parser.add_argument( '-d', '--debug', default = False, action = 'store_true',
        help = 'Reread the outputted ENDF file to see if it is okay. This is for debugging.' )
parser.add_argument( '-g', '--gnds', default = False, action = 'store_true',
        help = 'Write out GNDS file. This is for debugging.' )

args = parser.parse_args()

if( args.verbose ) :
    for key in dir( args ) :
        if( key[0] == '_' ) : continue
        print( "%s = %s" % ( key, getattr( args, key ) ) )

evaluatedStyle = 'eval'
rce = endfFileToGNDS.endfFileToGNDS(args.inputFile, toStdOut = args.verbose, skipBadData = args.skipBadData,
                                    doCovariances = False, verboseWarnings = False,
                                    printBadNK14 = args.printBadNK14, continuumSpectraFix = args.continuumSpectraFix)


reactionSuite = rce['reactionSuite']
if( args.gnds ) : reactionSuite.saveToFile( 'test.endf6.xml' )

if( args.verbose ) : print( 'Performing PREPRO operations.' )
temperature = PQUModule.PQU( args.temperature, "%s / k" % reactionSuite.reactions[0].crossSection.domainUnit )
reconStyle = 'recon'

def PREPRO( reaction ) :

    if( args.verbose ) : print( '    %s' % reaction )

    crossSection = reaction.crossSection
    evaluated = crossSection[evaluatedStyle]
    if( isinstance( evaluated, crossSectionModule.Reference ) ) : return

    if( reconStyle in crossSection ) : 
        crossSection.remove( evaluatedStyle )
        recon = crossSection[reconStyle]
        recon.label = evaluatedStyle

    evaluated = crossSection[evaluatedStyle]
    linear = None
    if isinstance(evaluated, crossSectionModule.Regions1d):
        linear = evaluated.toPointwise_withLinearXYs( accuracy = 1e-3, upperEps = 1e-8 )
    elif evaluated.interpolation != xDataEnumsModule.Interpolation.linlin:
        linear = evaluated.toPointwise_withLinearXYs( accuracy = 1e-3, upperEps = 1e-8 )
    if( linear is not None ) :
        crossSection.remove( evaluatedStyle )
        linear.label = evaluatedStyle
        crossSection.add( linear )

    if( args.temperature > 0 ) :
        heated = crossSection.heat( temperature, crossSection.domainMin )
        crossSection.remove( evaluatedStyle )
        heated.label = evaluatedStyle
        crossSection.add( heated )

for reaction in reactionSuite.reactions : PREPRO( reaction )
for sum in reactionSuite.sums.crossSectionSums : PREPRO( sum )

if( args.verbose ) : print( 'Writing PREPRO ENDF file.' )
flags = { 'verbosity' : 0 }
if( args.verbose ) : flags['verbosity'] = 31
outputFileName = os.path.basename( args.inputFile ) + '.prepro'
with open( outputFileName, 'w' ) as fout :
        endfList = reactionSuite.toENDF6( evaluatedStyle, flags, useRedsFloatFormat = args.useRedsFloatFormat ).split('\n')
        endfTxt='\n'.join([FIRSTLINE] + endfList[1:])
        fout.write( endfTxt )

if( args.debug ) : 
    if( args.verbose ) : print( 'Debugging.' )
    endfFileToGNDS.endfFileToGNDS(args.inputFile, toStdOut = args.verbose, skipBadData = args.skipBadData,
                                  doCovariances = False, verboseWarnings = False,
                                  printBadNK14 = args.printBadNK14, continuumSpectraFix = args.continuumSpectraFix)
