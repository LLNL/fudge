#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

import sys, os, shutil
binDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert( 0, os.path.dirname( binDir ) )

from pqu import PQU as PQUModule
from xData import standards as standardsModule

import site_packages.legacy.toENDF6.toENDF6     # this import adds 'toENDF6' methods to many GNDS classes
import site_packages.legacy.toENDF6.endfFormats as endfFormatsModule

from fudge.legacy.converting import endfFileToGNDS
from fudge.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc

from fudge.gnds import sums as sumsModule
from fudge.gnds.reactionData import crossSection as crossSectionModule

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
        print "%s = %s" % ( key, getattr( args, key ) )

evaluatedStyle = 'eval'
rce = endfFileToGNDS.endfFileToGNDS( args.inputFile, toStdOut = args.verbose, skipBadData = args.skipBadData,
                                   doCovariances = False, verboseWarnings = False, 
                                   printBadNK14 = args.printBadNK14, continuumSpectraFix = args.continuumSpectraFix )


reactionSuite = rce['reactionSuite']
if( args.gnds ) : reactionSuite.saveToFile( 'test.endf6.xml' )

if( args.verbose ) : print 'Performing PREPRO operations.'
temperature = PQUModule.PQU( args.temperature, "%s / k" % reactionSuite.reactions[0].crossSection.domainUnit )
reconStyle = 'recon'

def PREPRO( reaction ) :

    if( args.verbose ) : print '    %s' % reaction

    crossSection = reaction.crossSection
    evaluated = crossSection[evaluatedStyle]
    if( isinstance( evaluated, crossSectionModule.reference ) ) : return

    if( reconStyle in crossSection ) : 
        crossSection.remove( evaluatedStyle )
        recon = crossSection[reconStyle]
        recon.label = evaluatedStyle

    evaluated = crossSection[evaluatedStyle]
    linear = None
    if( isinstance( evaluated, crossSectionModule.regions1d ) ) :
        linear = evaluated.toPointwise_withLinearXYs( accuracy = 1e-3, upperEps = 1e-8 )
    elif( evaluated.interpolation != standardsModule.interpolation.linlinToken ) :
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

if( args.verbose ) : print 'Writing PREPRO ENDF file.'
flags = { 'verbosity' : 0 }
if( args.verbose ) : flags['verbosity'] = 31
outputFileName = os.path.basename( args.inputFile ) + '.prepro'
with open( outputFileName, 'w' ) as fout :
        endfList = reactionSuite.toENDF6( evaluatedStyle, flags, useRedsFloatFormat = args.useRedsFloatFormat ).split('\n')
        endfTxt='\n'.join([FIRSTLINE] + endfList[1:])
        fout.write( endfTxt )

if( args.debug ) : 
    if( args.verbose ) : print 'Debugging.'
    endfFileToGNDS.endfFileToGNDS( args.inputFile, toStdOut = args.verbose, skipBadData = args.skipBadData,
                                 doCovariances = False, verboseWarnings = False,
                                 printBadNK14 = args.printBadNK14, continuumSpectraFix = args.continuumSpectraFix )
