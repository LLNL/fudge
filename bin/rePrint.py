#! /usr/bin/env python

# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
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
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

import sys, os, shutil
binDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert( 0, os.path.dirname( binDir ) )

import site_packages.legacy.toENDF6.toENDF6     # this import adds 'toENDF6' methods to many GND classes
import site_packages.legacy.toENDF6.endfFormats as endfFormatsModule

from fudge.legacy.converting import endfFileToGND
from fudge.legacy.converting.ENDFToGND import endfFileToGNDMisc
from fudge.processing import processingInfo
from fudge.core.utilities import subprocessing

def process_args( ) :       # see https://docs.python.org/2/howto/argparse.html

    from argparse import ArgumentParser
    usage = """
>rePrint.py ENDF-formatted-file.endf [MT]

rePrint translates an ENDF file to the new GND format. If MT is given, only translate the specified MT"""

    parser = ArgumentParser( description=usage )
    parser.add_argument( 'inputFile', type=str, help='ENDF-6 file to translate to/from GND')
    parser.add_argument( 'MT', type=int, nargs="?", default=None )
    parser.add_argument( '-v', '--verbose', default = False, action="store_true",
                        help = "enable verbose output")
    parser.add_argument( '--verboseWarnings', default = False, action="store_true",
                        help = "increase detail for warnings" )
    parser.add_argument( '-x', '--xmlOnly', default=False, action="store_true",
                        help = "only translate ENDF to GND/Python and write xml files. Do not re-write endf file" )
    parser.add_argument( '-t', '--translateOnly', default=False, action="store_true",
                        help = "only translate ENDF to GND/Python, do not write any files" )
    parser.add_argument( "-c", '--checkCovars', default = False, action = "store_true",
                        help = "enable covariance checking" )
    parser.add_argument( "-s", "--skipBadData", default = False, action = "store_true",
                        help = "Recover from format errors if possible" )
    parser.add_argument( "-o", "--outline", default = False, action = "store_true",
                        help = "Set outline option to True for GND/XML output" )
    parser.add_argument( "--skipCovariances", default = False, action = "store_true",
                        help = "Do not translate covariance data" )

    # deprecated options:
    deprecated = parser.add_argument_group('deprecated')
    deprecated.add_argument( "--ignoreBadNK14", default = False, action = "store_true",
                        help = "Ignore NK mismatch between MF 12 and 14" )
    deprecated.add_argument( "--continuumSpectraFix", default = False, action = "store_true",
                        help = "Skip unnormalizeable continuum gamma distributions" )

    return parser.parse_args()

subprocessing.deleteFilesUsingGlob( 'test.endf6*' )

args = process_args()

counter = 0
def removeMFMT( MFMT, MFAlso = True, MTAlso = True ) :

    global counter
    if( counter == 0 ) : shutil.copy2( 'test.endf6.orig.noLineNumbers', 'test.endf6.orig.noLineNumbers%d' % counter )
    if( MTAlso and MFAlso ) :
        cmd = 'removeMFMT.py'
    elif( MTAlso ) :
        cmd = 'removeMTs.py'
    else :
        cmd = 'removeMFs.py'
    subprocessing.executeCommand( [ 'python', os.path.join( binDir, cmd ), 'test.endf6.orig.noLineNumbers%d' % counter,
            'test.endf6.orig.noLineNumbers%d' % ( counter + 1 ) ] + MFMT.split( ) )
    counter += 1


if( args.MT is None ) :
    shutil.copy2( args.inputFile, 'test.endf6.orig2' )
else :
    header, MAT, MTDatas = endfFileToGNDMisc.parseENDFByMT_MF( args.inputFile, stripMATMFMTCount = False )
    f = open( 'test.endf6.orig2', 'w' )
    f.write( header + '\n' )
    MTs = [ 451, args.MT ]
    if( args.MT == 18 ) :
        MTsp = [ 452, 455, 456, 458 ]
        for MTp in MTsp :
            if( MTp in MTDatas ) : MTs.append( MTp )
    MFLines = {}
    for MTp in MTs :
        if( MTp not in MTDatas ) : continue
        MFs = sorted( MTDatas[MTp] )
        for MF in MFs : 
            if( MF != 1 ) : continue
            if( MF not in MFLines ) : MFLines[MF] = []
            MFLines[MF] += MTDatas[MTp][MF] + [ endfFormatsModule.endfSENDLine( MAT, MF ) ]
    MTs.sort( )
    for MTp in MTs :
        if( MTp not in MTDatas ) : continue
        MFs = sorted( MTDatas[MTp] )
        for MF in MFs : 
            if( MF == 1 ) : continue
            if( MF not in MFLines ) : MFLines[MF] = []
            MFLines[MF] += MTDatas[MTp][MF] + [ endfFormatsModule.endfSENDLine( MAT, MF ) ]
    MFs = sorted( MFLines.keys( ) )
    for MF in MFs :
        f.write( '\n'.join( MFLines[MF] ) + '\n' )
        f.write( endfFormatsModule.endfFENDLine( MAT ) + '\n' )
    f.write( endfFormatsModule.endfMENDLine( ) + '\n' )
    f.write( '%s' % endfFormatsModule.endfTENDLine( ) + '\n' )
    f.close( )

flags = processingInfo.tempInfo( )
flags['verbosity'] = 0
subprocessing.executeCommand( [ 'python', os.path.join( binDir, 'reForm.py' ), 'test.endf6.orig2', 'test.endf6.orig' ] )
subprocessing.deleteFilesUsingGlob( 'test.endf6.orig2' )
style = 'eval'
deprecated = {'ignoreBadNK14': args.ignoreBadNK14,  'continuumSpectraFix': args.continuumSpectraFix}
rce = endfFileToGND.endfFileToGND( args.inputFile, singleMTOnly = args.MT, toStdOut = args.verbose,
                                   skipBadData = args.skipBadData, doCovariances = not( args.skipCovariances ),
                                   verboseWarnings = args.verboseWarnings, deprecatedOptions = deprecated )
errs = rce['errors']
if( 'reactionSuite' in rce ) :
    reactions, covariances = rce['reactionSuite'], rce['covarianceSuite']
    if( args.checkCovars ) :
        sys.stderr.write( ''.join(
            [ "%s: %s" % ( os.path.split( args.inputFile )[-1], warning ) for warning in covariances.check( ) ]
        ) )
    if( args.translateOnly ) : sys.exit( len( errs ) )

    reactions.saveToFile( 'test.endf6.xml', outline = args.outline )
    if( covariances is not None ) :
        covariances.saveToFile( 'test.endf6-covar.xml', outline = args.outline )
    if( args.xmlOnly ) : sys.exit( len( errs ) )

    if( args.verbose ) : flags['verbosity'] = 31
    with open( 'test.endf6', 'w' ) as fout:
        fout.write( reactions.toENDF6( style, flags, covarianceSuite = covariances ) )

    subprocessing.executeCommand( [ 'python', os.path.join( binDir, 'noLineNumbers.py' ), 'test.endf6', 'test.endf6.noLineNumbers' ] )
    subprocessing.executeCommand( [ 'python', os.path.join( binDir, 'noLineNumbers.py' ), 'test.endf6.orig', 'test.endf6.orig.noLineNumbers' ] )
    removeMFMT( '29 30', MTAlso = False )
    removeMFMT( '201 202 203 204 205 206 207', MFAlso = False )
    removeMFMT( '1 460' )
    shutil.copy2( 'test.endf6.orig.noLineNumbers%d' % counter, 'test.endf6.orig.noLineNumbers.clean' )
    subprocessing.executeCommand( [ 'python', os.path.join( binDir, 'fixZeroLeaders.py' ), 'test.endf6.orig.noLineNumbers.clean',
        'test.endf6.orig.noLineNumbers.cleanAndFixed' ] )
    subprocessing.deleteFilesUsingGlob( 'test.endf6.orig.noLineNumbers[0-9]*' )
    if( True ) : subprocessing.executeCommand( [ 'rm', '-f', 'test.endf6', 'test.endf6.orig', 'test.endf6.orig.noLineNumbers.clean' ] )
elif( 'element' in rce ) :
    element = rce['element']
    f = open( 'test.endf6.xml', 'w' )
    f.write( '\n'.join( element.toXMLList( outline = args.outline ) + [ '' ] ) )
    f.close( )
else :
    raise Exception( 'Unsupported return from endfFileToGND.endfFileToGND' )

sys.exit( len( errs ) )
