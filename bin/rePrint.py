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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

import sys, os, shutil
binDir = os.path.dirname( os.path.abspath( __file__ ) )
sys.path.insert(0, os.path.dirname( binDir ) )

from fudge.legacy.converting import endfFileToGND, endfFileToGNDMisc, endfFormats
from fudge.processing import processingInfo
from fudge.core.utilities import subprocessing

def process_args( ) :       # see http://docs.python.org/lib/optparse-tutorial.html

    from optparse import OptionParser
    usage = """
>rePrint.py ENDF-formatted-file.endf [MT]

rePrint translates an ENDF file to the new GND format. If MT is given, only translate the specified MT"""

    parser = OptionParser( usage )
    parser.set_defaults( verbose = True )
    parser.add_option( "-v", action = "store_true",  dest = "verbose", help = "enable verbose output" )
    parser.add_option( "-q", action = "store_false", dest = "verbose", help = "disable verbose output" )
    parser.add_option( "-c", action = "store_true",  dest = "checkCovars", help = "enable covariance checking" )
    parser.add_option( "--skipBadData", action = "store_true", dest = "skipBadData", help = "Recover from format errors if possible" )
    # for backwards compatibility, still keep the --recoverOnError option:
    parser.add_option( "--recoverOnError", action = "store_true", dest = "skipBadData", help = "Recover from format errors if possible" )

    return parser.parse_args()

subprocessing.deleteFilesUsingGlob( 'test.endf6*' )

Options, argv = process_args()

try :
    MT = int( argv[1] )
except :
    MT = None

if( MT is None ) :
    shutil.copy2( argv[0], 'test.endf6.orig2' )
else :
    header, MAT, MTDatas = endfFileToGNDMisc.parseENDFByMT_MF( argv[0], stripMATMFMTCount = False )
    f = open( 'test.endf6.orig2', 'w' )
    f.write( header + '\n' )
    MTs = [ 451, MT ]
    if( MT == 18 ) :
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
            MFLines[MF] += MTDatas[MTp][MF] + [ endfFormats.endfSENDLine( MAT, MF ) ]
    MTs.sort( )
    for MTp in MTs :
        if( MTp not in MTDatas ) : continue
        MFs = sorted( MTDatas[MTp] )
        for MF in MFs : 
            if( MF == 1 ) : continue
            if( MF not in MFLines ) : MFLines[MF] = []
            MFLines[MF] += MTDatas[MTp][MF] + [ endfFormats.endfSENDLine( MAT, MF ) ]
    MFs = sorted( MFLines.keys( ) )
    for MF in MFs :
        f.write( '\n'.join( MFLines[MF] ) + '\n' )
        f.write( endfFormats.endfFENDLine( MAT ) + '\n' )
    f.write( endfFormats.endfMENDLine( ) + '\n' )
    f.write( '%s' % endfFormats.endfTENDLine( ) + '\n' )
    f.close( )
subprocessing.executeCommand( [ 'python', os.path.join( binDir, 'reForm.py' ), 'test.endf6.orig2', 'test.endf6.orig' ] )
subprocessing.deleteFilesUsingGlob( 'test.endf6.orig2' )
rce = endfFileToGND.endfFileToGND( argv[0], singleMTOnly = MT, toStdOut = Options.verbose, skipBadData = Options.skipBadData )
r, c = rce['reactionSuite'], rce['covarianceSuite']
errs = rce['errors']
if( Options.checkCovars ) :
    sys.stderr.write( ''.join( [ "%s: %s" % ( os.path.split( argv[0] )[-1], warning ) for warning in c.checkCovariances( ) ] ) )

flags = processingInfo.tempInfo( )
flags['verbosity'] = 0
f = open( 'test.endf6.xml', 'w' )
f.write( '\n'.join( r.toXMLList( flags ) + [ '' ] ) )
f.close( )
if c:
    f = open( 'test.endf6-covar.xml', 'w' )
    f.write( '\n'.join( c.toXMLList( flags ) + [ '' ] ) )
    f.close()

if( Options.verbose ) : flags['verbosity'] = 31
f = open( 'test.endf6', 'w' )
f.write( r.toENDF6( flags, covarianceSuite=c ) )
f.close( )

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
sys.exit( len( errs ) )
