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

import sys

import fudge
from pqu import PQU
from fudge.gnd import alias

import xData.standards as standardsModule

from fudge.gnd.productData.distributions import unspecified as unspecifiedModule

from fudge.legacy.converting.ENDFToGND import endfFileToGNDMisc, ENDF_ITYPE_0
from fudge.legacy.converting.ENDFToGND import ENDF_ITYPE_3, ENDF_ITYPE_6
from fudge.legacy.converting import endf_endl as endf_endlModule
import toGNDMisc

class logFiles :

    def __init__( self, toStdOut = False, toStdErr = False, logFile = None, defaultIsStderrWriting = True ) :

        def addIfDifferent( newOutput ) :

            if( newOutput.isatty( ) ) :
                for output in self.outputs :
                    if( output.isatty( ) ) : return( False )
            return( True )

        self.defaultIsStderrWriting = defaultIsStderrWriting
        self.outputs = []
        if( toStdOut ) : self.outputs.append( sys.stdout )
        if( toStdErr and addIfDifferent( sys.stderr ) ) : self.outputs.append( sys.stderr )
        if( logFile ) :
            if( type( logFile ) == type( '' ) ) :
                logFile = open( logFile, 'w' )
            elif( type( logFile ) != type( sys.stdout ) ) :
                raise Exception( 'logFile type = %s is no supported' % ( type( logFile ) ) )
            if( addIfDifferent( logFile ) ) : self.outputs.append( logFile )

    def write( self, msg, stderrWriting = None ) :

        if( stderrWriting is None ) : stderrWriting = self.defaultIsStderrWriting
        for output in self.outputs :
            if( ( output == sys.stderr ) and ( not( stderrWriting ) ) ) : continue
            output.write( msg )

def endfFileToGND( fileName, xenslIsotopes = None, useFilesQAlways = True, singleMTOnly = None,
        MTs2Skip = None, parseCrossSectionOnly = False,
        toStdOut = True, toStdErr = True, logFile = None, skipBadData = False, doCovariances = True,
        verboseWarnings = False, deprecatedOptions = None ) :

    logs = logFiles( toStdOut = toStdOut, toStdErr = toStdErr, logFile = logFile, defaultIsStderrWriting = False )
    header, MAT, MTDatas = endfFileToGNDMisc.parseENDFByMT_MF( fileName, logFile = logs )

    styleName = 'eval'
    reconstructedStyleName = 'recon'

    if MTs2Skip is None: MTs2Skip = []
    if deprecatedOptions is None: deprecatedOptions = {}

    targetZA, targetMass, LRP, LFI, NLIB, NMOD = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MTDatas[451][1][0], logFile = logs )
    targetZA = int( targetZA )      # Target's ZA
    LRP = int( LRP )            # Resonance parameter data info
    LFI = int( LFI )            # Is fission present
    NLIB = int( NLIB )          # What library (e.g., 0 = ENDF/B
    NMOD = int( NMOD )          # Version modification flag
    isNaturalTarget = ( targetZA % 1000 ) == 0

    targetExcitationEnergy, STA, LIS, LISO, dummy, NFOR = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MTDatas[451][1][1], logFile = logs )
    STA = int( STA )            # Is nucleus unstable
    LIS = int( LIS )            # Excitation number
    LISO = int( LISO )          # Isomeric state number
    NFOR = int( NFOR )          # Must be 6 for ENDF/B6 format
    if (NFOR != 6):
        print ("    WARNING: endfFileToGND only supports ENDF-6 format. This file has unsupported NFOR=%d" % NFOR)

    projectileMass, dummy, LREL, dummy, NSUB, NVER = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MTDatas[451][1][2], logFile = logs )
    NSUB = int( NSUB )          # 10 * ZA + iType for projectile
    NVER = int( NVER )          # Evaluation version number
    LREL = int( LREL )          # Evaluation sub-version number
    IPART, ITYPE = NSUB / 10, NSUB % 10
    projectileZA = IPART
    if( projectileZA == 11 ) : projectileZA = 9

    targetTemperature, dummy, LDRZ, dummy, NWD, NXC = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MTDatas[451][1][3], logFile = logs )
    LDRZ = int( LDRZ )          # Primary or special evaluation of this material
    NWD = int( NWD )            # 
    NXC = int( NXC )            #

    library = {
        0: "ENDF/B",
        1: "ENDF/A",
        2: "JEFF",
        3: "EFF",
        4: "ENDF/B (HE)",
        5: "CENDL",
        6: "JENDL",
        21: "SG-23",
        31: "INDL/V",
        32: "INDL/A",
        33: "FENDL",
        34: "IRDF",
        35: "BROND (IAEA version)",
        36: "INGDB-90",
        37: "FENDL/A",
        41: "BROND",
    }.get( NLIB, 'Unknown' )
    libraryVersion = "%d.%d.%d" % ( NVER, LREL,NMOD )

    transportables = [ 'n', 'gamma' ]                          # ???? Check this with Gerry?
    info = toGNDMisc.infos( styleName, xenslIsotopes, transportables = transportables )
    info.doRaise = []
    try : 
        Date = endfFileToGNDMisc.getENDFDate( MTDatas[451][1][4][22:33] )
    except Exception as e :
        info.doRaise.append( str(e) )
        import datetime
        Date = datetime.datetime.today().strftime("%Y-%m-%d")
    author = MTDatas[451][1][4][33:66]

    evaluatedStyle = fudge.gnd.styles.evaluated( styleName, Date,
            PQU.PQU( PQU.pqu_float.surmiseSignificantDigits( targetTemperature ), 'K' ),
            library, libraryVersion )

    info.Date = Date
    info.verboseWarnings = verboseWarnings
    info.ignoreBadNK14 = False
    info.continuumSpectraFix = False
    deprecatedOptionList = [ 'ignoreBadNK14', 'continuumSpectraFix' ]
    for options in deprecatedOptions :
        if( options not in deprecatedOptionList ) : raise Exception( 'invalid deprecated option "%s"' % options )
        setattr( info, options, deprecatedOptions[options] )
    info.logs = logs
    info.MAT = MAT
    info.LRP = LRP
    info.NLIB = NLIB
    info.evaluation = library
    info.NMOD = NMOD
    info.NVER = NVER
    info.LREL = LREL
    info.evaluatedStyle = evaluatedStyle
    info.reconstructedStyle = reconstructedStyleName
    info.reconstructedAccuracy = 0.001

    info.projectile = { 0 : 'g', 1 :  'n', 11 : 'e-', 1001 : 'H1', 1002 : 'H2', 1003 : 'H3', 2003 : 'He3', 2004 : 'He4' }[IPART]
    info.projectileZA = projectileZA

    info.ZA_AWRMasses = {}
    info.ZA_AWRMasses[projectileZA] = { projectileMass : 1 }
    info.ZA_AWRMasses[targetZA] = { targetMass : 1 }
    info.ZAMasses[projectileZA] = info.masses.getMassFromZA( projectileZA )
    info.ZAMasses[targetZA] = targetMass * info.masses.getMassFromZA( 1 )
    info.ZAMasses[1] = info.masses.getMassFromZA( 1 ) # always need neutron mass in table
    
    info.MF12_LO2 = {}
    info.AWR_mode = None
#    info.AWR_mode = open( 'AWR_mode.out', 'w' )  # Used by BRB to print masses for testing.

    if( ITYPE == 6 )  :
        if( NSUB not in [ 6 ] ) : raise Exception( 'For ITYPE = %d, invalid NSUB = %s' % ( ITYPE, NSUB ) )
        return( ENDF_ITYPE_6.ITYPE_6( targetZA / 1000, MTDatas, info, verbose = True ) )

    projectile = toGNDMisc.getTypeNameGamma( info, projectileZA )
    info.level = targetExcitationEnergy
    levelIndex = None
    if( LIS  != 0 ) : levelIndex = LIS
    if( ITYPE in [ 0, 9 ] )  :
        target = toGNDMisc.getTypeNameGamma( info, targetZA, level = info.level, levelIndex = levelIndex )
        if( ( STA != 0 ) and not isinstance( target, fudge.gnd.xParticle.nuclearLevel ) ) : target.attributes['unstable'] = True
    elif( ITYPE == 3 )  :
        elementSymbol = fudge.particles.nuclear.elementSymbolFromZ( targetZA / 1000 )
        target = fudge.gnd.xParticle.element( elementSymbol )
    else :
        raise Exception( "Unsupported ITYPE = %s" % ITYPE )

    ZA2, MAT2 = endf_endlModule.ZAAndMATFromParticleName( target.name )
    MAT2 += LISO
    if( MAT2 != MAT ) : info.logs.write( "       WARNING: ENDF MAT = %s not as expected (i.e., %s).\n" % \
            ( MAT, MAT2 ), stderrWriting = True )

    documentation = fudge.gnd.documentation.documentation( 'endfDoc', '\n'.join( MTDatas[451][1][4:4+NWD] ) )
    reactionSuite = fudge.gnd.reactionSuite.reactionSuite( projectile, target, 
            particleList = info.particleList, style = evaluatedStyle, documentation = documentation, MAT = MAT )
    if( LISO != 0 ) :
        targetBaseName = target.name.split( '_' )[0]
        aliasName = alias.aliases.nuclearMetaStableName( targetBaseName, LISO )
        reactionSuite.addNuclearMetaStableAlias( targetBaseName, target.name, LISO )
    info.setReactionSuite( reactionSuite )
    info.target = reactionSuite.target
    info.targetZA = targetZA
    info.targetLevel = LIS

    MTDatas[451][1] = MTDatas[451][1][:4+NWD]

    covarianceSuite = None
    if( ( ITYPE == 0 ) or ( ITYPE == 9 ) ) :
        doRaise = NSUB not in { 0 : [ 0, 10, 10010, 10020, 10030, 20030, 20040 ], 9 : [ 19 ] }[ITYPE]
        if( doRaise ) : raise Exception( 'For ITYPE = %d, invalid NSUB = %s' % ( ITYPE, NSUB ) )
        covarianceSuite = ENDF_ITYPE_0.ITYPE_0( MTDatas, info, reactionSuite, singleMTOnly, MTs2Skip, parseCrossSectionOnly, doCovariances )
    elif( ITYPE == 3 )  :
        if( NSUB not in [ 3, 113 ] ) :
            raise Exception( 'For ITYPE = %d, invalid NSUB = %s' % ( ITYPE, NSUB ) )
        ENDF_ITYPE_3.ITYPE_3( MTDatas, info, reactionSuite, singleMTOnly, parseCrossSectionOnly, verbose = True )

    if( len( info.doRaise ) > 0 and not skipBadData ) :
        info.logs.write( '\nRaising due to following errors:\n' )
        for err in info.doRaise : info.logs.write( err + '\n' )
        raise Exception( 'len( info.doRaise ) > 0' )

    for reaction in reactionSuite.reactions : addUnspecifiedDistributions( info, reaction.outputChannel )
    for production in reactionSuite.productions : addUnspecifiedDistributions( info, production.outputChannel )

    return( { 'reactionSuite' : reactionSuite, 'covarianceSuite' : covarianceSuite, 'errors' : info.doRaise, 'info':info } )

def addUnspecifiedDistributions( info, outputChannel ) :

    if( outputChannel is None ) : return
    for product in outputChannel :
        if( len( product.distribution ) == 0 ) :
            product.distribution.add( unspecifiedModule.form( info.style, productFrame = standardsModule.frames.labToken ) )
        addUnspecifiedDistributions( info, product.decayChannel )

if( __name__ == '__main__' ) :

    rce = endfFileToGND( fileName = sys.argv[1] )
    x, c = rce['reactionSuite'], rce['covarianceSuite']
    f = open( 'test.xml', 'w' )
    f.write( '\n'.join( x.toXMLList( ) + [ '' ] ) )
    f.close( )
    if( c is not None ) : # covariances
        f = open( 'test-covar.xml', 'w' )
        f.write( '\n'.join( c.toXMLList( ) + [ '' ] ) )
        f.close()
