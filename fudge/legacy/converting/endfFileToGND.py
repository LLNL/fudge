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

import sys

import fudge
from pqu import PQU as PQUModule
from fudge.gnd import alias

from xData import standards as standardsModule

from PoPs import database as databasePoPsModule
from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import halflife as halflifeModule
from PoPs.quantities import mass as massModule
from PoPs.families import nuclearLevel as nuclearLevelModule
from PoPs.groups import isotope as isotopeModule
from PoPs.groups import chemicalElement as chemicalElementModule
from PoPs import alias as PoPsAliasModule

from fudge.gnd import physicalQuantity as physicalQuantityModule
from fudge.gnd import reactionSuite as reactionSuiteModule
from fudge.gnd import styles as stylesModule
from fudge.gnd.productData.distributions import unspecified as unspecifiedModule

from fudge.legacy.converting.ENDFToGND import endfFileToGNDMisc, ENDF_ITYPE_0
from fudge.legacy.converting.ENDFToGND import ENDF_ITYPE_2, ENDF_ITYPE_3, ENDF_ITYPE_4, ENDF_ITYPE_5, ENDF_ITYPE_6
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


def readMF1MT451( _MAT, _MTDatas, styleName='eval', logFile=None, verboseWarnings=False, **kwargs ):

    # We're going to save everything in the info instance
    info = toGNDMisc.infos(styleName)
    info.doRaise = []
    info.logs = logFile
    info.verboseWarnings = verboseWarnings

    # Line #1
    targetZA, targetMass, LRP, LFI, NLIB, NMOD = \
        endfFileToGNDMisc.sixFunkyFloatStringsToFloats( _MTDatas[451][1][0], logFile=logFile )
    targetZA = int(targetZA)  # Target's ZA
    LRP = int(LRP)            # Resonance parameter data info
    LFI = int(LFI)            # Is fission present
    NLIB = int(NLIB)          # What library (e.g., 0 = ENDF/B
    NMOD = int(NMOD)          # Version modification flag
    isNaturalTarget = (targetZA % 1000) == 0

    # Line #2
    targetExcitationEnergy, STA, LIS, LISO, dummy, NFOR = \
        endfFileToGNDMisc.sixFunkyFloatStringsToFloats( _MTDatas[451][1][1], logFile=logFiles )
    STA = int(STA)            # Is nucleus unstable
    LIS = int(LIS)            # Excitation number
    LISO = int(LISO)          # Isomeric state number
    NFOR = int(NFOR)          # Must be 6 for ENDF-6 format
    if (NFOR != 6):
        print ("    WARNING: endfFileToGND only supports ENDF-6 format. This file has unsupported NFOR=%d" % NFOR)

    # Line #3
    projectileMass, EMAX, LREL, dummy, NSUB, NVER = \
        endfFileToGNDMisc.sixFunkyFloatStringsToFloats( _MTDatas[451][1][2], logFile=logFiles )
    NSUB = int(NSUB)  # 10 * ZA + iType for projectile
    NVER = int(NVER)  # Evaluation version number
    LREL = int(LREL)  # Evaluation sub-version number
    IPART, ITYPE = NSUB / 10, NSUB % 10
    projectileZA = IPART
    if (projectileZA == 11): projectileZA = 9

    # Line #4
    targetTemperature, dummy, LDRZ, dummy, NWD, NXC = \
        endfFileToGNDMisc.sixFunkyFloatStringsToFloats( _MTDatas[451][1][3], logFile=logFiles )
    LDRZ = int(LDRZ)  # Primary or special evaluation of this material
    NWD = int(NWD)    #
    NXC = int(NXC)    #

    # Save the library name and version
    info.library = {
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
    }.get(NLIB, 'Unknown')
    info.evaluation=info.library # why do I need this?
    info.libraryVersion = "%d.%d.%d" % (NVER, LREL, NMOD)

    # Save the evaluation date
    try:
        info.Date = endfFileToGNDMisc.getENDFDate(_MTDatas[451][1][4][22:33])
    except Exception as e:
        if not kwargs.get('ignoreBadDate'):
            info.doRaise.append(str(e))
        import datetime
        info.Date = datetime.datetime.today().strftime("%Y-%m-%d")

    # Save general information
    info.author = _MTDatas[451][1][4][33:66]
    info.documentation = fudge.gnd.documentation.documentation('endfDoc', '\n'.join(_MTDatas[451][1][4:4 + NWD]))
    info.MAT  = _MAT
    info.LISO = LISO
    info.LIS  = LIS
    info.EMAX = EMAX
    info.levelIndex = LIS
    info.NSUB = NSUB
    info.LFI  = LFI
    info.STA  = STA
    info.targetTemperature=targetTemperature
    info.materialName = _MTDatas[451][1][4][:11].strip()
    info.targetMass = targetMass
    info.isNaturalTarget = isNaturalTarget
    info.LRP  = LRP
    info.NLIB = NLIB
    info.NMOD = NMOD
    info.NVER = NVER
    info.LREL = LREL
    info.NXC  = NXC
    info.NWD  = NWD
    info.LDRZ = LDRZ
    info.ITYPE= ITYPE

    # All evaluations need a particle database and we fill in a lot of that information from MF1MT451
    info.PoPsLabel = 'default'
    info.PoPs = databasePoPsModule.database('protare_internal', '1.0')

    # Let's tell info about the target
    info.targetZA = targetZA
    info.targetLevel = LIS
    info.level = targetExcitationEnergy
    info.addMassAWR(info.targetZA, targetMass)

    # If there is a projectile, tell info about it too
    if ITYPE in [4,5]:
        info.projectile=None
        info.projectileZA=None
    else:
        info.projectile = {0: 'g', 1: 'n', 11: 'e-', 1001: 'H1', 1002: 'H2', 1003: 'H3', 2003: 'He3', 2004: 'He4'}[IPART]
        info.projectileZA = projectileZA
        info.addMassAWR(projectileZA, projectileMass, asTarget=False)

    return info


def endfFileToGND( fileName, useFilesQAlways = True, singleMTOnly = None, evaluation = None,
        MTs2Skip = None, parseCrossSectionOnly = False,
        toStdOut = True, toStdErr = True, logFile = None, skipBadData = False, doCovariances = True,
        verboseWarnings = False, verbose = 1, reconstructResonances=True, **kwargs ) :

    logs = logFiles( toStdOut = toStdOut, toStdErr = toStdErr, logFile = logFile, defaultIsStderrWriting = False )
    header, MAT, MTDatas = endfFileToGNDMisc.parseENDFByMT_MF( fileName, logFile = logs )

    styleName = 'eval'
    reconstructedStyleName = 'recon'

    if MTs2Skip is None: MTs2Skip = []

    # Parse the ENDF documentation section
    info = readMF1MT451( MAT, MTDatas, styleName = styleName, logFile = logs, verboseWarnings = verboseWarnings, **kwargs )

    # OK, now decide what to do
    # First batch of ITYPE's are the special cases
    if( info.ITYPE == 2 ) : # Thermal scattering law data
        return( ENDF_ITYPE_2.ITYPE_2( MTDatas, info, verbose = verbose ) )
    elif( info.ITYPE == 4 ) : # Decay data, including SFPY
        return( ENDF_ITYPE_4.ITYPE_4( MTDatas, info, verbose = verbose ) )
    elif( info.ITYPE == 5 ) : # Decay data, including SFPY
        return( ENDF_ITYPE_5.ITYPE_5( MTDatas, info, verbose = verbose ) )
    elif( info.ITYPE == 6 )  : # Atomic relaxation data
        if( info.NSUB not in [ 6 ] ) : raise ValueError( 'For ITYPE = %d, invalid NSUB = %s' % ( info.ITYPE, info.NSUB ) )
        return( ENDF_ITYPE_6.ITYPE_6( info.targetZA / 1000, MTDatas, info, verbose = verbose ) )

    # Second batch of ITYPE's are general transport data
    if( info.LIS  != 0 ) : levelIndex = info.LIS
    if( info.ITYPE in [ 0, 1, 9 ] )  : # Nuclear transport data (g, n, p, d, t, h, a), including FPY
        target = toGNDMisc.getTypeNameGamma( info, info.targetZA, level = info.level, levelIndex = info.levelIndex )
        targetID = target.id
        if( ( info.STA != 0 ) ) :
            halflife = halflifeModule.string( info.PoPsLabel, halflifeModule.UNSTABLE, halflifeModule.baseUnit )
            target = info.PoPs[targetID]
            if( isinstance( target, nuclearLevelModule.particle ) ) :
                target = target.nucleus
            if( isinstance( target, isotopeModule.suite ) ) :
                target = target[0].nucleus
            if( len( target.halflife ) == 0 ) : target.halflife.add( halflife )
    elif( info.ITYPE == 3 ) :  # Atomic transport data (e, x)
        targetZ = info.targetZA / 1000
        targetID = chemicalElementModule.symbolFromZ[targetZ]
        targetName = chemicalElementModule.nameFromZ[targetZ]
        info.PoPs.add( chemicalElementModule.suite( targetID, targetZ, targetName ) )
    else :
        raise ValueError( "Unsupported ITYPE = %s" % info.ITYPE )

    # Add alias if needed for metastable target
    ZA2, MAT2 = endf_endlModule.ZAAndMATFromParticleName( targetID )
    MAT2 += info.LISO
    if( MAT2 != info.MAT ) :
        info.logs.write( "       WARNING: ENDF MAT = %s not as expected (i.e., %s).\n" % ( info.MAT, MAT2 ), stderrWriting = True )
    if( info.LISO != 0 ) :
        targetBaseName = targetID.split( '_' )[0]
        aliasName = alias.aliases.nuclearMetaStableName( targetBaseName, info.LISO )
        info.PoPs.add( PoPsAliasModule.metaStable( aliasName, targetID, info.LISO ) )

    # Compute the reconstructed and evaluated Styles
    evaluatedStyle = fudge.gnd.styles.evaluated(styleName, '',
                                                physicalQuantityModule.temperature(
                                                    PQUModule.pqu_float.surmiseSignificantDigits(info.targetTemperature),
                                                    'K'),
                                                info.library, info.libraryVersion, date=info.Date)
    if (evaluation is None): evaluation = "%s-%d.%d" % (info.library, info.NVER, info.LREL)
    info.reconstructedStyle = stylesModule.crossSectionReconstructed(reconstructedStyleName,
                                                                     derivedFrom=evaluatedStyle.label,
                                                                     date="2016-11-06")

    # Stuff for handling transport data
    info.printBadNK14 = True
    info.continuumSpectraFix = False
    info.ignoreMF10Fission = False
    options = ['printBadNK14', 'continuumSpectraFix', 'ignoreBadDate', 'ignoreMF10Fission']
    for option in kwargs:
        if (option not in options): raise DeprecationWarning('invalid deprecated option "%s"' % option)
        setattr(info, option, kwargs[option])
    info.evaluatedStyle = evaluatedStyle
    info.reconstructedAccuracy = 0.001
    info.MF12_LO2 = {}
    info.AWR_mode = None
    info.particleSpins = {}
    info.missingTwoBodyMasses = {}
    info.MF4ForNonNeutrons = []
    projectile = toGNDMisc.getTypeNameGamma(info, info.projectileZA)

    # Set up the reactionSuite
    reactionSuite = reactionSuiteModule.reactionSuite( projectile.id, targetID, evaluation,
            style = info.evaluatedStyle, documentation = info.documentation, MAT = MAT, PoPs = info.PoPs )
    MTDatas[451][1] = MTDatas[451][1][:4+info.NWD]
    info.reactionSuite = reactionSuite
    info.PoPs = reactionSuite.PoPs
    info.target = reactionSuite.target

    # Set up the covarianceSuite, if needed, and do other ITYPE-specific calculations
    covarianceSuite = None
    if( ( info.ITYPE == 0 ) or ( info.ITYPE == 9 ) ) :
        doRaise = info.NSUB not in { 0 : [ 0, 10, 10010, 10020, 10030, 20030, 20040 ], 9 : [ 19 ] }[info.ITYPE]
        if( doRaise ) : raise ValueError( 'For ITYPE = %d, invalid NSUB = %s' % ( info.ITYPE, info.NSUB ) )
        covarianceSuite = ENDF_ITYPE_0.ITYPE_0( MTDatas, info, reactionSuite, singleMTOnly, MTs2Skip,
                                                parseCrossSectionOnly, doCovariances, verbose,
                                                reconstructResonances=reconstructResonances )
    elif( info.ITYPE in [ 3, 1 ] ) :
        if( info.NSUB not in [ 3, 11, 113 ] ) :
            raise ValueError( 'For ITYPE = %d, invalid NSUB = %s' % ( info.ITYPE, info.NSUB ) )
        ENDF_ITYPE_3.ITYPE_3( MTDatas, info, reactionSuite, singleMTOnly, parseCrossSectionOnly, verbose = verbose )

    # Fill up the reactionSuite
    for reaction in reactionSuite.reactions : addUnspecifiedDistributions( info, reaction.outputChannel )
    for production in reactionSuite.productions : addUnspecifiedDistributions( info, production.outputChannel )
    PoPs = reactionSuite.PoPs

    lightMasses = { # atomic masses in amu, taken from AME2012 (http://www.nndc.bnl.gov/masses/mass.mas12)
        'n':  1.00866491585,
        'H1': 1.00782503223,
        'H2': 2.01410177812,
        'H3': 3.01604927791,
        'He3': 3.01602932008,
        'He4': 4.00260325413,
    }
    for product1ID in info.missingTwoBodyMasses :
        try:
            product1Mass = PoPs[product1ID].getMass( 'amu' )
        except Exception as exc: # FIXME currently getMass raises Exception, should probably be KeyError instead
            if product1ID in lightMasses:
                product1Mass = lightMasses[product1ID]
                PoPs[product1ID].mass.add(
                    massModule.double( info.PoPsLabel, product1Mass, quantityModule.stringToPhysicalUnit('amu')) )
            else:
                raise exc

        massPTmP1 = PoPs[reactionSuite.projectile].getMass( 'amu' ) + PoPs[reactionSuite.target].getMass( 'amu' ) - \
                    product1Mass
        summedQM = 0
        for product2ID, QM in info.missingTwoBodyMasses[product1ID] : summedQM += QM
        mass = massPTmP1 - PQUModule.PQU( summedQM / len( info.missingTwoBodyMasses[product1ID] ), 'eV/c**2' ).getValueAs( 'amu' )
        mass = massModule.double( info.PoPsLabel, mass, quantityModule.stringToPhysicalUnit( 'amu' ) )
        residual = reactionSuite.PoPs[info.missingTwoBodyMasses[product1ID][0][0]]
        if( len( residual.mass ) == 0 ) : residual.mass.add( mass )

    # Did we have trouble?
    if( len( info.MF4ForNonNeutrons ) > 0 ) :
        print '    WARNING: MF=4 data for non-neutron product: MTs = %s' % sorted( info.MF4ForNonNeutrons )
    if( len( info.doRaise ) > 0 and not skipBadData ) :
        info.logs.write( '\nRaising due to following errors:\n' )
        for err in info.doRaise : info.logs.write( err + '\n' )
        raise Exception( 'len( info.doRaise ) > 0' )

    return( { 'reactionSuite' : reactionSuite, 'covarianceSuite' : covarianceSuite, 'errors' : info.doRaise, 'info':info } )

def addUnspecifiedDistributions( info, outputChannel ) :

    if( outputChannel is None ) : return
    for product in outputChannel :
        if( len( product.distribution ) == 0 ) :
            product.distribution.add( unspecifiedModule.form( info.style, productFrame = standardsModule.frames.labToken ) )
        addUnspecifiedDistributions( info, product.outputChannel )

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
