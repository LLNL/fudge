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

from xData import standards as standardsModule

from PoPs import database as databasePoPsModule
from PoPs import alias as PoPsAliasModule
from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import halflife as halflifeModule
from PoPs.quantities import mass as massModule
from PoPs.families import nuclide as nuclideModule
from PoPs.groups import isotope as isotopeModule
from PoPs.groups import chemicalElement as chemicalElementModule
from PoPs.groups import misc as chemicalElementMiscModule

from fudge.gnds import institution as institutionModule
from fudge.gnds import physicalQuantity as physicalQuantityModule
from fudge.gnds import reactionSuite as reactionSuiteModule
from fudge.gnds import styles as stylesModule
from fudge.gnds.productData.distributions import unspecified as unspecifiedModule

from .ENDFToGNDS import endfFileToGNDSMisc, ENDF_ITYPE_0
from .ENDFToGNDS import ENDF_ITYPE_2, ENDF_ITYPE_3, ENDF_ITYPE_4, ENDF_ITYPE_6
from .ENDFToGNDS import ENDF_ITYPE_1 as ENDF_ITYPE_1Module
from .ENDFToGNDS import ENDF_ITYPE_5 as ENDF_ITYPE_5Module
from . import endf_endl as endf_endlModule
from . import toGNDSMisc

from site_packages.legacy.toENDF6 import ENDFconversionFlags as ENDFconversionFlagsModule

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
    info = toGNDSMisc.infos(styleName)
    info.doRaise = []
    info.logs = logFile
    info.verboseWarnings = verboseWarnings
    info.ENDFconversionFlags = ENDFconversionFlagsModule.ENDFconversionFlags()

    # Line #1
    targetZA, targetMass, LRP, LFI, NLIB, NMOD = \
        endfFileToGNDSMisc.sixFunkyFloatStringsToFloats( _MTDatas[451][1][0], logFile=logFile )
    targetZA = int(targetZA)  # Target's ZA
    LRP = int(LRP)            # Resonance parameter data info
    LFI = int(LFI)            # Is fission present
    NLIB = int(NLIB)          # What library (e.g., 0 = ENDF/B
    NMOD = int(NMOD)          # Version modification flag
    isNaturalTarget = (targetZA % 1000) == 0

    # Line #2
    targetExcitationEnergy, STA, LIS, LISO, dummy, NFOR = \
        endfFileToGNDSMisc.sixFunkyFloatStringsToFloats( _MTDatas[451][1][1], logFile=logFiles )
    STA = int(STA)            # Is nucleus unstable
    LIS = int(LIS)            # Excitation number
    LISO = int(LISO)          # Isomeric state number
    if( ( targetExcitationEnergy != 0 ) or ( LIS != 0 ) or ( LISO != 0 ) ) :
        if( ( targetExcitationEnergy == 0 ) or ( LIS == 0 ) or ( LISO == 0 ) ) :
            print( "    WARNING: target excitation energy = %s, LIS = %s and LISO = %s are not consistent." % ( targetExcitationEnergy, LIS, LISO ) )
    NFOR = int(NFOR)          # Must be 6 for ENDF-6 format
    if (NFOR != 6):
        print ("    WARNING: endfFileToGNDS only supports ENDF-6 format. This file has unsupported NFOR=%d" % NFOR)

    # Line #3
    projectileMass, EMAX, LREL, dummy, NSUB, NVER = \
        endfFileToGNDSMisc.sixFunkyFloatStringsToFloats( _MTDatas[451][1][2], logFile=logFiles )
    NSUB = int(NSUB)  # 10 * ZA + iType for projectile
    NVER = int(NVER)  # Evaluation version number
    LREL = int(LREL)  # Evaluation sub-version number
    IPART, ITYPE = NSUB / 10, NSUB % 10
    projectileZA = IPART
    if (projectileZA == 11): projectileZA = 9

    # Line #4
    targetTemperature, dummy, LDRZ, dummy, NWD, NXC = \
        endfFileToGNDSMisc.sixFunkyFloatStringsToFloats( _MTDatas[451][1][3], logFile=logFiles )
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
        info.Date = endfFileToGNDSMisc.getENDFDate(_MTDatas[451][1][4][22:33])
    except Exception as e:
        if not kwargs.get('ignoreBadDate'):
            info.doRaise.append(str(e))
        import datetime
        info.Date = datetime.datetime.today().strftime("%Y-%m-%d")

    # Save general information
    info.author = _MTDatas[451][1][4][33:66]
    info.documentation = fudge.gnds.documentation.documentation('endfDoc', '\n'.join(_MTDatas[451][1][4:4 + NWD]))
    info.MAT  = _MAT
    info.LISO = LISO
    info.LIS  = LIS
    info.EMAX = EMAX
    info.levelIndex = LIS
    info.NSUB = NSUB
    info.LFI  = LFI
    info.STA  = STA
    info.targetTemperature = targetTemperature
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
    info.PoPsLabel = styleName
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


def endfFileToGNDS( fileName, useFilesQAlways = True, singleMTOnly = None, evaluation = None,
                    MTs2Skip = None, parseCrossSectionOnly = False,
                    toStdOut = True, toStdErr = True, logFile = None, skipBadData = False, doCovariances = True,
                    verboseWarnings = False, verbose = 1, reconstructResonances=True, **kwargs ) :

    logs = logFiles( toStdOut = toStdOut, toStdErr = toStdErr, logFile = logFile, defaultIsStderrWriting = False )
    header, MAT, MTDatas = endfFileToGNDSMisc.parseENDFByMT_MF( fileName, logFile = logs )

    styleName = 'eval'
    reconstructedStyleName = 'recon'

    if MTs2Skip is None: MTs2Skip = []

    # Parse the ENDF documentation section
    info = readMF1MT451( MAT, MTDatas, styleName = styleName, logFile = logs, verboseWarnings = verboseWarnings, **kwargs )

    if( info.ITYPE == 1 ) :                                                         # NFY
        return( ENDF_ITYPE_1Module.ITYPE_1( MTDatas, info, verbose = verbose ) )
    elif( info.ITYPE == 2 ) :                                                       # Thermal scattering law data
        return( ENDF_ITYPE_2.ITYPE_2( MTDatas, info, verbose = verbose ) )
    elif( info.ITYPE == 4 ) :                                                       # Decay data
        return( ENDF_ITYPE_4.ITYPE_4( MTDatas, info, verbose = verbose ) )
    elif( info.ITYPE == 5 ) :                                                       # SFY
        return( ENDF_ITYPE_5Module.ITYPE_5( MTDatas, info, verbose = verbose ) )
    elif( info.ITYPE == 6 )  :                                                      # Atomic relaxation data
        if( info.NSUB not in [ 6 ] ) : raise ValueError( 'For ITYPE = %d, invalid NSUB = %s' % ( info.ITYPE, info.NSUB ) )
        return( ENDF_ITYPE_6.ITYPE_6( info.targetZA / 1000, MTDatas, info, verbose = verbose ) )

    # Second batch of ITYPE's are general transport data
    if( info.LIS  != 0 ) : levelIndex = info.LIS
    if( info.ITYPE in [ 0, 1, 9 ] )  : # Nuclear transport data (g, n, p, d, t, h, a), including FPY
        target = toGNDSMisc.getTypeNameGamma( info, info.targetZA, level = info.level, levelIndex = info.levelIndex )
        targetID = target.id
        if( ( info.STA != 0 ) ) :
            halflife = halflifeModule.string( info.PoPsLabel, halflifeModule.UNSTABLE, halflifeModule.baseUnit )
            target = info.PoPs[targetID]
            if( isinstance( target, nuclideModule.particle ) ) : target = target.nucleus
            if( len( target.halflife ) == 0 ) : target.halflife.add( halflife )
    elif( info.ITYPE == 3 ) :  # Atomic transport data (e, x)
        targetZ = info.targetZA / 1000
        targetID = chemicalElementMiscModule.symbolFromZ[targetZ]
        targetName = chemicalElementMiscModule.nameFromZ[targetZ]
        info.PoPs.add( chemicalElementModule.chemicalElement( targetID, targetZ, targetName ) )
    else :
        raise ValueError( "Unsupported ITYPE = %s" % info.ITYPE )

    # Add alias if needed for metastable target
    ZA2, MAT2 = endf_endlModule.ZAAndMATFromParticleName( targetID )
    MAT2 += info.LISO
    if( MAT2 != info.MAT ) :
        info.logs.write( "       WARNING: ENDF MAT = %s not as expected (i.e., %s).\n" % ( info.MAT, MAT2 ), stderrWriting = True )
    info.targetID = targetID
    if( info.LISO != 0 ) :
        aliasName = PoPsAliasModule.metaStable.metaStableNameFromNuclearLevelNameAndMetaStableIndex( targetID, info.LISO )
        info.PoPs.add( PoPsAliasModule.metaStable( aliasName, targetID, info.LISO ) )
        info.targetID = aliasName

    # Compute the reconstructed and evaluated Styles
    projectileDomain = stylesModule.projectileEnergyDomain( 1e-5, 2e+7, 'eV' )  # will be overwritten after parsing all reactions
    evaluatedStyle = fudge.gnds.styles.evaluated(styleName, '',
                                                physicalQuantityModule.temperature(
                                                    PQUModule.pqu_float.surmiseSignificantDigits(info.targetTemperature),
                                                    'K'),
                                                projectileDomain,
                                                info.library, info.libraryVersion, date=info.Date)
    if (evaluation is None): evaluation = "%s-%d.%d" % (info.library, info.NVER, info.LREL)
    info.reconstructedStyle = stylesModule.crossSectionReconstructed(reconstructedStyleName,
                                                                     derivedFrom=evaluatedStyle.label,
                                                                     date="2016-11-06")

    # Stuff for handling transport data
    info.printBadNK14 = True
    info.continuumSpectraFix = False
    info.acceptBadMF10FissionZAP = False
    options = ['printBadNK14', 'continuumSpectraFix', 'ignoreBadDate', 'acceptBadMF10FissionZAP']
    for option in kwargs:
        if (option not in options): raise DeprecationWarning('invalid or deprecated option "%s"' % option)
        setattr(info, option, kwargs[option])
    info.evaluatedStyle = evaluatedStyle
    info.reconstructedAccuracy = 0.001
    info.MF12_LO2 = {}
    info.AWR_mode = None
    info.particleSpins = {}
    info.missingTwoBodyMasses = {}
    info.MF4ForNonNeutrons = []
    projectile = toGNDSMisc.getTypeNameGamma(info, info.projectileZA)

    # Set up the reactionSuite
    reactionSuite = reactionSuiteModule.reactionSuite( projectile.id, info.targetID, evaluation,
            style = info.evaluatedStyle, documentation = info.documentation, MAT = MAT, PoPs = info.PoPs )
    MTDatas[451][1] = MTDatas[451][1][:4+info.NWD]
    info.reactionSuite = reactionSuite
    info.PoPs = reactionSuite.PoPs
    info.target = reactionSuite.target

    # Set up the covarianceSuite, if needed, and do other ITYPE-specific calculations
    covarianceSuite = None
    if( info.ITYPE in [ 0, 9 ] ) :
        doRaise = info.NSUB not in { 0 : [ 0, 10, 10010, 10020, 10030, 20030, 20040 ], 9 : [ 19 ] }[info.ITYPE]
        if( doRaise ) : raise ValueError( 'For ITYPE = %d, invalid NSUB = %s' % ( info.ITYPE, info.NSUB ) )
        covarianceSuite = ENDF_ITYPE_0.ITYPE_0( MTDatas, info, reactionSuite, singleMTOnly, MTs2Skip,
                                                parseCrossSectionOnly, doCovariances, verbose,
                                                reconstructResonances=reconstructResonances )
    elif( info.ITYPE in [ 3 ] ) :
        if( info.NSUB not in [ 3, 11, 113 ] ) :
            raise ValueError( 'For ITYPE = %d, invalid NSUB = %s' % ( info.ITYPE, info.NSUB ) )
        ENDF_ITYPE_3.ITYPE_3( MTDatas, info, reactionSuite, singleMTOnly, parseCrossSectionOnly, verbose = verbose )

    # Fix incident energy domain limits:
    reactionList = reactionSuite.reactions
    if len(reactionList) == 0:
        reactionList = reactionSuite
    nonThresholdList = [reac for reac in reactionSuite.reactions if reac.outputChannel.Q.getConstant() >= 0]

    if len(nonThresholdList) > 0:
        domainMin = max( [reac.crossSection.domainMin for reac in nonThresholdList] )
    else:
        domainMin = min( [reac.crossSection.domainMin for reac in reactionList] )
    domainMax = min( [reac.crossSection.domainMax for reac in reactionList] )

    info.evaluatedStyle.projectileEnergyDomain.min = domainMin
    info.evaluatedStyle.projectileEnergyDomain.max = domainMax

    # Fill up the reactionSuite
    for reaction in reactionSuite.reactions : addUnspecifiedDistributions( info, reaction.outputChannel )
    for production in reactionSuite.productions : addUnspecifiedDistributions( info, production.outputChannel )

    if len(info.ENDFconversionFlags.flags) > 0:
        flagsToPurge = []
        for conversion in info.ENDFconversionFlags.flags:
            if conversion.link.getRootAncestor() is covarianceSuite:
                conversion.root = '$covariances'    # FIXME need better way of handling this
            elif conversion.link.getRootAncestor() is not reactionSuite:
                flagsToPurge.append( conversion )   # flag doesn't point to any GNDS node
        for flag in flagsToPurge:
            info.ENDFconversionFlags.flags.remove( flag )
        info.ENDFconversionFlags.flags.sort(key = lambda val: val['flags'])
        LLNLdata = institutionModule.institution("LLNL")
        LLNLdata.data.append(info.ENDFconversionFlags)
        reactionSuite.applicationData.add( LLNLdata )

    PoPs = reactionSuite.PoPs

    lightMasses = { # atomic masses in amu, taken from AME2012 (http://www.nndc.bnl.gov/masses/mass.mas12)
        'n':  1.00866491585,
        'H1': 1.00782503223,
        'H2': 2.01410177812,
        'H3': 3.01604927791,
        'He3': 3.01602932008,
        'He4': 4.00260325413,
    }
    for pid in lightMasses:
        if pid in info.PoPs:
            if len(info.PoPs[pid].mass) == 0:
                info.PoPs[pid].mass.add(
                    massModule.double( info.PoPsLabel, lightMasses[pid], quantityModule.stringToPhysicalUnit('amu') )
                )

    for product1ID in info.missingTwoBodyMasses :
        if len(info.PoPs[product1ID].mass) != 0:    # mass was added to PoPs, likely in previous block
            product1Mass = info.PoPs[product1ID].getMass("amu")
        else:
            raise ValueError("Missing mass for 2-body reaction product %s" % product1ID)

        massPTmP1 = PoPs[reactionSuite.projectile].getMass( 'amu' ) + PoPs[targetID].getMass( 'amu' ) - \
                    product1Mass
        summedQM = 0
        for product2ID, QM in info.missingTwoBodyMasses[product1ID] : summedQM += QM
        mass = massPTmP1 - PQUModule.PQU( summedQM / len( info.missingTwoBodyMasses[product1ID] ), 'eV/c**2' ).getValueAs( 'amu' )
        mass = massModule.double( info.PoPsLabel, mass, quantityModule.stringToPhysicalUnit( 'amu' ) )
        residual = reactionSuite.PoPs[info.missingTwoBodyMasses[product1ID][0][0]]
        if( len( residual.mass ) == 0 ) : residual.mass.add( mass )

    # Did we have trouble?
    if( len( info.MF4ForNonNeutrons ) > 0 ) :
        print('    WARNING: MF=4 data for non-neutron product: MTs = %s' % sorted( info.MF4ForNonNeutrons ))
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
            info.ENDFconversionFlags.add( product, 'implicitProduct' )
        addUnspecifiedDistributions( info, product.outputChannel )

if( __name__ == '__main__' ) :

    rce = endfFileToGNDS( fileName = sys.argv[1])
    x, c = rce['reactionSuite'], rce['covarianceSuite']
    f = open( 'test.xml', 'w' )
    f.write( '\n'.join( x.toXMLList( ) + [ '' ] ) )
    f.close( )
    if( c is not None ) : # covariances
        f = open( 'test-covar.xml', 'w' )
        f.write( '\n'.join( c.toXMLList( ) + [ '' ] ) )
        f.close()
