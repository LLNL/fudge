# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
For translating ENDF ITYPE=5 data (spontaneous fission product yields)
"""

from xData import xDataArray as arrayModule

from PoPs import alias as aliasModule
from PoPs import styles as stylesPoPsModule
from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.quantities import halflife as halflifeModule
from PoPs.groups import misc as chemicalElementMiscPoPsModule
from PoPs.fissionFragmentData import productYield as productYieldModule
from PoPs.fissionFragmentData import duration as durationModule
from PoPs.fissionFragmentData import time as timeModule
from PoPs.fissionFragmentData import yields as yieldsModule

from brownies.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc as endfFileToGNDSMiscModule

from .. import toGNDSMisc as toGNDSMiscModule

energyUnit = 'eV'

def parseMF454_459( info, fissionFragmentData, MT, MTDatas ) :

    MFData = MTDatas[MT]
    MF8Data = MFData[8]

    ZA, AWR, LE_plus1 = endfFileToGNDSMiscModule.sixFunkyFloatStringsToIntsAndFloats( MF8Data[0], intIndices = [ 0, 2 ], logFile = info.logs )[0:3]
    if( LE_plus1 != 1 ) : raise Exception( 'This is not spontaneous fission yield data.' )

    line = 1
    fissionProducts = []
    line, ENDFList = endfFileToGNDSMiscModule.getList( line, MF8Data, logFile = info.logs )
    NFP = ENDFList ['N2']
    data = ENDFList['data']
    yields = []
    uncertainties = []
    for iFP in range( NFP ):
        ZAFP = int( data[iFP*4] )
        FPS =  int( data[iFP*4+1] )
        FPID = chemicalElementMiscPoPsModule.idFromZA( ZAFP )
        if( FPS > 0 ) : FPID = aliasModule.metaStable.metaStableNameFromNuclearLevelNameAndMetaStableIndex( FPID, FPS )
        if( FPID in fissionProducts ) : raise Exception( 'FPID = "%s" already in fissionProducts.' % FPID )
        fissionProducts.append( FPID )

        yields.append( data[iFP*4+2] )
        uncertainties.append( data[iFP*4+3] )

    if( len( fissionFragmentData.productYields ) == 0 ) :
        fissionFragmentData.productYields.add( productYieldModule.productYield( info.style, fissionProducts ) )
    else :
        if( fissionFragmentData.productYields[0].nuclides != fissionProducts ) :
            raise Exception( 'Fission products differ between independent and cumulative.' )

    duration = durationModule.duration( )
    if( MT == 454 ) :
        duration.time.add( timeModule.double( info.style, 0.0, 's' ) )
    else :
        duration.time.add( timeModule.string( info.style, 'unspecified', 's' ) )

    duration.yields.values = yieldsModule.values( yields )

    diagonal = arrayModule.diagonal( shape = ( len( uncertainties ), len( uncertainties  ) ), data = uncertainties )
    covariance = yieldsModule.covariance( diagonal )
    uncertainty = yieldsModule.uncertainty( covariance )
    duration.yields.uncertainty = uncertainty

    fissionFragmentData.productYields[0].durations.add( duration )

def ITYPE_5( MTDatas, info, verbose = 0 ) :
    """Parses ENDF SFY data."""

    errors = []

    evalStyle = stylesPoPsModule.evaluated( 'eval', '', info.library, info.libraryVersion, info.Date )
    evalStyle.documentation.endfCompatible.body = '\n'.join( MTDatas[451][1][4:-info.NXC] )     # Add ENDF documentation
    info.PoPs.styles.add(evalStyle)

    parentAtom = toGNDSMiscModule.getPoPsParticle( info, info.targetZA, name = None, levelIndex = info.levelIndex, level = info.level, levelUnit = energyUnit )
    if( len( parentAtom.mass ) == 0 ) :
        parentAtom.mass.add( massModule.double( info.PoPsLabel, info.massTracker.neutronMass * info.targetMass, quantityModule.stringToPhysicalUnit( 'amu' ) ) )
    if( len( parentAtom.nucleus.halflife ) == 0 ) :
        parentAtom.nucleus.halflife.add( halflifeModule.string( info.PoPsLabel, halflifeModule.UNSTABLE, halflifeModule.baseUnit ) )

    fissionFragmentData = parentAtom.fissionFragmentData

    if( 454 in MTDatas ) : parseMF454_459( info, fissionFragmentData, 454, MTDatas )
    if( 459 in MTDatas ) : parseMF454_459( info, fissionFragmentData, 459, MTDatas )

    return { 'PoPs' : info.PoPs, 'errors' : errors, 'info' : info }
