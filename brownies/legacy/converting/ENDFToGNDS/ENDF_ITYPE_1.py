# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
For translating ENDF ITYPE=1 data (induced fission product yields)
"""

from xData import xDataArray as arrayModule

from PoPs import alias as aliasModule
from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.quantities import halflife as halflifeModule
from PoPs.groups import misc as chemicalElementMiscPoPsModule
from PoPs.fissionFragmentData import time as timeModule

from fudge.channelData.fissionFragmentData import fissionFragmentData as fissionFragmentDataModule
from fudge.channelData.fissionFragmentData import productYield as productYieldModule
from fudge.channelData.fissionFragmentData import duration as durationModule, incidentEnergy as incidentEnergyModule
from fudge.channelData.fissionFragmentData import yields as yieldsModule

from brownies.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc as endfFileToGNDSMiscModule

from .. import toGNDSMisc as toGNDSMiscModule

energyUnit = 'eV'

def parseMF454_459( info, fissionFragmentData, MT, MTDatas ) :

    MFData = MTDatas[MT]
    MF8Data = MFData[8]

    ZA, AWR, LE_plus1 = endfFileToGNDSMiscModule.sixFunkyFloatStringsToIntsAndFloats( MF8Data[0], intIndices = [ 0, 2 ], logFile = info.logs )[0:3]
    line = 1

    if( MT == 454 ) :
        duration = durationModule.duration( 'initial' )
        duration.time.add( timeModule.double( 'initial', 0.0, 's' ) )
    else :
        duration = durationModule.duration( 'unspecified' )
        duration.time.add( timeModule.string( 'unspecified', 'unspecified', 's' ) )

    for iLE in range( LE_plus1 ) :
        line, ENDFList = endfFileToGNDSMiscModule.getList( line, MF8Data, logFile = info.logs )
        NFP = ENDFList['N2']
        data = ENDFList['data']
        interpolation = ENDFList['L1']                      # Currently ignored.

        fissionProducts = []
        yields = []
        uncertainties = []
        for iFP in range( NFP ):
            ZAFP = int( data[iFP*4] )
            FPS =  int( data[iFP*4+1] )
            FPID = chemicalElementMiscPoPsModule.idFromZA( ZAFP )
            if( FPS > 0 ) : FPID = aliasModule.metaStable.metaStableNameFromNuclearLevelNameAndMetaStableIndex( FPID, FPS )
            if( FPID in fissionProducts ) : raise Exception( 'FPID = %s already in fissionProducts' % FPID )
            fissionProducts.append( FPID )

            yields.append( data[iFP*4+2] )
            uncertainties.append( data[iFP*4+3] )

        if( ( MT == 454 ) and ( iLE == 0 ) ) :
            productYield = productYieldModule.productYield( info.style )
            fissionFragmentData.productYields.add( productYield )
        else :
            productYield = fissionFragmentData.productYields[0]

        incidentEnergy = incidentEnergyModule.incidentEnergy( str( iLE ), fissionProducts )
        incidentEnergy.energy.add( incidentEnergyModule.double( str( iLE ), ENDFList[ 'C1' ], energyUnit ) )

        incidentEnergy.values = yieldsModule.values( yields )

        diagonal = arrayModule.diagonal( shape = ( len( uncertainties ), len( uncertainties  ) ), data = uncertainties )
        covariance = yieldsModule.covariance( diagonal )
        incidentEnergy.uncertainty = yieldsModule.uncertainty( covariance )

        duration.incidentEnergies.add( incidentEnergy )

    productYield.durations.add( duration )

def ITYPE_1( MTDatas, info, verbose = 0 ) :
    """Parses ENDF NFY data."""

    errors = []

    info.PoPs.documentation.endfCompatible.body = '\n'.join( MTDatas[451][1][4:-info.NXC] )         # Add ENDF documentation

    parentAtom = toGNDSMiscModule.getPoPsParticle( info, info.targetZA, name = None, levelIndex = info.levelIndex, level = info.level, levelUnit = energyUnit )
    if( len( parentAtom.mass ) == 0 ) :
        parentAtom.mass.add( massModule.double( info.PoPsLabel, info.massTracker.neutronMass * info.targetMass, quantityModule.stringToPhysicalUnit( 'amu' ) ) )
    if( len( parentAtom.nucleus.halflife ) == 0 ) :
        parentAtom.nucleus.halflife.add( halflifeModule.string( info.PoPsLabel, halflifeModule.UNSTABLE, halflifeModule.baseUnit ) )

    fissionFragmentData = fissionFragmentDataModule.fissionFragmentData( )

    if( 454 in MTDatas ) : parseMF454_459( info, fissionFragmentData, 454, MTDatas )
    if( 459 in MTDatas ) : parseMF454_459( info, fissionFragmentData, 459, MTDatas )

    return { 'fissionFragmentData' : fissionFragmentData, 'PoPs' : info.PoPs, 'errors' : errors, 'info' : info }
