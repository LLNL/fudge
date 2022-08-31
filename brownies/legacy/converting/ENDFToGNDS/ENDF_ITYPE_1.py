# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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
from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule

from PoPs.fissionFragmentData import nuclides as nuclidesModule
from PoPs.fissionFragmentData import time as timeModule
from PoPs.fissionFragmentData import yields as yieldsModule

from fudge.outputChannelData.fissionFragmentData import fissionFragmentData as fissionFragmentDataModule
from fudge.outputChannelData.fissionFragmentData import elapsedTime as elapsedTimeModule
from fudge.outputChannelData.fissionFragmentData import productYield as productYieldModule
from fudge.outputChannelData.fissionFragmentData import incidentEnergy as incidentEnergyModule

from brownies.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc as endfFileToGNDSMiscModule

from .. import toGNDSMisc as toGNDSMiscModule

energyUnit = 'eV'

def parseMF454_459( info, fissionFragmentData, MT, MTDatas ) :

    MFData = MTDatas[MT]
    MF8Data = MFData[8]

    ZA, AWR, LE_plus1 = endfFileToGNDSMiscModule.sixFunkyFloatStringsToIntsAndFloats( MF8Data[0], intIndices = [ 0, 2 ], logFile = info.logs )[0:3]
    line = 1

    if( MT == 454 ) :
        elapsedTime = elapsedTimeModule.ElapsedTime( 'initial' )
        elapsedTime.time.add( timeModule.Double( 'initial', 0.0, 's' ) )
    else :
        elapsedTime = elapsedTimeModule.ElapsedTime( 'unspecified' )
        elapsedTime.time.add( timeModule.String( 'unspecified', 'unspecified', 's' ) )

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
            if( FPS > 0 ) : FPID = aliasModule.MetaStable.metaStableNameFromNuclearLevelNameAndMetaStableIndex(FPID, FPS)
            if( FPID in fissionProducts ) : raise Exception( 'FPID = %s already in fissionProducts' % FPID )
            fissionProducts.append( FPID )

            yields.append( data[iFP*4+2] )
            uncertainties.append( data[iFP*4+3] )

        if( ( MT == 454 ) and ( iLE == 0 ) ) :
            productYield = productYieldModule.ProductYield( info.style )
            fissionFragmentData.productYields.add( productYield )
        else :
            productYield = fissionFragmentData.productYields[0]

        incidentEnergy = incidentEnergyModule.IncidentEnergy( str( iLE ) )
        incidentEnergy.energy.add( incidentEnergyModule.Double( str( iLE ), ENDFList[ 'C1' ], energyUnit ) )

        incidentEnergy.yields.nuclides = nuclidesModule.Nuclides(fissionProducts)
        incidentEnergy.yields.values = yieldsModule.Values( yields )

        variances = [u**2 for u in uncertainties]
        diagonal = arrayModule.Diagonal( shape = ( len( variances ), len( variances  ) ), data = variances )
        covariance = yieldsModule.Covariance( diagonal )
        incidentEnergy.yields.uncertainty = yieldsModule.Uncertainty( covariance )

        elapsedTime.incidentEnergies.add( incidentEnergy )

    productYield.elapsedTimes.add( elapsedTime )

def ITYPE_1( MTDatas, info, verbose = 0 ) :
    """Parses ENDF NFY data."""

    errors = []

    info.PoPs.documentation.endfCompatible.body = '\n'.join( MTDatas[451][1][4:-info.NXC] )         # Add ENDF documentation
    endfFileToGNDSMiscModule.completeDocumentation(info, info.PoPs.documentation)

    parentAtom = toGNDSMiscModule.getPoPsParticle( info, info.targetZA, name = None, levelIndex = info.levelIndex, level = info.level, levelUnit = energyUnit )
    if( len( parentAtom.mass ) == 0 ) :
        parentAtom.mass.add( massModule.Double( info.PoPsLabel, info.massTracker.neutronMass * info.targetMass, quantityModule.stringToPhysicalUnit( 'amu' ) ) )
    if( len( parentAtom.nucleus.halflife ) == 0 ) :
        parentAtom.nucleus.halflife.add( halflifeModule.String( info.PoPsLabel, halflifeModule.UNSTABLE, halflifeModule.baseUnit ) )

    fissionFragmentData = fissionFragmentDataModule.FissionFragmentData( )

    if( 454 in MTDatas ) : parseMF454_459( info, fissionFragmentData, 454, MTDatas )
    if( 459 in MTDatas ) : parseMF454_459( info, fissionFragmentData, 459, MTDatas )

    # check if all elapsed times / incident energies use same list of nuclides
    nuclides = fissionFragmentData.productYields[info.style].elapsedTimes[0].incidentEnergies[0].yields.nuclides
    allIdentical = True
    for et in fissionFragmentData.productYields[info.style].elapsedTimes:
        for ie in et.incidentEnergies:
            if ie.yields.nuclides.data != nuclides.data:
                allIdentical = False
                break

    if allIdentical:
        # move list of nuclides up to productYield
        oldProductYield = fissionFragmentData.productYields.pop(info.style)
        newProductYield = productYieldModule.ProductYield(info.style, nuclides=nuclides)
        for et in oldProductYield.elapsedTimes:
            for ie in et.incidentEnergies:
                ie.yields.nuclides = yieldsModule.NuclidesLink(newProductYield.nuclides, relative=True)
            newProductYield.elapsedTimes.add(et)

        fissionFragmentData.productYields.add(newProductYield)

    return { 'fissionFragmentData' : fissionFragmentData, 'PoPs' : info.PoPs, 'errors' : errors, 'info' : info }
