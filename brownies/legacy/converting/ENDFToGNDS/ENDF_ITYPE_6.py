# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
For translating ENDF ITYPE=6 data  (atomic relaxation sub-library)
"""


from PoPs import IDs as IDsPoPsModule
from PoPs import styles as stylesPoPsModule
from PoPs import database as databasePoPsModule
from PoPs.groups import misc as chemicalElementMiscModule
from PoPs.groups import chemicalElement as chemicalElementModule
from PoPs.atomic import atomic as atomicModule
from PoPs.atomic import configuration as configurationModule
from PoPs.decays import decayData as decayDataModule
from PoPs.decays import product as productModule
from PoPs.decays import probability as probabilityModule
from PoPs.quantities import bindingEnergy as bindingEnergyModule

from . import endfFileToGNDSMisc
from .ENDF_ITYPE_3_6_Misc import MT_AtomicConfigurations

def ITYPE_6( Z, MTDatas, info, verbose = 0 ) :
    """Parses ENDF atomic relaxation data."""

    MT = 533
    errors = []
    elementSymbol = chemicalElementMiscModule.symbolFromZ[Z]

    MTs = sorted( MTDatas.keys( ) )
    if( 451 in MTs ) : MTs.remove( 451 )
    if( MTs != [ MT ] ) : raise Exception( 'Only allowed MT is %s (and 451): %s' % ( MT, MTs ) )
    warningList = []

    MFData = MTDatas[MT]
    info.logs.write( '    %3d %s' % ( MT, sorted( MFData.keys( ) ) ) )

    if( list( MFData.keys( ) ) != [ 28 ] ) :
        raise Exception( 'For MT %s data, only allowed MF is 28: %s' % ( MT, list( MFData.keys( ) ) ) )
    MF28 = MFData[28]

    PoPs = databasePoPsModule.database( "ENDF atomic relaxation for Z = %s" % Z, '1.0', formatVersion = info.formatVersion )
    evalStyle = stylesPoPsModule.evaluated( 'eval', '', info.library, info.libraryVersion, info.Date )
    evalStyle.documentation.endfCompatible.body = info.documentation

    PoPs.styles.add( evalStyle )
    chemicalElement = chemicalElementModule.chemicalElement( elementSymbol, Z, chemicalElementMiscModule.nameFromZ[Z] )
    PoPs.add( chemicalElement )

    atomicData = atomicModule.atomic( )
    chemicalElement.atomicData = atomicData

    offset = 0
    ZA, AWP, dummy, dummy, NSS, dummy = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats(MF28[offset], intIndices = [4], logFile = info.logs)
    offset += 1
    for subshell in range( NSS ) :
        SUBI, dummy, dummy, dummy, dummy, NTR = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats(MF28[offset], intIndices = [0, 5], logFile = info.logs)
        offset += 1
        EBI, ELN, dummy, dummy, dummy, dummy = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats(MF28[offset], intIndices = [ ], logFile = info.logs)
        offset += 1

        configuration = configurationModule.configuration( MT_AtomicConfigurations[534+SUBI-1], ELN )
        atomicData.configurations.add( configuration )
        configuration.bindingEnergy.add( bindingEnergyModule.double( 'eval', EBI, 'eV' ) )

        decayData = configuration.decayData

        probabilitySum = 0.
        info.logs.write( ' : %s' % MT_AtomicConfigurations[534+SUBI-1] )
        for i1 in range( NTR ) :
            SUBJ, SUBK, ETR, FTF, dummy, dummy = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats(MF28[offset], intIndices = [0, 1], logFile = info.logs)
            offset += 1

            decayMode = decayDataModule.decayMode( str( i1 ), 'atomicRelaxation' )
            decayMode.probability.add( probabilityModule.double( "BR", FTF, '' ) )
            decay = decayDataModule.decay( str( 0 ), decayDataModule.decayModesParticle )

            if( SUBK == 0 ) :
                decay.products.add( productModule.product( IDsPoPsModule.photon, IDsPoPsModule.photon ) )
                elementID = '%s{%s}' % ( elementSymbol, MT_AtomicConfigurations[534+SUBJ-1] )
            else :
                decay.products.add( productModule.product( IDsPoPsModule.electron, IDsPoPsModule.electron ) )
                elementID = '%s{%s,%s}' % ( elementSymbol, MT_AtomicConfigurations[534+SUBJ-1], MT_AtomicConfigurations[534+SUBK-1] )
            decay.products.add( productModule.product( elementID, elementID ) )
            probabilitySum += FTF
            decayMode.decayPath.add( decay )
            decayData.decayModes.add( decayMode )
        if( NTR > 0 ) :
            if( abs( 1 - probabilitySum ) > 1e-4 ) : warningList.append( 'decay probabilities sum to %.6f not 1.' % probabilitySum )

    if( offset != len( MF28 ) ) : raise Exception( 'Not allow lines of MF 28 data processed: offset = %d, len( MF28 ) = %d' % ( offset, len( MF28 ) ) )
    info.logs.write( '\n' )
    for warning in warningList : info.logs.write( "       WARNING: %s\n" % warning, stderrWriting = True )

    return( { 'PoPs' : PoPs, 'errors' : errors } )
