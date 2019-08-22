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

from __future__ import print_function, division

from pqu.PQU import PQU

from xData import values as valuesModule
from xData import XYs as XYsModule
from xData.uncertainty.physicalQuantity import uncertainty as uncertaintyModule
from xData.uncertainty.physicalQuantity import standard as standardModule

from PoPs import IDs as IDsPoPsModule
from PoPs import styles as stylesPoPsModule
from PoPs import alias as PoPsAliasModule

from PoPs.decays import decayData as decayDataModule
from PoPs.decays import probability as probabilityModule
from PoPs.decays import product as productModule
from PoPs.decays import Q as QModule
from PoPs.decays import spectrum as spectrumModule
from PoPs.decays import averageEnergy as averageEnergyModule

from PoPs.quantities import halflife as halflifeModule
from PoPs.quantities import mass as massModule

from PoPs.families import nucleus as nucleusModule

from PoPs.groups import isotope as isotopeModule
from PoPs.groups import misc as chemicalElementMiscPoPsModule

from fudge.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc
from fudge.legacy.converting import toGNDSMisc

energyUnit = 'eV'
alphaID = chemicalElementMiscPoPsModule.nucleusIDFromZAndA( 2, 4 )

decayType = { 0 : 'gamma',
              1 : 'beta-',
              2 : 'beta+ or e.c.', # or electron capture...
              3 : 'IT',
              4 : 'alpha',
              5 : 'n',
              6 : 'SF',
              7 : 'p',
              8 : 'e-',
              9 : 'xray',
             10 : 'unknown' }

decayParticles = {  0 : IDsPoPsModule.photon,
                    1 : IDsPoPsModule.electron,
                    2 : IDsPoPsModule.positron,
                    3 : None,
                    4 : alphaID,
                    5 : IDsPoPsModule.neutron,
                    6 : None,
                    7 : IDsPoPsModule.proton,
                    8 : IDsPoPsModule.electron,
                    9 : IDsPoPsModule.photon,
                   10 : None }

STYPProduct = { 0 : spectrumModule.gammaProduct,
                1 : spectrumModule.betaMinusProduct,
                2 : spectrumModule.betaPlusOrElectronCaptureProduct,
                4 : spectrumModule.alphaProduct,
                5 : spectrumModule.neutronProduct,
                6 : spectrumModule.SFProduct,
                7 : spectrumModule.protronProduct,
                8 : spectrumModule.discreteElectronProduct,
                9 : spectrumModule.xRayProduct }

def addProductToDecayMode( info, decayMode, decayParticleID, RFS = 0, lastDecayMode = False ) :

    if( lastDecayMode and ( RFS != 0 ) ) :
        isotopeNuclideID = chemicalElementMiscPoPsModule.nuclideIDFromIsotopeSymbolAndIndex( decayParticleID, RFS )
        metaStableID = PoPsAliasModule.metaStable.metaStableNameFromNuclearLevelNameAndMetaStableIndex( isotopeNuclideID, RFS )
        info.PoPs.add( PoPsAliasModule.metaStable( metaStableID, isotopeNuclideID, RFS ) )
        decayParticleID = metaStableID
    decayMode.products.add( decayMode.products.uniqueLabel( productModule.product( decayParticleID, decayParticleID ) ) )

def addDecayMode( info, decayMode, decayParticleIndex, parentZA, RFS, lastDecayMode ) :

    decayParticleID = decayParticles[decayParticleIndex]
    decayModeType = decayDataModule.decayModesParticle
    if( decayParticleID is None ) :
        if( decayParticleIndex == 3 ) :
            decayModeType = decayDataModule.decayModesIT
        elif( decayParticleIndex == 6 ) :
            decayModeType = decayDataModule.decayModesSF
        else :
            print( '    Oops', decayParticleIndex )
            return

    complete = decayParticleIndex not in [ 2, 6 ]
    _decay = decayDataModule.decay( str( len( decayMode.decayPath ) ), decayModeType, complete = complete )

    decayParticleZ = 0
    decayParticleA = 0
    if( decayParticleID is not None ) :
        if( decayParticleID != IDsPoPsModule.positron ) : addProductToDecayMode( info, _decay, decayParticleID )

        if( decayParticleID == IDsPoPsModule.electron ) :
            decayParticleZ = -1
            addProductToDecayMode( info, _decay, IDsPoPsModule.electronAntiNeutrino )
        elif( decayParticleID == IDsPoPsModule.positron ) :
            decayParticleZ = +1
        elif( decayParticleID == IDsPoPsModule.photon ) :
            pass
        elif( decayParticleID == IDsPoPsModule.neutron ) :
            decayParticleZ, decayParticleA = 0, 1
        elif( decayParticleID == IDsPoPsModule.proton ) :
            decayParticleZ, decayParticleA = 1, 1
        else :
            decayParticleZ, decayParticleA = 2, 4

    residialProductZA = parentZA - ( 1000 * decayParticleZ + decayParticleA )
    if( residialProductZA == 1 ) :
        residialProductID = IDsPoPsModule.neutron
    else :
        residialProductID = chemicalElementMiscPoPsModule.idFromZAndA( residialProductZA // 1000, residialProductZA % 1000 )
    if( decayParticleIndex != 6 ) : addProductToDecayMode( info, _decay, residialProductID, RFS, lastDecayMode )

    decayMode.decayPath.add( _decay )
    return( residialProductZA )

def getRTYP( line, offset = 0 ) :

    RTYP = line[offset:offset+11]
    RTYP = RTYP.split( '+' )[0]
    RTYP = RTYP.split( 'e' )[0]
    RTYP = RTYP.split( 'E' )[0]
    RTYP = RTYP.strip( '0' )
    initialDecay, otherDecays = RTYP.split( '.' )
    initialDecay = int( initialDecay )
    otherDecays = [ int( digit ) for digit in otherDecays ]

    RTYP_key = "%d.%s" % ( initialDecay, ''.join( [ "%s" % other for other in otherDecays ] ) )

    return( initialDecay, otherDecays, RTYP_key )

def addDiscreteSpectra( decayModes, RTYP_key, decayParticleIndex, discreteSpectra ) :

    def addInternalConversionCoefficients( discrete, RIC, dRIC, label ) :

        internalConversionCoefficients = spectrumModule.shell( RIC, label = label )
        addUncertainty( internalConversionCoefficients, dRIC )
        discrete.internalConversionCoefficients.add( internalConversionCoefficients )

    decayParticleID = decayParticles[decayParticleIndex]
    if( decayParticleID is None ) : raise Exception( 'decayParticleID is None, what is to be done?' )

    decayParticleLabel = STYPProduct[decayParticleIndex]
    spectra = decayModes[RTYP_key].spectra
    if( decayParticleLabel not in spectra ) :
        spectrum = spectrumModule.spectrum( decayParticleLabel, decayParticleID )
        spectra.add( spectrum )
    else :
        spectrum = spectra[decayParticleLabel]

    RTYP, TYPE, RI, dRI, RIS, dRIS = discreteSpectra['RTYP_plus']
    if( RIS != 0 ) : print( "STYP, RIS =", decayParticleIndex, RIS )

    if( TYPE == 1.0 ) :
        type = spectrumModule.transitionType.allowed
    elif( TYPE == 2.0 ) :
        type = spectrumModule.transitionType.firstForbidden
    elif( TYPE == 3.0 ) :
        type = spectrumModule.transitionType.secondForbidden
    else :
        type = None

    intensity = spectrumModule.intensity( RI )
    addUncertainty( intensity, dRI )

    energy = spectrumModule.energy( discreteSpectra['ER'], energyUnit )
    addUncertainty( energy, discreteSpectra['dER'] )

    positronIntensity = None
    if( decayParticleIndex == 0 ) :
        RIS = spectrumModule.internalPairFormationCoefficient( RIS )
    elif( decayParticleIndex == 2 ) :   # STYP=2 means interpret RIS differently
        positronIntensity = spectrumModule.positronEmissionIntensity( RIS )
        addUncertainty( positronIntensity, dRIS )
        RIS = None
    else  :
        RIS = None

    discrete = spectrumModule.discrete( intensity, energy, type = type, _internalPairFormationCoefficient = RIS,
            _positronEmissionIntensity = positronIntensity )
    if( 'RICC' in discreteSpectra ) :
        addInternalConversionCoefficients( discrete, discreteSpectra['RICC'], discreteSpectra['dRICC'], spectrumModule.shell.total )
    if( 'RICK' in discreteSpectra ) :
        addInternalConversionCoefficients( discrete, discreteSpectra['RICK'], discreteSpectra['dRICK'], spectrumModule.shell.KShell )
    if( 'RICL' in discreteSpectra ) :
        addInternalConversionCoefficients( discrete, discreteSpectra['RICL'], discreteSpectra['dRICL'], spectrumModule.shell.LShell )
    spectrum.append( discrete )

def addContinuumSpectrum( decayModes, RTYP_key, decayParticleIndex, regions, covariance ) :

    decayParticleID = decayParticles[decayParticleIndex]
    if( decayParticleID is None ) : raise TypeError( 'decayParticleID is None, what is to be done?' )

    decayParticleLabel = STYPProduct[decayParticleIndex]
    if RTYP_key not in decayModes: raise KeyError( "RTYP %s not found is list of decay modes, RTYP's in this evaluation are %s"%(RTYP_key,str(decayModes.keys())) )
    spectra = decayModes[RTYP_key].spectra
    if( decayParticleLabel not in spectra ) :
        spectrum = spectrumModule.spectrum( decayParticleLabel, decayParticleID )
        spectra.add( spectrum )
    else :
        spectrum = spectra[decayParticleLabel]

    if( len( regions ) == 1 ) :
        data = XYsModule.XYs1d( regions[0], axes = spectrumModule.continuum.defaultAxes( energyUnit ) )
    else :
        print( "len( regions ) > 1" )
        return
    continuum = spectrumModule.continuum( data )
    spectrum.append( continuum )

def addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, cls ) :

    if( len( aveDecayEnergies ) > 0 ) :
        if( aveDecayEnergies[0] >= 0 ) :
            averageEnergy = cls( aveDecayEnergies[0], energyUnit )
            addUncertainty( averageEnergy, aveDecayEnergies[1] )
            decayData.averageEnergies.add( averageEnergy )
        aveDecayEnergies = aveDecayEnergies[2:]

    return( aveDecayEnergies )

def addUncertainty( quantity, value ) :

    if( value == 0 ) : return
    quantity.uncertainty = uncertaintyModule.uncertainty( standardModule.standard( uncertaintyModule.double( value ) ) )

def printInfo( verbose, dataIndex, MT457MF8Data ) :

    if( verbose > 4 ) : print( '    dataIndex = ', dataIndex, MT457MF8Data[dataIndex], verbose )
    return( 1 )

def ITYPE_4( MTDatas, info, verbose = 0 ) :
    """Convert ENDF ITYPE 4 data (decay data) into PoPs."""

    errors = []

    info.PoPs.styles.add( stylesPoPsModule.evaluated( 'eval', '', info.library, info.libraryVersion, info.Date ) )

    # Add ENDF documentation
    info.PoPs.documentation = '\n'.join( MTDatas[451][1][4:-info.NXC] )

    # The atom containing the parent nucleus
    parentAtom = toGNDSMisc.getPoPsParticle( info, info.targetZA, name = None, levelIndex = info.levelIndex, level = info.level,
                                           levelUnit = energyUnit )
    parentAtom.mass.add(
        massModule.double( info.style, info.massTracker.amuMasses[info.targetZA], 'amu' )
    )
    decayParticle = parentAtom
    if(chemicalElementMiscPoPsModule.hasNucleus( parentAtom)) : decayParticle = parentAtom.nucleus

    if( 457 in MTDatas ) :
        dataIndex = 0
        MT457MF8Data = MTDatas[457][8]

        decayData = decayParticle.decayData

        # Read parent info block
        ZA, AWR, LIS, LISO, NST, NSP = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats( MT457MF8Data[dataIndex], intIndices = [ 2, 3, 4, 5 ] )
        dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )

        if( LIS != 0 ) :
            aliasID = "%s_m%s" % ( parentAtom.isotope.symbol, LISO )
            info.PoPs.add( PoPsAliasModule.metaStable( aliasID, parentAtom.id, LISO ) )

        HL, dHL, dummy, dummy, NC2, dummy = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats( MT457MF8Data[dataIndex], intIndices=[4] )
        dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )

        NC = NC2 // 2               # NC: number of decay energies
        if( NC not in [ 3, 17 ] ) : raise ValueError( "Number of decay energies must be 3 or 17, found %d" % NC )
        nlines = ( NC2 + 5 ) // 6
        aveDecayEnergies = []
        for index in range( nlines ) :
            aveDecayEnergies.extend( endfFileToGNDSMisc.sixFunkyFloatStringsToFloats( MT457MF8Data[dataIndex] ) )
            dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )

        if( len( [ aveDecayEnergy for aveDecayEnergy in aveDecayEnergies if( aveDecayEnergy != 0 ) ] ) > 0 ) :
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.lightParticles )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.electroMagneticRadiation )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.heavyParticles )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.betaMinus )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.betaPlus )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.AugerElectron )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.conversionElectron )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.gamma )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.xRay )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.internalBremsstrahlung )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.annihilation )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.alpha )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.recoil )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.spontaneousFission )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.fissionNeutrons )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.proton )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.neutrino )

        SPI, PAR, dum, dum, NDK6, NDK = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats( MT457MF8Data[dataIndex], intIndices = [1,2,3,4,5] )
        dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )
        if( -77.8 < SPI < -77.7 ) : SPI = None  # see page 180 of ENDF manual  # FIXME: Need "UNKNOWN" markup

        halflifeValue = None
        if( NST == 1 ) : halflifeValue = 'stable'
        if( PAR == 0 ) : PAR = None                     # BRB, Happens for dec-013_Al_040.endf
        toGNDSMisc.addParticleData( decayParticle, info, massValue = None, spinValue = SPI, parityValue = PAR,
                chargeValue = None, halflifeValue = halflifeValue )

        if( NST == 1 ) :                                # Nucleus is stable, nothing to do
            pass
        elif( NST == 0 ) :                              # Nucleus is unstable
            halflife = halflifeModule.double( info.PoPsLabel, HL, halflifeModule.baseUnit )
            addUncertainty( halflife, dHL )
            decayParticle.halflife.add( halflife )

            # Iterate over decay chains
            decayModes = {}
            for i1 in range( NDK ) :  # iterate over decay modes
                initialDecay, otherDecays, RTYP_key = getRTYP( MT457MF8Data[dataIndex] )
                RTYP, RFS, Q, dQ, BR, dBR = endfFileToGNDSMisc.sixFunkyFloatStringsToFloats( MT457MF8Data[dataIndex] )
                dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )
                if( verbose > 3 ) : print( 'RTYP = ', RTYP )
                RFS = int( RFS )                    # BRB, Isometic state of daughter

                allDecays = [ initialDecay ] + otherDecays
                decayMode = decayDataModule.decayMode( str( i1 ), ','.join( [ decayType[r_id] for r_id in allDecays ] ) )
                decayModes[RTYP_key] = decayMode

                Q = QModule.double( '', Q, energyUnit )
                addUncertainty( Q, dQ )
                decayMode.Q.add( Q )

                probability = probabilityModule.double( "BR", BR, '' )
                addUncertainty( probability, dBR )
                decayMode.probability.add( probability )
                                                                                                    # FIXME: what should we do for really small probabilities, so small ENDF gives 0, like spontaneous fission?
                n1 = len( otherDecays ) - 1
                residualZA = addDecayMode( info, decayMode, initialDecay, info.targetZA, RFS, -1 == n1 )
                for i2, otherDecay in enumerate( otherDecays ) :
                    residualZA = addDecayMode( info, decayMode, otherDecay, residualZA, RFS, i2 == n1 )
                decayData.decayModes.add( decayMode )

                                                    # LCON determines whether continuum spectra given: LCON 0: discrete only, 1: continuous only, 2: both
            for i1 in range( NSP ) :                # Iterate over spectra
                initialDecay, otherDecays, STYP_key = getRTYP( MT457MF8Data[dataIndex], 11 )
                dum, STYP, LCON, zero, six, NER = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats( MT457MF8Data[dataIndex], intIndices = [2,3,4,5] )
                dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )
                STYP = initialDecay
                specType = decayType[STYP]

                FD, dFD, ERave, dERave, FC, dFC = endfFileToGNDSMisc.sixFunkyFloatStringsToFloats( MT457MF8Data[dataIndex] ) # FD: discrete normalization, FC: continuum norm., ERave: average decay energy
                dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )

                if( LCON != 1 ) :                   # Discrete spectra given.
                    for i2 in range( NER ) :
                        ER, dER, dum, dum, NT, dum = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats( MT457MF8Data[dataIndex], intIndices = [2,3,4,5] )
                        dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )
                        discreteSpectra = { 'ER': ER, 'dER': dER }
                        if( LCON != 1 ) :
                            if( NT not in ( 6, 8, 12 ) ) : raise ValueError( "Unexpected NT=%d encountered!" % NT )
                            initialDecay, otherDecays, RTYP_key = getRTYP( MT457MF8Data[dataIndex] )
                            discreteSpectra['RTYP_plus'] = endfFileToGNDSMisc.sixFunkyFloatStringsToFloats( MT457MF8Data[dataIndex] )
                            dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )
                            if( verbose > 3 ) : print( 'STYP, RTYP = ', STYP, RTYP_key )
                            if( NT > 6 ) :
                                RICC, dRICC, RICK, dRICK, RICL, dRICL = endfFileToGNDSMisc.sixFunkyFloatStringsToFloats( MT457MF8Data[dataIndex] )
                                dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )
                                discreteSpectra['RICC'], discreteSpectra['dRICC'] = RICC, dRICC
                                if( NT > 8 ) :
                                    discreteSpectra['RICK'], discreteSpectra['dRICK'] = RICK, dRICK
                                    discreteSpectra['RICL'], discreteSpectra['dRICL'] = RICL, dRICL
                        addDiscreteSpectra( decayModes, RTYP_key, STYP, discreteSpectra )

                if( LCON != 0 ) :                   # Continuum spectra given.
                    initialDecay, otherDecays, RTYP_key = getRTYP( MT457MF8Data[dataIndex] )
                    try:
                        dataIndex, TAB1, regions = endfFileToGNDSMisc.getTAB1Regions( dataIndex, MT457MF8Data )
                        RTYP = TAB1['C1']
                        if( verbose > 3 ) : print( 'STYP, RTYP = ', STYP, RTYP )
                        LCOV = int( TAB1['L2'] )
                        if( len( regions ) != 1 ) : raise IndexError( 'len( regions ) = %s' % len( regions ) )
                        covariance = None
                        if( LCOV != 0 ) :
                            print( "    LCOV != 0" )
                            dataIndex, covariance = getList( dataIndex, MT457MF8Data )
                        addContinuumSpectrum( decayModes, RTYP_key, STYP, regions, covariance )
                    except Exception as err:
                        import sys
                        print('ERROR: '+str(err),file=sys.stderr)
        else :
            raise ValueError( "In decay data (MT=457,MF=8), NST should be 0 or 1. Found %d" % NST )

        # Add documentation to PoPs

    return( { 'PoPs' : info.PoPs, 'errors' : errors, 'info' : info } )
