# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
For translating ENDF ITYPE=4 data (radioactive decay data)
"""

from xData import XYs1d as XYs1dModule
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

from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule

from brownies.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc
from brownies.legacy.converting import toGNDSMisc

energyUnit = 'eV'
alphaID = chemicalElementMiscPoPsModule.nucleusIDFromZAndA( 2, 4 )

decayType = { 0 : 'gamma',
              1 : 'beta-',
              2 : 'beta+ or e.c.', # or electron capture...
              3 : 'IT',
              4 : 'alpha',
              5 : IDsPoPsModule.neutron,
              6 : 'SF',
              7 : IDsPoPsModule.proton,
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
        metaStableID = PoPsAliasModule.MetaStable.metaStableNameFromNuclearLevelNameAndMetaStableIndex(isotopeNuclideID, RFS)
        info.PoPs.add(PoPsAliasModule.MetaStable(metaStableID, isotopeNuclideID, RFS))
        decayParticleID = metaStableID
    decayMode.products.add( decayMode.products.uniqueLabel( productModule.Product( decayParticleID, decayParticleID ) ) )

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
    _decay = decayDataModule.Decay( str( len( decayMode.decayPath ) ), decayModeType, complete = complete )

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

        internalConversionCoefficients = spectrumModule.Shell(RIC, label=label, unit='')
        addUncertainty( internalConversionCoefficients, dRIC )
        discrete.internalConversionCoefficients.add( internalConversionCoefficients )

    decayParticleID = decayParticles[decayParticleIndex]
    if( decayParticleID is None ) : raise Exception( 'decayParticleID is None, what is to be done?' )

    decayParticleLabel = STYPProduct[decayParticleIndex]
    spectra = decayModes[RTYP_key].spectra
    if( decayParticleLabel not in spectra ) :
        spectrum = spectrumModule.Spectrum( decayParticleLabel, decayParticleID )
        spectra.add( spectrum )
    else :
        spectrum = spectra[decayParticleLabel]

    RTYP, TYPE, RI, dRI, RIS, dRIS = discreteSpectra['RTYP_plus']
    if( RIS != 0 ) : print( "STYP, RIS =", decayParticleIndex, RIS )

    if( TYPE == 1.0 ) :
        type = spectrumModule.TransitionType.allowed
    elif( TYPE == 2.0 ) :
        type = spectrumModule.TransitionType.firstForbidden
    elif( TYPE == 3.0 ) :
        type = spectrumModule.TransitionType.secondForbidden
    else :
        type = spectrumModule.TransitionType.none

    intensity = spectrumModule.Intensity( RI )
    addUncertainty( intensity, dRI )

    energy = spectrumModule.Energy( discreteSpectra['ER'], energyUnit )
    addUncertainty( energy, discreteSpectra['dER'] )

    positronIntensity = None
    if( decayParticleIndex == 0 and RIS != 0 ) :
        RIS = spectrumModule.InternalPairFormationCoefficient( RIS )
    elif( decayParticleIndex == 2 and RIS != 0 ) :   # STYP=2 means interpret RIS differently
        positronIntensity = spectrumModule.PositronEmissionIntensity( RIS )
        addUncertainty( positronIntensity, dRIS )
        RIS = None
    else  :
        RIS = None

    discrete = spectrumModule.Discrete( intensity, energy, type = type, _internalPairFormationCoefficient = RIS,
            _positronEmissionIntensity = positronIntensity )
    if( 'RICC' in discreteSpectra ) :
        addInternalConversionCoefficients( discrete, discreteSpectra['RICC'], discreteSpectra['dRICC'], spectrumModule.Shell.total )
    if( 'RICK' in discreteSpectra ) :
        addInternalConversionCoefficients( discrete, discreteSpectra['RICK'], discreteSpectra['dRICK'], spectrumModule.Shell.KShell )
    if( 'RICL' in discreteSpectra ) :
        addInternalConversionCoefficients( discrete, discreteSpectra['RICL'], discreteSpectra['dRICL'], spectrumModule.Shell.LShell )
    spectrum.append( discrete )

def addContinuumSpectrum( decayModes, RTYP_key, decayParticleIndex, regions, covariance ) :

    decayParticleID = decayParticles[decayParticleIndex]
    if( decayParticleID is None ) : raise TypeError( 'decayParticleID is None, what is to be done?' )

    decayParticleLabel = STYPProduct[decayParticleIndex]
    if RTYP_key not in decayModes: raise KeyError( "RTYP %s not found is list of decay modes, RTYP's in this evaluation are %s"%(RTYP_key,str(decayModes.keys())) )
    spectra = decayModes[RTYP_key].spectra
    if( decayParticleLabel not in spectra ) :
        spectrum = spectrumModule.Spectrum( decayParticleLabel, decayParticleID )
        spectra.add( spectrum )
    else :
        spectrum = spectra[decayParticleLabel]

    if( len( regions ) == 1 ) :
        data = XYs1dModule.XYs1d( regions[0], axes = spectrumModule.Continuum.defaultAxes( energyUnit ) )
    else :
        print( "len( regions ) > 1" )
        return
    continuum = spectrumModule.Continuum( data )
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
    quantity.uncertainty = uncertaintyModule.Uncertainty( standardModule.Standard( uncertaintyModule.Double( value ) ) )

def printInfo( verbose, dataIndex, MT457MF8Data ) :

    if( verbose > 4 ) : print( '    dataIndex = ', dataIndex, MT457MF8Data[dataIndex], verbose )
    return( 1 )

def ITYPE_4( MTDatas, info, verbose = 0 ) :
    """Convert ENDF ITYPE 4 data (decay data) into PoPs."""

    errors = []

    evalStyle = stylesPoPsModule.Evaluated( 'eval', '', info.library, info.libraryVersion, info.Date )
    evalStyle.documentation.endfCompatible.body = '\n'.join( MTDatas[451][1][4:-info.NXC] )     # Add ENDF documentation
    endfFileToGNDSMisc.completeDocumentation(info, evalStyle.documentation)

    info.PoPs.styles.add( evalStyle )

    # The atom containing the parent nucleus
    parentAtom = toGNDSMisc.getPoPsParticle(info, info.targetZA, name = None, levelIndex = info.levelIndex, level = info.level,
                                            levelUnit = energyUnit)
    parentAtom.mass.add(
        massModule.Double( info.style, info.massTracker.amuMasses[info.targetZA], 'amu' )
    )
    decayParticle = parentAtom
    if(chemicalElementMiscPoPsModule.hasNucleus( parentAtom)) : decayParticle = parentAtom.nucleus

    if( 457 in MTDatas ) :
        dataIndex = 0
        MT457MF8Data = MTDatas[457][8]

        decayData = decayParticle.decayData

        # Read parent info block
        ZA, AWR, LIS, LISO, NST, NSP = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats(MT457MF8Data[dataIndex], intIndices = [2, 3, 4, 5])
        dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )

        if( LIS != 0 ) :
            aliasID = "%s_m%s" % ( parentAtom.isotope.symbol, LISO )
            info.PoPs.add(PoPsAliasModule.MetaStable(aliasID, parentAtom.id, LISO))

        HL, dHL, dummy, dummy, NC2, dummy = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats(MT457MF8Data[dataIndex], intIndices=[4])
        dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )

        NC = NC2 // 2               # NC: number of decay energies
        if( NC not in [ 3, 17 ] ) : raise ValueError( "Number of decay energies must be 3 or 17, found %d" % NC )
        nlines = ( NC2 + 5 ) // 6
        aveDecayEnergies = []
        for index in range( nlines ) :
            aveDecayEnergies.extend(endfFileToGNDSMisc.sixFunkyFloatStringsToFloats(MT457MF8Data[dataIndex]))
            dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )

        if( len( [ aveDecayEnergy for aveDecayEnergy in aveDecayEnergies if( aveDecayEnergy != 0 ) ] ) > 0 ) :
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.LightParticles )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.ElectroMagneticRadiation )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.HeavyParticles )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.BetaMinus )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.BetaPlus )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.AugerElectron )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.ConversionElectron )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.Gamma )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.XRay )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.InternalBremsstrahlung )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.Annihilation )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.Alpha )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.Recoil )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.SpontaneousFission )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.FissionNeutrons )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.Proton )
            aveDecayEnergies = addAverageDecayEnergyIfPresent( decayData, aveDecayEnergies, averageEnergyModule.Neutrino )

        SPI, PAR, dum, dum, NDK6, NDK = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats(MT457MF8Data[dataIndex], intIndices = [1, 2, 3, 4, 5])
        dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )
        if( -77.8 < SPI < -77.7 ) : SPI = None  # see page 180 of ENDF manual  # FIXME: Need "UNKNOWN" markup

        halflifeValue = None
        if( NST == 1 ) : halflifeValue = 'stable'
        if( PAR == 0 ) : PAR = None                     # BRB, Happens for dec-013_Al_040.endf
        toGNDSMisc.addParticleData(decayParticle, info, massValue = None, spinValue = SPI, parityValue = PAR,
                                   chargeValue = None, halflifeValue = halflifeValue)

        if( NST == 1 ) :                                # Nucleus is stable, nothing to do
            pass
        elif( NST == 0 ) :                              # Nucleus is unstable
            halflife = halflifeModule.Double( info.PoPsLabel, HL, halflifeModule.baseUnit )
            addUncertainty( halflife, dHL )
            decayParticle.halflife.add( halflife )

            # Iterate over decay chains
            decayModes = {}
            for i1 in range( NDK ) :  # iterate over decay modes
                initialDecay, otherDecays, RTYP_key = getRTYP( MT457MF8Data[dataIndex] )
                RTYP, RFS, Q, dQ, BR, dBR = endfFileToGNDSMisc.sixFunkyFloatStringsToFloats(MT457MF8Data[dataIndex])
                dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )
                if( verbose > 3 ) : print( 'RTYP = ', RTYP )
                RFS = int( RFS )                    # BRB, Isometic state of daughter

                allDecays = [ initialDecay ] + otherDecays
                decayMode = decayDataModule.DecayMode( str( i1 ), ','.join( [ decayType[r_id] for r_id in allDecays ] ) )
                decayModes[RTYP_key] = decayMode

                Q = QModule.Double( info.PoPsLabel, Q, energyUnit )
                addUncertainty( Q, dQ )
                decayMode.Q.add( Q )

                probability = probabilityModule.Double( "BR", BR, '' )
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
                dum, STYP, LCON, zero, six, NER = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats(MT457MF8Data[dataIndex], intIndices = [2, 3, 4, 5])
                dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )
                STYP = initialDecay
                specType = decayType[STYP]

                FD, dFD, ERave, dERave, FC, dFC = endfFileToGNDSMisc.sixFunkyFloatStringsToFloats(MT457MF8Data[dataIndex]) # FD: discrete normalization, FC: continuum norm., ERave: average decay energy
                dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )

                if( LCON != 1 ) :                   # Discrete spectra given.
                    for i2 in range( NER ) :
                        ER, dER, dum, dum, NT, dum = endfFileToGNDSMisc.sixFunkyFloatStringsToIntsAndFloats(MT457MF8Data[dataIndex], intIndices = [2, 3, 4, 5])
                        dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )
                        discreteSpectra = { 'ER': ER, 'dER': dER }
                        if( LCON != 1 ) :
                            if( NT not in ( 6, 8, 12 ) ) : raise ValueError( "Unexpected NT=%d encountered!" % NT )
                            initialDecay, otherDecays, RTYP_key = getRTYP( MT457MF8Data[dataIndex] )
                            discreteSpectra['RTYP_plus'] = endfFileToGNDSMisc.sixFunkyFloatStringsToFloats(MT457MF8Data[dataIndex])
                            dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )
                            if( verbose > 3 ) : print( 'STYP, RTYP = ', STYP, RTYP_key )
                            if( NT > 6 ) :
                                RICC, dRICC, RICK, dRICK, RICL, dRICL = endfFileToGNDSMisc.sixFunkyFloatStringsToFloats(MT457MF8Data[dataIndex])
                                dataIndex += printInfo( verbose, dataIndex, MT457MF8Data )
                                discreteSpectra['RICC'], discreteSpectra['dRICC'] = RICC, dRICC
                                if( NT > 8 ) :
                                    discreteSpectra['RICK'], discreteSpectra['dRICK'] = RICK, dRICK
                                    discreteSpectra['RICL'], discreteSpectra['dRICL'] = RICL, dRICL
                        addDiscreteSpectra( decayModes, RTYP_key, STYP, discreteSpectra )

                if( LCON != 0 ) :                   # Continuum spectra given.
                    initialDecay, otherDecays, RTYP_key = getRTYP( MT457MF8Data[dataIndex] )
                    try:
                        dataIndex, TAB1, regions = endfFileToGNDSMisc.getTAB1Regions(dataIndex, MT457MF8Data)
                        RTYP = TAB1['C1']
                        if( verbose > 3 ) : print( 'STYP, RTYP = ', STYP, RTYP )
                        LCOV = int( TAB1['L2'] )
                        if( len( regions ) != 1 ) : raise IndexError( 'len( regions ) = %s' % len( regions ) )
                        covariance = None
                        if( LCOV != 0 ) :
                            print( "    LCOV != 0" )
                            dataIndex, covariance = endfFileToGNDSMisc.getList(dataIndex, MT457MF8Data)
                        addContinuumSpectrum( decayModes, RTYP_key, STYP, regions, covariance )
                    except Exception as err:
                        import sys
                        print('ERROR: '+str(err),file=sys.stderr)
        else :
            raise ValueError( "In decay data (MT=457,MF=8), NST should be 0 or 1. Found %d" % NST )

    return( { 'PoPs' : info.PoPs, 'errors' : errors, 'info' : info } )
