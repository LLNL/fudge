# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from brownies.legacy.toENDF6 import endfFormats as endfFormatsModule
from PoPs import database as databaseModule
from PoPs.families import nuclide as nuclideModule

from fudge import documentation as documentationModule

from brownies.legacy.converting import massTracker as massTrackerModule
from brownies.legacy.endl import misc as miscENDLModule


def toENDF6( self, style, flags, verbosityIndent = '', useRedsFloatFormat = False,
             lineNumbers = True, **kwargs ) :

    _useRedsFloatFormat = endfFormatsModule.useRedsFloatFormat
    endfFormatsModule.useRedsFloatFormat = useRedsFloatFormat

    if( flags == {} ) : flags['verbosity'] = 0
    info = {}
    info['massTracker'] = massTrackerModule.massTracker()
    info['PoPs'] = self

    endfMFList = { 1 : { 451 : [] }, 8 : {}, 28 : {} }

    isAtomicRelaxation = False
    isSFY = False
    if( ( len( self.chemicalElements ) > 0 ) and ( self.chemicalElements[0].atomicData is not None ) ) :
        isAtomicRelaxation = True
    else :
        if( len( self.chemicalElements ) == 1 ) :
            isotopes = self.chemicalElements[0].isotopes
            if( len( isotopes ) == 1 ) :
                isotope = isotopes[0]
                if( len( isotope ) == 1 ) :
                    isSFY = len( isotope[0].fissionFragmentData.productYields ) > 0

    if( isAtomicRelaxation ) :
        Z = self.chemicalElements[0].Z
        MAT = 100 * Z
        info['ZA'] = 1000 * Z
        info['AWR'] = info['massTracker'].getElementalMassAMU( info['ZA'] )

        EMAX, LREL, LFI = 1e+11, 0, 0    # FIXME hard-coded
        levelEnergy_eV = STA = LIS = info['LISO'] = 0
        NSUB, NVER = 6, 8
    else :                                      # Decay or spontanepous fission yield sub-library 
        nuclideID = None
        for key in self.keys():
            if( isinstance( self[key], nuclideModule.particle ) ) :
                nuclideID = key
                nucleus = self[nuclideID].nucleus
                levelEnergy_eV = nucleus.energy.float( 'eV' )
                LIS = nucleus.index
                break
        if( ( nuclideID is None ) and ( list( self.keys( ) ) == [ 'n' ] ) ) :  # neutron is a special case
            nuclideID = 'n'
            nucleus = self[nuclideID]
            levelEnergy_eV = LIS = 0

        info['LISO'] = 0
        if (nuclideID in self.aliases):
            info['LISO'] = self[nuclideID].metaStableIndex
            nuclideID = self[nuclideID].pid
        Z, A, suffix, info['ZA'] = miscENDLModule.getZ_A_suffix_andZAFromName( nuclideID )
        info['AWR'] = self[nuclideID].mass.float('amu') / info['massTracker'].neutronMass

        STA = 0
        if( ( len( nucleus.decayData.decayModes ) > 0 ) or isSFY ) : STA = 1

        LFI = 0
        for decayMode in nucleus.decayData.decayModes:
            if decayMode.mode == 'SF': LFI = 1

        EMAX, LREL = 0, 0   # FIXME hard-coded
        NSUB, NVER = 4, 8
        if isSFY:
            NSUB = 5

    endfDoc = self.styles.getEvaluatedStyle().documentation.endfCompatible.body
    if( endfDoc is None ) :
        endfDoc = []
    else :
        endfDoc = endfDoc.split( '\n' )

    NLIB, NMOD = 0, 1   # FIXME hard-coded
    NFOR = 6            # ENDF-6

    docHeader = [ endfFormatsModule.endfHeadLine( info['ZA'], info['AWR'], -1, LFI, NLIB, NMOD ),
            endfFormatsModule.endfHeadLine( levelEnergy_eV, STA, LIS, info['LISO'], 0, NFOR ),
            endfFormatsModule.endfHeadLine( 0, EMAX, LREL, 0, NSUB, NVER ),
            endfFormatsModule.endfHeadLine( 0, 0, 0, 0, len( endfDoc ), 2 ) ]
    new_doc = documentationModule.documentation( 'endf', '\n'.join( docHeader + endfDoc ) )
    endfMFList[1][451] += endfFormatsModule.toEndfStringList( new_doc )

    if( isAtomicRelaxation ) :
        self.chemicalElements[0].atomicData.toENDF6( 533, endfMFList, flags, info, verbosityIndent )
    elif( isSFY ) :                             # Spontanepous fission yield sub-library
        self.chemicalElements[0].isotopes[0].nuclides[0].fissionFragmentData.productYields[0].toENDF6( endfMFList, flags, info, verbosityIndent )
    else:                                       # decay sub-library
        NST = 1 - STA
        NSP = sum([len(decayMode.spectra) for decayMode in nucleus.decayData.decayModes])

        if nucleus.halflife[0].value == 'stable':
            T_half = 0
        else:
            T_half = nucleus.halflife.float('s')
        T_half_uncert = 0
        if nucleus.halflife[0].uncertainty is not None:
            T_half_uncert = nucleus.halflife[0].uncertainty.form.value.float('s')

        SPI, PAR = -77.777, 0   # hacky way to indicate 'unknown'
        if len(nucleus.spin) > 0:
            SPI = nucleus.spin.float('hbar')
        if len(nucleus.parity) > 0:
            PAR = nucleus.parity[0].value

        NC = len(nucleus.decayData.averageEnergies)
        if NC==0: NC=3  # stable isotopes still list '0' for average energies

        endfMFList[8][457] = [
            endfFormatsModule.endfHeadLine(info['ZA'], info['AWR'], LIS, info['LISO'], NST, NSP),
            endfFormatsModule.endfContLine(T_half, T_half_uncert, 0, 0, 2*NC, 0)
        ]

        dataList = []
        for averageEnergy in nucleus.decayData.averageEnergies:
            dataList.append(averageEnergy.value)
            if averageEnergy.uncertainty is not None:
                dataList.append(averageEnergy.uncertainty.form.value.float('eV'))
            else:
                dataList.append(0)
        if len(dataList) == 0:
            dataList.append( [0]*6 )
        endfMFList[8][457] += endfFormatsModule.endfDataList(dataList)

        NDK = len(nucleus.decayData.decayModes)
        if NDK == 0:    # stick in a row of zeros
            endfMFList[8][457] += [
                endfFormatsModule.endfContLine(SPI, PAR, 0, 0, 6, 0),
                endfFormatsModule.endfDataLine([0]*6)
            ]
        else:
            endfMFList[8][457].append( endfFormatsModule.endfContLine(SPI, PAR, 0, 0, 6 * NDK, NDK) )

        nucleus.decayData.toENDF6( 457, endfMFList, flags, info, verbosityIndent )

    # FIXME currently parsing MAT out of the documentation since it's non-standard in decay sub-lib.
    # Better would be to either make it standard (same as rest of ENDF) or explicitly store it somewhere...
    MATLine = endfDoc[2].upper( )
    MAT = int( MATLine.split( 'MATERIAL' )[1].split( )[0] )

    return( endfFormatsModule.endfMFListToFinalFile( endfMFList, MAT, lineNumbers = lineNumbers ) )

databaseModule.database.toENDF6 = toENDF6
