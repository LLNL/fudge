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

from site_packages.legacy.toENDF6 import endfFormats as endfFormatsModule
from PoPs import database as databaseModule
from PoPs.families import nuclide as nuclideModule

from fudge.gnds import documentation as documentationModule

from fudge.legacy.converting import endf_endl as endf_endlModule
from fudge.legacy.converting import massTracker as massTrackerModule
from fudge.legacy.endl import misc as miscENDLModule

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
        if nuclideID is None and self.keys() == ['n']:  # neutron is a special case
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

    endfDoc = self.documentation.split('\n')    # FIXME documentation should be in a class, not a raw string

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
        pass
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
