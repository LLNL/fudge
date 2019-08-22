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

from xData import array as arrayModule

from PoPs import alias as aliasModule
from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.quantities import halflife as halflifeModule
from PoPs.groups import misc as chemicalElementMiscPoPsModule
from PoPs.fissionFragmentData import time as timeModule

from fudge.gnds.channelData.fissionFragmentData import fissionFragmentData as fissionFragmentDataModule
from fudge.gnds.channelData.fissionFragmentData import productYield as productYieldModule
from fudge.gnds.channelData.fissionFragmentData import duration as durationModule
from fudge.gnds.channelData.fissionFragmentData import incidentEnergy as incidentEnergyModule
from fudge.gnds.channelData.fissionFragmentData import yields as yieldsModule

from fudge.legacy.converting.ENDFToGNDS import endfFileToGNDSMisc as endfFileToGNDSMiscModule

from .. import toGNDSMisc as toGNDSMiscModule
from . import ENDF_ITYPE_0_Misc as ENDF_ITYPE_0_MiscModule

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

    allData = {}
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

def ITYPE_1( MTDatas, info, verbose = False ) :

    errors = []

    info.PoPs.documentation = '\n'.join( MTDatas[451][1][4:-info.NXC] )         # Add ENDF documentation

    parentAtom = toGNDSMiscModule.getPoPsParticle( info, info.targetZA, name = None, levelIndex = info.levelIndex, level = info.level, levelUnit = energyUnit )
    if( len( parentAtom.mass ) == 0 ) :
        parentAtom.mass.add( massModule.double( info.PoPsLabel, info.targetMass, quantityModule.stringToPhysicalUnit( 'amu' ) ) )
    if( len( parentAtom.nucleus.halflife ) == 0 ) :
        parentAtom.nucleus.halflife.add( halflifeModule.string( info.PoPsLabel, halflifeModule.UNSTABLE, halflifeModule.baseUnit ) )

    fissionFragmentData = fissionFragmentDataModule.fissionFragmentData( )

    if( 454 in MTDatas ) : parseMF454_459( info, fissionFragmentData, 454, MTDatas )
    if( 459 in MTDatas ) : parseMF454_459( info, fissionFragmentData, 459, MTDatas )

    return { 'fissionFragmentData' : fissionFragmentData, 'PoPs' : info.PoPs, 'errors' : errors, 'info' : info }
