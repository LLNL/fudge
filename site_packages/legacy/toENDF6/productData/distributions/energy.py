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

from pqu import PQU as PQUModule

from xData import standards as standardsModule
from xData import regions as regionsModule

from fudge.gnd.productData.distributions import energy as energyModule

from ... import gndToENDF6 as gndToENDF6Module
from ... import endfFormats as endfFormatsModule

__metaclass__ = type
#
# form
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    weight = targetInfo['delayedNubarWeight']
    energySubform = self.energySubform
    if( hasattr( energySubform, 'toENDF6' ) ) :
        if( isinstance( energySubform, ( energyModule.discreteGamma, energyModule.primaryGamma ) ) ) :
            energySubform.toENDF6( flags, targetInfo )
            return
        if( MT in [ 527, 528 ] ) :
            NK, MF5 = energySubform.toENDF6( flags, targetInfo, weight = weight, MT = MT )
        else :
            NK, MF5 = energySubform.toENDF6( flags, targetInfo, weight = weight )

        if( MT == 455 ) :
            endfMFList[5][MT].append( MF5 ) 
        elif( MT == 528 ) :
            endfMFList[26][MT] = [ endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], 0, 0, NK, 0 ) ] + \
                    MF5 + [ endfFormatsModule.endfSENDLineNumber( ) ]
        else :
            endfMFList[5][MT]  = [ endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], 0, 0, NK, 0 ) ] + \
                    MF5 + [ endfFormatsModule.endfSENDLineNumber( ) ]
    else :
        print 'WARNING: energy subform "%s" does not have method toENDF6' % energySubform.moniker 

energyModule.form.toENDF6 = toENDF6

#
# XYs2d subform
#
def toENDF6( self, flags, targetInfo, weight = None, MT = None ) :

    NE = len( self )
    EInFactor = PQUModule.PQU(1, self.axes[-1].unit ).getValueAs('eV')
    if( weight is None ) : weight = [ [ self[0].value * EInFactor, 1.0 ],
                                      [ self[-1].value * EInFactor, 1.0 ] ]
    EInInterpolation = gndToENDF6Module.gndToENDF2PlusDInterpolationFlag( self.interpolation, self.interpolationQualifier )
    C1, C2, LAW, LANG = 0, 0, 0, 0
    if( MT == 527 ) : C1, C2, LAW, LANG = 11, 5.438675e-4, 1, 2
    ENDFDataList = [ endfFormatsModule.endfContLine( C1, C2, 0, 1, 1, len( weight ) ) ] + \
        endfFormatsModule.endfInterpolationList( [ len( weight ), 2 ] ) + endfFormatsModule.endfNdDataList( weight ) + \
        [ endfFormatsModule.endfContLine( 0, 0, LAW, LANG, 1, NE ) ] + endfFormatsModule.endfInterpolationList( [ NE, EInInterpolation ] )
    for energy in self :
        if( MT == 527 ) :
            ENDFDataList.append( endfFormatsModule.endfContLine( 0, energy.value * EInFactor, 0, 0,
                    2 * len( energy ), len( energy ) ) )
            ENDFDataList += endfFormatsModule.endfNdDataList( energy, xUnit = 'eV', yUnit = '1/eV' )
        else :
            if( isinstance( energy, regionsModule.regions1d ) ) :
                interpolations, data = [], []
                for region in energy :
                    regionData = region.copyDataToXYs( )
                    if( len( data ) > 0 ) :
                        if( data[-1] == regionData[0] ) : regionData.pop( 0 )
                    data += regionData
                    interpolations.append( len( data ) )
                    interpolations.append( gndToENDF6Module.gndToENDFInterpolationFlag( region.interpolation ) )
                NR, NE = len( interpolations ) / 2, interpolations[-2]
                ENDFDataList.append( endfFormatsModule.endfContLine( 0, energy.value * EInFactor, 0, 0, NR, NE ) )
                ENDFDataList += endfFormatsModule.endfInterpolationList( interpolations )
                ENDFDataList += endfFormatsModule.endfNdDataList( data )
            else :
                interpolation = gndToENDF6Module.gndToENDFInterpolationFlag( energy.interpolation )
                ENDFDataList.append( endfFormatsModule.endfContLine( 0, energy.value * EInFactor, 0, 0, 1, len( energy ) ) )
                ENDFDataList += endfFormatsModule.endfInterpolationList( [ len( energy ), interpolation ] )
                ENDFDataList += endfFormatsModule.endfNdDataList( energy, xUnit = 'eV', yUnit = '1/eV' )
    return( 1, ENDFDataList )

energyModule.XYs2d.toENDF6 = toENDF6
#
# regions2d subform
#
def toENDF6_oneRegion( self, EInFactor, startingIndex = 0 ) :

    ENDFDataList = []
    for index, energy_in in enumerate( self[startingIndex:] ) :
        ENDFDataList += [ endfFormatsModule.endfContLine( 0., energy_in.value * EInFactor, 0, 0, 1, len( energy_in ) ) ]
        ENDFDataList += endfFormatsModule.endfInterpolationList(
                [ len( energy_in ), gndToENDF6Module.gndToENDFInterpolationFlag( energy_in.interpolation ) ] )
        ENDFDataList += endfFormatsModule.endfNdDataList( energy_in, xUnit = 'eV', yUnit = '1/eV' )
    EInInterpolation = gndToENDF6Module.gndToENDF2PlusDInterpolationFlag( self.interpolation, self.interpolationQualifier )
    return( len( self[startingIndex:] ), EInInterpolation, ENDFDataList )

def toENDF6( self, flags, targetInfo, weight = None ) :

    interpolations, ENDFDataList, length = [], [], 0
    EInFactor = PQUModule.PQU( 1, self.axes[-1].unit ).getValueAs( 'eV' )
    for ridx, region in enumerate( self ) :
        startingIndex = 0
        if ridx>0 and region[0].copyDataToXYs() == self[ridx-1][-1].copyDataToXYs(): startingIndex = 1
        lengthSub, EInInterpolation, ENDFDataListSub = toENDF6_oneRegion( region, EInFactor, startingIndex )
        length += lengthSub
        interpolations += [ length, EInInterpolation ]
        ENDFDataList += ENDFDataListSub
    if( weight is None ) : weight = [ [ self[0][0].value * EInFactor, 1.0 ], [ self[-1][-1].value * EInFactor, 1.0 ] ]
    NR = len( interpolations ) / 2
    NE = interpolations[-2]
    ENDFDataList = [ endfFormatsModule.endfContLine( 0, 0, 0, 1, 1, len( weight ) ) ] + \
            endfFormatsModule.endfInterpolationList( [ len( weight ), 2 ] ) + endfFormatsModule.endfNdDataList( weight ) + \
            [ endfFormatsModule.endfContLine( 0, 0, 0, 0, NR, NE ) ] + endfFormatsModule.endfInterpolationList( interpolations ) + \
            ENDFDataList
    return( 1, ENDFDataList )

energyModule.regions2d.toENDF6 = toENDF6

#
# functionalBase
#
def toENDF6( self, flags, targetInfo, weight = None ) :

    U, EFL, EFH = 0, 0, 0
    if( self.LF == 12 ) :
        EFL, EFH = self.EFL.getValueAs( 'eV' ), self.EFH.getValueAs( 'eV' )
    else :
        U = self.U.getValueAs( 'eV' )
    parameter1 = self.parameter1.data
    energyFactor = float( PQUModule.PQU( '1 eV' ) / PQUModule.PQU( 1, parameter1.axes[1].unit ) )
    if( weight is None ) :
        weight = [ [ energyFactor * parameter1[0][0], 1.0 ], [ energyFactor * parameter1[-1][0], 1.0 ] ]
        interpolation = 2
    elif( hasattr( weight, 'axes' ) ):
        interpolation = gndToENDF6Module.gndToENDFInterpolationFlag( weight.interpolation )
    else:
        interpolation = 2
    ENDFDataList = [ endfFormatsModule.endfContLine( U, 0, 0, self.LF, 1, len( weight ) ) ] + \
        endfFormatsModule.endfInterpolationList( [ len( weight ), interpolation ] ) + endfFormatsModule.endfNdDataList( weight )
    ENDFDataList += endfFormatsModule.toTAB1( parameter1, 'eV', 'eV', C1 = EFL, C2 = EFH )
    if( self.parameter2 is not None ) :
        parameter2 = self.parameter2.data
        yUnit = ''
        if( self.LF in [ 5, 11 ] ) : yUnit = '1/eV'
        ENDFDataList += endfFormatsModule.toTAB1( parameter2, 'eV', yUnit )
    return( 1, ENDFDataList )

energyModule.functionalBase.toENDF6 = toENDF6

#
# NBodyPhaseSpace
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    mass = self.numberOfProductsMasses.getValueAs( 'amu' ) / targetInfo['massTracker'].neutronMass
    ENDFDataList = [ endfFormatsModule.endfContLine( mass, 0, 0, 0, 0, self.numberOfProducts ) ]
    gndToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, 6, standardsModule.frames.centerOfMassToken, ENDFDataList )

energyModule.NBodyPhaseSpace.toENDF6 = toENDF6

#
# weighted
#
def toENDF6( self, flags, targetInfo ) :

    return( self.functional.toENDF6( flags, targetInfo, weight = self.weight )[1] )

energyModule.weighted.toENDF6 = toENDF6

#
# weightedFunctionals
#
def toENDF6( self, flags, targetInfo, weight = None ) :                 # The weight is ignored in this method, only for compatibility with other toENDF6's.

    ENDFDataList = []
    for weight in self.weights : ENDFDataList += weight.toENDF6( flags, targetInfo )
    return( len( self.weights ), ENDFDataList )

energyModule.weightedFunctionals.toENDF6 = toENDF6
