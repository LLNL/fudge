# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from xData import enums as xDataEnumsModule
from xData import regions as regionsModule

from fudge.productData.distributions import energy as energyModule

from ... import gndsToENDF6 as gndsToENDF6Module
from ... import endfFormats as endfFormatsModule

#
# form
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    weight = targetInfo['delayedNubarWeight']
    energySubform = self.energySubform
    if( hasattr( energySubform, 'toENDF6' ) ) :
        if( isinstance( energySubform, ( energyModule.DiscreteGamma, energyModule.PrimaryGamma ) ) ) :
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
        print( 'WARNING: energy subform "%s" does not have method toENDF6' % energySubform.moniker )

energyModule.Form.toENDF6 = toENDF6

#
# XYs2d subform
#
def toENDF6( self, flags, targetInfo, weight = None, MT = None ) :

    NE = len( self )
    EInFactor = PQUModule.PQU(1, self.axes[-1].unit ).getValueAs('eV')
    if( weight is None ) : weight = [ [ self[0].outerDomainValue * EInFactor, 1.0 ],
                                      [ self[-1].outerDomainValue * EInFactor, 1.0 ] ]
    EInInterpolation = gndsToENDF6Module.gndsToENDF2PlusDInterpolationFlag( self.interpolation, self.interpolationQualifier )
    C1, C2, LAW, LANG = 0, 0, 0, 0
    if( MT == 527 ) : C1, C2, LAW, LANG = 11, 5.438675e-4, 1, 2
    ENDFDataList = [ endfFormatsModule.endfContLine( C1, C2, 0, 1, 1, len( weight ) ) ] + \
        endfFormatsModule.endfInterpolationList( [ len( weight ), 2 ] ) + endfFormatsModule.endfNdDataList( weight ) + \
        [ endfFormatsModule.endfContLine( 0, 0, LAW, LANG, 1, NE ) ] + endfFormatsModule.endfInterpolationList( [ NE, EInInterpolation ] )
    for energy in self :
        if( MT == 527 ) :
            ENDFDataList.append( endfFormatsModule.endfContLine( 0, energy.outerDomainValue * EInFactor, 0, 0,
                    2 * len( energy ), len( energy ) ) )
            ENDFDataList += endfFormatsModule.endfNdDataList( energy, xUnit = 'eV', yUnit = '1/eV' )
        else :
            if( isinstance( energy, regionsModule.Regions1d ) ) :
                interpolations, data = [], []
                for region in energy :
                    regionData = region.copyDataToXYs( )
                    if( len( data ) > 0 ) :
                        if( data[-1] == regionData[0] ) : regionData.pop( 0 )
                    data += regionData
                    interpolations.append( len( data ) )
                    interpolations.append( gndsToENDF6Module.gndsToENDFInterpolationFlag( region.interpolation ) )
                NR, NE = len( interpolations ) / 2, interpolations[-2]
                ENDFDataList.append( endfFormatsModule.endfContLine( 0, energy.outerDomainValue * EInFactor, 0, 0, NR, NE ) )
                ENDFDataList += endfFormatsModule.endfInterpolationList( interpolations )
                ENDFDataList += endfFormatsModule.endfNdDataList( data )
            else :
                interpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag( energy.interpolation )
                ENDFDataList.append( endfFormatsModule.endfContLine( 0, energy.outerDomainValue * EInFactor, 0, 0, 1, len( energy ) ) )
                ENDFDataList += endfFormatsModule.endfInterpolationList( [ len( energy ), interpolation ] )
                ENDFDataList += endfFormatsModule.endfNdDataList( energy, xUnit = 'eV', yUnit = '1/eV' )
    return( 1, ENDFDataList )

energyModule.XYs2d.toENDF6 = toENDF6
#
# Regions2d subform
#
def toENDF6_oneRegion( self, EInFactor, startingIndex = 0 ) :

    ENDFDataList = []
    for index, energy_in in enumerate( self[startingIndex:] ) :
        ENDFDataList += [ endfFormatsModule.endfContLine( 0., energy_in.outerDomainValue * EInFactor, 0, 0, 1, len( energy_in ) ) ]
        ENDFDataList += endfFormatsModule.endfInterpolationList(
                [ len( energy_in ), gndsToENDF6Module.gndsToENDFInterpolationFlag( energy_in.interpolation ) ] )
        ENDFDataList += endfFormatsModule.endfNdDataList( energy_in, xUnit = 'eV', yUnit = '1/eV' )
    EInInterpolation = gndsToENDF6Module.gndsToENDF2PlusDInterpolationFlag( self.interpolation, self.interpolationQualifier )
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
    if( weight is None ) : weight = [ [ self[0][0].outerDomainValue * EInFactor, 1.0 ], [ self[-1][-1].outerDomainValue * EInFactor, 1.0 ] ]
    NR = len( interpolations ) / 2
    NE = interpolations[-2]
    ENDFDataList = [ endfFormatsModule.endfContLine( 0, 0, 0, 1, 1, len( weight ) ) ] + \
            endfFormatsModule.endfInterpolationList( [ len( weight ), 2 ] ) + endfFormatsModule.endfNdDataList( weight ) + \
            [ endfFormatsModule.endfContLine( 0, 0, 0, 0, NR, NE ) ] + endfFormatsModule.endfInterpolationList( interpolations ) + \
            ENDFDataList
    return( 1, ENDFDataList )

energyModule.Regions2d.toENDF6 = toENDF6

#
# FunctionalBase
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
        interpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag( weight.interpolation )
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

energyModule.FunctionalBase.toENDF6 = toENDF6

#
# NBodyPhaseSpace
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    mass = self.mass.getValueAs( 'amu' ) / targetInfo['massTracker'].neutronMass
    ENDFDataList = [ endfFormatsModule.endfContLine( mass, 0, 0, 0, 0, self.numberOfProducts ) ]
    gndsToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, 6, xDataEnumsModule.Frame.centerOfMass, ENDFDataList )

energyModule.NBodyPhaseSpace.toENDF6 = toENDF6

#
# Weighted
#
def toENDF6( self, flags, targetInfo ) :

    return( self.functional.toENDF6( flags, targetInfo, weight = self.weight )[1] )

energyModule.Weighted.toENDF6 = toENDF6

#
# WeightedFunctionals
#
def toENDF6( self, flags, targetInfo, weight = None ) :                 # The weight is ignored in this method, only for compatibility with other toENDF6's.

    ENDFDataList = []
    for weight in self.weights : ENDFDataList += weight.toENDF6( flags, targetInfo )
    return( len( self.weights ), ENDFDataList )

energyModule.WeightedFunctionals.toENDF6 = toENDF6
