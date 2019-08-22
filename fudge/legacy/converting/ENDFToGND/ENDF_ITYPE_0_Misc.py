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

epsilonExponent = -10

import sys
import math

from pqu import PQU as PQUModule
from fudge.core.math import matrix as gndMatrix
from fudge.core.math import linearAlgebra as linearAlgebraModule

from fudge.core.utilities import brb

import xData.standards as standardsModule
import xData.axes as axesModule
import xData.values as valuesModule
import xData.XYs as XYsModule
import xData.regions as regionsModule
import xData.link as linkModule
import xData.multiD_XYs as multiD_XYsModule
import xData.array as arrayModule
import xData.gridded as griddedModule
import xData.uncertainties as uncertaintiesModule

import fudge.gnd.alias as aliasModule
import fudge.gnd.xParticle as xParticleModule
import fudge.gnd.tokens as tokensModule
import fudge.gnd.resonances as resonancesModule
import fudge.gnd.sums as sumsModule

import fudge.gnd.covariances.base as covarianceBaseModule
import fudge.gnd.covariances.covarianceSuite as covarianceSuiteModule
import fudge.gnd.covariances.section as covarianceSectionModule
import fudge.gnd.covariances.summed as covarianceSummedModule
import fudge.gnd.covariances.mixed as covarianceMixedModule
import fudge.gnd.covariances.modelParameters as covarianceModelParametersModule
import fudge.gnd.covariances.distributions as covarianceDistributionsModule

import fudge.gnd.channels as channelsModule

import fudge.gnd.channelData.fissionEnergyReleased as fissionEnergyReleasedModule

import fudge.gnd.reactions.reaction as reactionModule
import fudge.gnd.reactions.production as productionModule
import fudge.gnd.reactionData.crossSection as crossSectionModule

import fudge.gnd.productData.multiplicity as multiplicityModule
import fudge.gnd.productData.distributions.base as distributionsBaseModule
import fudge.gnd.productData.distributions.unknown as unknownModule
import fudge.gnd.productData.distributions.unspecified as unspecifiedModule
import fudge.gnd.productData.distributions.angular as angularModule
import fudge.gnd.productData.distributions.energy as energyModule
import fudge.gnd.productData.distributions.uncorrelated as uncorrelatedModule
import fudge.gnd.productData.distributions.angularEnergy as angularEnergyModule
import fudge.gnd.productData.distributions.energyAngular as energyAngularModule
import fudge.gnd.productData.distributions.KalbachMann as KalbachMannModule
import fudge.gnd.productData.distributions.reference as referenceModule

from fudge.legacy.converting import endf_endl, toGNDMisc
import endfFileToGNDMisc

frames = { 1 : standardsModule.frames.labToken, 2 : standardsModule.frames.centerOfMassToken }
ENDF_Accuracy = endfFileToGNDMisc.ENDF_Accuracy
FUDGE_EPS = endfFileToGNDMisc.FUDGE_EPS
productNameToZA = { 'n' : 1, 'H1' : 1001, 'H2' : 1002, 'H3' : 1003, 'He3' : 2003, 'He4' : 2004, 'gamma' : 0 }
lightIsotopeNames = [ 'n', 'H1', 'H2', 'H3', 'He3', 'He4' ]

# The following is a kludge for ENDF/B-VII.1 files that are missing some MF 8 data. Namely,
#    n-020_Ca_040, n-020_Ca_042, n-020_Ca_043, n-020_Ca_044, n-020_Ca_046, n-020_Ca_048,
#    n-082_Pb_204, n-082_Pb_206, n-082_Pb_207 and n-082_Pb_208.
# Data for the following meta-stables are bogus as I could not find values:
# Re182_m1, Ir186_m1, Au185_m1, Hg187_m1, Hg189_m1, Hg191_m1, Tl186_m1, Tl188_m1, Tl190_m1, Tl191_m1, Tl193_m1, Tl194_m1, Pb193_m1 and Pb199_m1.
# Meta-stable energies in eV. Data as ( levelIndex, levelEnergy )
metaStableData = {
    "Al26_m1" : ( 1, 228305. ),     "Cl34_m1" : ( 1, 146360. ),     "W179_m1" : ( 2, 221926. ),         "W183_m1" : ( 7, 309493. ),
    "W185_m1" : ( 6, 197383. ),     "Re182_m1" : ( 1, 1000. ),      "Re184_m1" : ( 5, 188010. ),        "Re186_m1" : ( 4, 149000. ),
    "Re188_m1" : ( 7, 172069. ),    "Re190_m1" : ( 3, 210000. ),    "Os181_m1" : ( 1, 49200. ),         "Os183_m1" : ( 2, 170710. ),
    "Os189_m1" : ( 1, 30812. ),     "Os191_m1" : ( 1, 74382. ),     "Ir186_m1" : ( 1, 1000. ),          "Ir190_m1" : ( 2, 26100. ),
    "Ir190_m2" : ( 37, 376400. ) ,  "Ir191_m1" : ( 3, 171240. ),    "Ir192_m1" : ( 3, 56720. ),         "Ir192_m2" : ( 16, 168140. ),
    "Ir193_m1" : ( 2, 80239. ),     "Ir194_m1" : ( 7, 147072. ),    "Ir195_m1" : ( 2, 1e5 ),            "Ir196_m1" : ( 4, 410000. ),
    "Ir197_m1" : ( 2, 115000. ),    "Pt183_m1" : ( 1, 34500. ),     "Pt185_m1" : ( 2, 103410. ),        "Pt193_m1" : ( 5, 149780. ),
    "Pt195_m1" : ( 7, 259300. ),    "Pt197_m1" : ( 9, 399590. ),    "Pt199_m1" : ( 8, 424000. ),        "Au185_m1" : ( 1, 1000. ),
    "Au187_m1" : ( 2, 120510. ),    "Au189_m1" : ( 3, 247230. ),    "Au193_m1" : ( 4, 290190. ),        "Au195_m1" : ( 4, 318580. ),
    "Au196_m1" : ( 3, 84656. ),     "Au197_m1" : ( 4, 409150. ),    "Au200_m1" : ( 11, 962000. ),       "Hg185_m1" : ( 4, 99300. ),
    "Hg187_m1" : ( 1, 1000. ),      "Hg189_m1" : ( 2, 1000. ),      "Hg191_m1" : ( 1, 1000. ),          "Hg193_m1" : ( 3, 140760. ),
    "Hg195_m1" : ( 3, 176070. ),    "Hg197_m1" : ( 4, 298930. ),    "Hg199_m1" : ( 7, 532480. ),        "Tl186_m1" : ( 1, 1000. ),
    "Tl187_m1" : ( 2, 335009. ),    "Tl188_m1" : ( 1, 1000. ),      "Tl189_m1" : ( 1, 257600. ),        "Tl190_m1" : ( 1, 1000. ),
    "Tl191_m1" : ( 1, 1000. ),      "Tl192_m1" : ( 1, 156000. ),    "Tl193_m1" : ( 1, 1000. ),          "Tl194_m1" : ( 1, 1000. ),
    "Tl195_m1" : ( 2, 482630. ),    "Tl196_m1" : ( 6, 394200. ),    "Tl198_m1" : ( 7, 543500. ),        "Pb187_m1" : ( 1, 81000. ),
    "Pb191_m1" : ( 1, 138000. ),    "Pb193_m1" : ( 1, 1000. ),      "Pb195_m1" : ( 2, 202900. ),        "Pb197_m1" : ( 2, 319310. ),
    "Pb199_m1" : ( 1, 1000. ),      "Pb201_m1" : ( 4, 629100. ),    "Pb202_m1" : ( 13, 2169830. ),      "Pb203_m1" : ( 6, 825200. ),
    "Sc42_m1"  : ( 2, 616300. ),    "Sc44_m1"  : ( 1, 271240. ),    "Sc46_m1"  : ( 1, 142528. ),
    "Mn50_m1"  : ( 1, 225280. ),    "Mn52_m1"  : ( 1, 377749. ),
    "Se73_m1"  : ( 1,  25710. ),    "Se77_m1"  : ( 1, 161922.3 ),   "Se79_m1"  : ( 1, 95770. ),         "Se81_m1"  : ( 1, 103000. ),
    "Br74_m1"  : ( 1,  13580. ),    "Br76_m1"  : ( 1, 102580. ),    "Br77_m1"  : ( 1, 105860. ),        "Br79_m1"  : ( 1, 207610. ),
    "Br80_m1"  : ( 1,  85843. ),    "Br82_m1"  : ( 1,  45949.2 ),   "Br84_m1"  : ( 1, 320000. ),
    "Kr79_m1"  : ( 1, 129770. ),    "Kr81_m1"  : ( 1, 190640. ),    "Kr83_m1"  : ( 1,  41557.5 ),       "Kr85_m1"  : ( 1, 304871. ),
    "Rb81_m1"  : ( 1, 86317.1 ),    "Rb82_m1"  : ( 1,  69000. ),    "Rb84_m1"  : ( 1, 463591. ),        "Rb86_m1"  : ( 1, 556051. ),
    "Rb90_m1"  : ( 1, 106901. ),
    "Sr83_m1"  : ( 1, 259151. ),    "Sr85_m1"  : ( 1, 238791. ),    "Sr87_m1"  : ( 1, 388533. ),
    "Y83_m1"   : ( 1,  62047.1 ),   "Y85_m1"   : ( 1, 196800. ),    "Y86_m1"   : ( 1, 218211. ),        "Y87_m1"   : ( 1, 380821. ),
    "Y89_m1"   : ( 1, 908971. ),    "Y90_m1"   : ( 1, 681671. ),    "Y91_m1"   : ( 1, 555580. ),
    "Zr85_m1"  : ( 1, 292281. ),    "Zr87_m1"  : ( 1, 335841. ),    "Zr89_m1"  : ( 1, 587821. ),
    "Nb90_m1"  : ( 1, 124670. ),    "Nb91_m1"  : ( 1, 104601. ),    "Nb92_m1"  : ( 1, 135501. ),        "Nb93_m1"  : ( 1,  30770.1 ),
    "Mo91_m1"  : ( 1, 653011. ),
    "Tc90_m1"  : ( 1, 100000. ),    "Tc91_m1"  : ( 1, 139300. ),    "Tc93_m1"  : ( 1, 391840. ),        "Tc94_m1"  : ( 1,  76000. ),
    "Tc95_m1"  : ( 1,  38900. ),    "Tc96_m1"  : ( 1,  34230. ),
    "Bi192_m1" : ( 1, 147000. ),    "Bi193_m1" : ( 1, 308000. ),    "Bi194_m1" : ( 1, 109000. ),        "Bi194_m2"  : ( 2, 200000. ),
    "Bi195_m1" : ( 1, 401000. ),    "Bi196_m1" : ( 1, 271000. ),    "Bi197_m1" : ( 1, 500000. ),        "Bi198_m1"  : ( 1, 100000. ),
    "Bi198_m2" : ( 2,  248500. ),   "Bi199_m1" : ( 1, 667000. ),    "Bi200_m1" : ( 1, 428200. ),        "Bi201_m1"  : ( 1, 846350. ) }

class myIter:
    """Iterator that keeps track of line number."""

    def __init__( self, iterable ):

        self.index = 0
        self.iterable = iter( iterable )
        self.length = len( iterable )

    def next( self ) :

        next_ = self.iterable.next()
        self.index += 1
        return next_

# Two useful shortcuts for reading ENDF data.
funkyF = endfFileToGNDMisc.sixFunkyFloatStringsToFloats
def funkyFI( a, logFile = sys.stderr ) :   # read ENDF line with 2 floats and 4 ints

    return( endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( a, [ 2, 3, 4, 5 ], logFile = logFile ) )

# create some custom Exceptions:
class BadResonances( Exception ) : pass
class BadCovariance( Exception ) : pass

def calculateZA( ZACompound, ZAOther, minus = True ) :
    """This function handles the removal (or addition) of ZAOther to ZACompound include natural compound (but not a natural other)."""

    if( ( ZACompound % 1000 ) == 0 ) : ZAOther = 1000 * ( ZAOther // 1000 )
    if( minus ) : return( ZACompound - ZAOther )
    return( ZACompound + ZAOther )

def countZAMasses( info, ZA, AWR ) :

    if( ZA != 0 ) :
        if( ZA not in info.ZA_AWRMasses ) : info.ZA_AWRMasses[ZA] = {}
        ZAData = info.ZA_AWRMasses[ZA]
        if( AWR not in ZAData ) : ZAData[AWR] = 0
        ZAData[AWR] += 1

def printAWR_mode( info, MT, MF, ZA, AWR ) :

    if( info.AWR_mode is not None ) : info.AWR_mode.write( "AWR_mode:: MT: %s: MF: %s: ZA: %s:  AWR: %s::\n" % ( MT, MF, ZA, AWR ) )

def nudgeValue( value, sign ) :

    if( value == 0 ) : raise ValueError( 'FIXME' )
    valueEpsilon = 10**( math.floor( math.log10( abs( value ) ) ) + epsilonExponent )
    return( value + sign * valueEpsilon )

def getCrossSectionForm( info, crossSectionRegions ) :

    axes = crossSectionModule.defaultAxes( )
    if( len( crossSectionRegions ) == 1 ) :         # Store as XYs1d.
        crossSectionForm = crossSectionModule.XYs1d( data = crossSectionRegions[0], label = info.style,
                axes = axes, accuracy = ENDF_Accuracy, interpolation = crossSectionRegions[0].interpolation )
    else :
        crossSectionForm = crossSectionModule.regions1d( label = info.style, axes = axes )
        for region in crossSectionRegions :
            if( len( region ) > 1 ) : crossSectionForm.append( region )
    return( crossSectionForm )

def getMultiplicity( multiplicity, EPrior, Ein ) :

    if(   isinstance( multiplicity, XYsModule.XYs1d ) ) :
        return( multiplicity.evaluate( Ein ) )
    elif( isinstance( multiplicity, regionsModule.regions1d ) ) :
        for region in multiplicity :
            if( Ein <= region.domainMax( ) ) :
                if( region.domainMax( ) == EPrior == Ein ) : continue          # Next region is the one we want.
                return( region.evaluate( Ein ) )
        return( 0 )
    raise Exception( 'unsupported multiplicity form "%s"' % multiplicity.moniker )

def uncorrelated( style, frame, angularSubform, energySubform ) :

    _angularSubform = uncorrelatedModule.angularSubform( angularSubform )
    _energySubform = uncorrelatedModule.energySubform( energySubform )
    return( uncorrelatedModule.form( style, frame, _angularSubform, _energySubform ) )

def getMultiplicityPointwiseOrPieceWise( info, data, warningList ) :

    axes = multiplicityModule.XYs1d.defaultAxes( )
    if( len( data  ) == 1 ) :
        multiplicity = multiplicityModule.XYs1d( data = data[0], label = info.style, axes = axes, 
                accuracy = ENDF_Accuracy, interpolation = data[0].interpolation )
# BRB : fix me.
#    elif( ( len( data ) == 2 ) and ( data[0][-1] == data[1][0] ) and
#            ( data[0].axes[0].interpolation == data[1].axes[0].interpolation ) ) :
#        xys = data[1].copyDataToXYs( )      # This is a XYs1d data masquerading as regions1d, convert to XYs1d.
#        xys[0][0] *= ( 1 + FUDGE_EPS )
#        xys = data[0].copyDataToXYs( ) + xys
#        multiplicity = multiplicityModule.XYs1d( data = xys, axes = axea, accuracy = ENDF_Accuracy )
    else :
        multiplicity = multiplicityModule.regions1d( label = info.style, axes = axes )
        for region in data :
            multiplicity.append( region.copy( axes=axes ) )
    return( multiplicity )

def getTotalOrPromptFission( info, MT456Data, totalOrPrompt, warningList ) :

    ZA, AWR, dummy, LNU, dummy, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MT456Data[0], logFile = info.logs )
    ZA = int( ZA )
    countZAMasses( info, ZA, AWR )
    LNU = int( LNU )
    info.logs.write( '     %s fission neutron data: LNU = %d\n' % ( totalOrPrompt, LNU ) )
    if( LNU == 1 ) :
        dataLine, polynomial = endfFileToGNDMisc.getList( 1, MT456Data, logFile = info.logs )
        axes = multiplicityModule.polynomial.defaultAxes( )
        domainMin, domainMax = 1e-5, 20e6        # BRB, these need to be set from data
        fissionMultiplicity = multiplicityModule.polynomial( coefficients = polynomial['data'], label = info.style, 
                axes = axes, domainMin = domainMin, domainMax = domainMax )
    else :
        axes = multiplicityModule.XYs1d.defaultAxes( )
        dataLine, TAB1, multiplicityRegions = endfFileToGNDMisc.getTAB1Regions( 1, MT456Data, axes = axes, logFile = info.logs )
        fissionMultiplicity = getMultiplicityPointwiseOrPieceWise( info, multiplicityRegions, warningList )
    return( fissionMultiplicity )

def getDelayedFission( info, MT455Data, warningList ) :

    info.logs.write( '     Delayed fission neutron data (MT=455)' )
    MT455DataMF1 = MT455Data[1]
    ZA, AWR, LDG, LNU, dummy, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MT455DataMF1[0], logFile = info.logs )
    ZA = int( ZA )
    countZAMasses( info, ZA, AWR )
    LDG, LNU = int( LDG ), int( LNU )
    info.logs.write( ' LDG=%s LNU=%s' % ( LDG, LNU ) )
    if( LDG != 0 ) : raise Exception( "Only energy-independent delayed fission neutrons are supported" )
    if( LNU != 2 ) : raise Exception( "Only tables of delayed fission neutrons are supported" )

    dataLine, decayRateData = endfFileToGNDMisc.getList( 1, MT455DataMF1, logFile = info.logs )
    NNF = int( decayRateData['NPL'] )
    decayRates = decayRateData['data']

    axes = axesModule.axes( labelsUnits = { 0 : ( 'multiplicty' , '' ), 1 : ( 'energy_in', 'eV' ) } )
    dataLine, TAB1, multiplicityRegions = endfFileToGNDMisc.getTAB1Regions( dataLine, MT455DataMF1, logFile = info.logs, axes = axes )

    interps = [region.interpolation for region in multiplicityRegions]
    if len(set(interps)) > 1:
        raise Exception("Currently only one interpolation flag is supported")
    nubarInterpolation = multiplicityRegions[0].interpolation

    if( 5 in MT455Data ) :
        delayedNeutronEnergies, weights = readMF5( info, 455, MT455Data[5], warningList, delayNeutrons = True )
        weightsSum = XYsModule.XYs1d( [], axes = weights[0].axes, interpolation = weights[0].interpolation )  # Renormalize weights to sum to 1.
        for weight in weights : weightsSum = weightsSum + weight
        for weight in weights :
            for i1, xy in enumerate( weight ) : XYsModule.pointwiseXY.__setitem__( weight, i1, xy[1] / weightsSum.evaluate( xy[0] ) )
        totalDelayedMultiplicity = getMultiplicityPointwiseOrPieceWise( info, multiplicityRegions, warningList )
    else :
        delayedNeutronEnergies = len( decayRates ) * [ None ]
        totalDelayedMultiplicity_ = [ region / len( decayRates ) for region in multiplicityRegions ] 
        nuBar = getMultiplicityPointwiseOrPieceWise( info, totalDelayedMultiplicity_, warningList )
        if( len( decayRates ) > 1 ) : warningList.append( 'More than one delayed fission neutron decay time but no MF = 5 data' )
    if( len( decayRates ) != len( delayedNeutronEnergies ) ) :
        warningList.append( "MF1 MT455 claims %d delayed neutron groups, but MF5 MT455 claims %d groups" % \
                ( len( decayRates ), len( delayedNeutronEnergies ) ) )
        info.doRaise.append( warningList[-1] )
        decayRates = []
    decayChannel = channelsModule.NBodyOutputChannel( )      # Dummy channel, so Q's value does not matter.
    decayChannel.Q.add( toGNDMisc.returnConstantQ( info.style, 0. ) )
    for i1, decayRate in enumerate( decayRates ) :
        energySubform = delayedNeutronEnergies[i1]
        if( energySubform is not None ) :
            weight = weights[i1]
            weightsInterpolation = weight.interpolation
            if( weightsInterpolation == standardsModule.interpolation.flatToken ) :
                if( ( len( weight ) == 2 ) and ( nubarInterpolation == standardsModule.interpolation.linlinToken ) ) :
                    weightsInterpolation = standardsModule.interpolation.linlinToken
                    weight = XYsModule.XYs1d( data = [ xy for xy in weight ], axes = weight.axes )
            if( weightsInterpolation != nubarInterpolation ) :
                raise Exception( 'For total nubar and energy interpolation differ which is currently not supported' )
            if( weightsInterpolation != standardsModule.interpolation.linlinToken ) :
                raise Exception( 'For energy only "lin-lin" interpolation is supported: %s' % weightsInterpolation )
            totalDelayedM = totalDelayedMultiplicity
            if( totalDelayedMultiplicity.domainMax( ) > weight.domainMax( ) ) :
                totalDelayedM = totalDelayedMultiplicity.domainSlice( domainMax = weight.domainMax( ) )
            if( weight.domainMax( ) > totalDelayedM.domainMax( ) ) : weight = weight.domainSlice( domainMax = totalDelayedM.domainMax( ) )
            if( weight.domainMin( ) < totalDelayedM.domainMin( ) ) :
                if( len( weight ) == 2 ) :
                    if( weight[0][1] == weight[1][1] ) : weight = weight.domainSlice( domainMin = totalDelayedM.domainMin( ) )
            multiplicity = totalDelayedM * weight
            nuBar = getMultiplicityPointwiseOrPieceWise( info, [ multiplicity ], warningList )

        # I don't remember why this was added  (round delayed nubar to 10 sig figs). Disabling it for now.
        """
        for xy in nuBar:
            m, e = math.frexp( xy[-1] )
            nuBar.setValue( xy[0], round( xy[-1], 10 - int( e / 3.321928095 ) ) )   # BRB, Caleb, what is this?
            #nuBar.setValue( xy[0], round(xy[-1], 10+int(abs(math.log10(xy[-1])))) )
        """

        particle = toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameGamma( info, 1 ), multiplicity = nuBar )
        particle.addAttribute( 'emissionMode', tokensModule.delayedToken )
        particle.addAttribute( 'decayRate', PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( decayRate ), "1/s" ) )
        if( energySubform is not None ) :
            angularSubform = angularModule.isotropic( )
            form = uncorrelated( info.style, frames[1], angularSubform, energySubform )
            particle.distribution.add( form )
        decayChannel.products.add( decayChannel.products.uniqueLabel( particle ) )
    info.logs.write( '\n' )
    return( decayChannel )

def getFissionEnergies( info, MF458Data, warningList ) :
    """
    For NPLY = 0 this data consists of pairs ( energy, standard deviation ).  The energies are:
    EFR: kinetic energy of all fission fragments
    ENP: energy of prompt fission neutrons
    END: energy of delayed fission neutrons
    EGP: energy of prompt fission gammas
    EGD: energy of delayed fission gammas
    EB: energy of delayed fission betas
    ENU: energy of fission neutrinos
    ER: EFR + ENP + END + EGP + EGD + EB, fission Q minus the neutrinos
    ET: ER + ENU, total fission Q
    For NPLY > 0 the fission energies are polynomials of degree NPLY in the incident energy,
    and the structure of the above table is repeated, once for each polynomial coefficient.
    """

    MF1Data = MF458Data[1]
    dataLine = 0
    ZA, AWR, dummy, dummy, dummy, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF1Data[dataLine], intIndices = [ 0 ], logFile = info.logs )
    ZA = int( ZA )
    countZAMasses( info, ZA, AWR )

    dataLine += 1
    dummy, dummy, dummy, NPLY, lengthData, numEnergies = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF1Data[ dataLine ], logFile = info.logs )
    NPLY = int( NPLY )               # order of the polynomial representation of energy produced
    lengthData = int( lengthData )   # how much data
    numEnergies = int( numEnergies ) # number of fission energies (9)

    dataLine += 1
    energies = endfFileToGNDMisc.nFunkyFloatStringsToFloats( numEnergies, dataLine, MF1Data, dimension = 2, logFile = info.logs )

    labels = fissionEnergyReleasedModule.polynomial.labels
    fissionEnergyData = {}
    for label in labels : fissionEnergyData[label] = []
    for i1 in xrange( 0, numEnergies, 9 ) :
        for j1, label in enumerate( labels ) : fissionEnergyData[label].append( energies[i1 + j1] )

    fissionEnergies = fissionEnergyReleasedModule.component( )
    fissionEnergies.addForm( fissionEnergyReleasedModule.polynomial( info.style, NPLY, fissionEnergyData, 'eV', hasUncertainties = True ) )
## BRB uncomment    warningList.append( 'fission energy release needs work (including processing)' )

    return( fissionEnergies )

def angularLegendrePiecewiseToPointwiseIfPossible( piecewiseForm ) :

    if( len( piecewiseForm ) != 1 ) : return( piecewiseForm )
    pointwise = angularModule.XYs2d( axes = piecewiseForm.axes )
    for series in piecewiseForm[0] : pointwise.append( series )
    return( pointwise )

def angularLegendreToPointwiseOrPiecewiseLegendre( MT, angularData, warningList, MF, msg, subformPointwise = None ) :

    regions = []            # list of regions.
    i1 = 0
    priorEnergy = -1
    lists = angularData['Lists']
    try :
        for i2, interpolationFlag in angularData['interpolationInfo'] :
            if( interpolationFlag > 5 ) : raise Exception( 'Unsupported interpolation = %s for MF=%s, MT=%s' % ( interpolationFlag, MF, MT ) )
            interpolationQualifier, interpolationEin = endfFileToGNDMisc.ENDFInterpolationToGND2plusd( interpolationFlag )
            region = [ interpolationQualifier, interpolationEin, [] ]
            if( i1 > 0 ) :
                if( lists[i1]['C2'] != lastRegion[0] ) : region[2].append( lastRegion )
            for i3 in xrange( i1, i2 ) :
                energy, coefficients = lists[i3]['C2'], [ 1.0 ] + lists[i3]['data']
                if( energy == priorEnergy ) :
                    if( coefficients == lastRegion[-1] ):
                        continue
                    if len(region[2]) > 1: regions.append( region )
                    region = [ interpolationQualifier, interpolationEin, [] ]
                lastRegion = ( energy, coefficients )
                region[2].append( lastRegion )
                priorEnergy = energy
            if len(region[2]) > 1: regions.append( region )
            i1 = i2
    except ValueError as err :
        warningList.append( 'ValueError occurred when constructing LegendreSeries: %s' % err )
        raise

    axes = angularModule.XYs2d.defaultAxes( asLegendre = True )
    if( ( len( regions ) == 1 ) and ( subformPointwise is None ) ) :
        interpolationQualifier, interpolationEin, region = regions[0]
        subformLegendre = angularModule.XYs2d( axes = axes, interpolation = interpolationEin,
                interpolationQualifier = interpolationQualifier )
        for i1, ( energy, coefficients ) in enumerate( region ) :
            subformLegendre.append( angularModule.Legendre( coefficients = coefficients, value = energy ) )
    else :
        subformLegendre = angularModule.regions2d( axes = axes )
        for i1, regionInfo in enumerate( regions ) :
            interpolationQualifier, interpolationEin, region = regionInfo
            LegendreRegion = angularModule.XYs2d( interpolation = interpolationEin,
                    interpolationQualifier = interpolationQualifier )
            for i2, ( energy, coefficients ) in enumerate( region ) :
                LegendreRegion.append( angularModule.Legendre( coefficients = coefficients, value = energy ), makeCopy = False )
            subformLegendre.append( LegendreRegion )
        if( subformPointwise is not None ) :
            if( isinstance( subformPointwise, angularModule.regions2d ) ) :
                    raise 'hell - need to implement this'
            else :
                region = angularModule.XYs2d( interpolation = subformPointwise.interpolation,
                    interpolationQualifier = subformPointwise.interpolationQualifier, axes = subformPointwise.axes )
                for data in subformPointwise : region.append( data )
                subformLegendre.append( region )
    return( subformLegendre )

def convertNuclearPlusInterferenceDataToPiecewise( MT, angularData, warningList, MF, msg, identicalParticles ) :
    """Return distribution.LegendreNuclearPlusCoulombInterference instance. This in turn contains
    Legendre sections for both the nuclear and interference contributions."""

    axes = angularModule.XYs2d.defaultAxes( asLegendre = True )
    nuclear = angularModule.regions2d( axes = axes )
    interferenceReal = angularModule.regions2d( axes = axes )
    interferenceImaginary = angularModule.regions2d( axes = axes )
    index, start, lastRegion = 0, 0, None
    lists = angularData['Lists']
    for end, interpolationFlag in angularData['interpolationInfo'] :
        interpolationQualifier, interpolationE_in = endfFileToGNDMisc.ENDFInterpolationToGND2plusd( interpolationFlag )
        if( interpolationFlag > 6 ) : raise Exception( 'Unsupported interpolation = %s for MF=%s, MT=%s' % ( interpolationFlag, MF, MT ) )
        #interpolationE_in, interpolationCl, interpolationQualifier = endfFileToGNDMisc.ENDFInterpolationToGNDAxes3plusd( interpolationFlag )
        region_Nuc = angularModule.XYs2d( interpolation = interpolationE_in )
        region_IntReal = angularModule.XYs2d( interpolation = interpolationE_in )
        region_IntImaginary = angularModule.XYs2d( interpolation = interpolationE_in )
        priorEnergy = -1
        if( lastRegion is not None ) :
            if( lists[start]['C2'] != lastRegion ) : # ensure no gaps between regions
                region_Nuc.append( angularModule.Legendre( coefficients = nuclear_term,
                        value = lastRegion ), makeCopy = False )
                region_IntReal.append( angularModule.Legendre( coefficients = interference_termReal,
                        value = lastRegion ), makeCopy = False )
                region_IntImaginary.append( angularModule.Legendre( coefficients = interference_termImaginary,
                        value = lastRegion ), makeCopy = False )
        for idx in xrange( start, end ) :
            list = lists[idx]
            energy = list['C2']
            if( energy == priorEnergy ) :                           # This fixes a problem with some data having two same energy values.
                energy += FUDGE_EPS * energy
                warningList.append( 'same energies, second one being incremented for MT = %d, MF = %d, %s' % ( MT, MF, msg ) )
            priorEnergy = energy
            if( identicalParticles ) :
                splitPoint = list['N2'] + 1
            else :
                splitPoint = list['N2'] * 2 + 1
            nuclear_term = list['data'][:splitPoint]
            interference_term, interference_termReal, interference_termImaginary = list['data'][splitPoint:], [], []
            for jdx in range( 0, len( interference_term ), 2 ) :
                interference_termReal.append( interference_term[jdx] )
                interference_termImaginary.append( interference_term[jdx+1] )
            region_Nuc.append( angularModule.Legendre( coefficients = nuclear_term,
                    value = energy ), makeCopy = False )
            region_IntReal.append( angularModule.Legendre( coefficients = interference_termReal,
                    value = energy ), makeCopy = False )
            region_IntImaginary.append( angularModule.Legendre( coefficients = interference_termImaginary,
                    value = energy ), makeCopy = False )
        lastRegion = energy
        nuclear[index] = region_Nuc
        interferenceReal[index] = region_IntReal
        interferenceImaginary[index] = region_IntImaginary
        index += 1
        start = end
    nuclear = angularLegendrePiecewiseToPointwiseIfPossible( nuclear )
    interferenceReal = angularLegendrePiecewiseToPointwiseIfPossible( interferenceReal )
    interferenceImaginary = angularLegendrePiecewiseToPointwiseIfPossible( interferenceImaginary )
    return( nuclear, interferenceReal, interferenceImaginary )

def convertAngularToPointwiseOrPiecewiseFromTAB2_TAB1( MT, angularTAB1, warningList ) :

    angularData = angularTAB1['TAB2s']
    axes = angularModule.XYs2d.defaultAxes( )
    if( len( angularData ) == 1 ) :
        interplation, angularData = angularData[0]      # BRB: need to check interpolation?
        subform = angularModule.XYs2d( axes = axes )
        for xys in angularData :
            if( len( xys ) > 1 ) : raise 'hell - need to support this'
            xys = xys[0]
            xys = angularModule.XYs1d( data = xys, interpolation = xys.interpolation,
                accuracy = ENDF_Accuracy, value = xys.value )
            subform.append( xys )
    else :
        raise 'hell - regions2d tabulated angular not currently supported'
    return( subform )

def convertAngularToPointwiseOrPiecewiseFromTAB2_List( MT, LANG, angularList, warningList ) :
    """
    Like convertAngularToPointwiseOrPiecewiseFromTAB2_TAB1 except mu,P given as LISTs instead of TAB1s.
    """

    angularData = angularList['Lists']
    axes = angularModule.XYs2d.defaultAxes( )
    try :
        interpolation = { 12 : 2, 14 : 4 }[LANG]
    except :
        print 'interpolation = %d' % interpolation
        raise 'hell - what is this and how does it differ from 12'
    interpolation = endfFileToGNDMisc.ENDFInterpolationToGND1d( interpolation )
    subform = angularModule.XYs2d( axes = axes )
    e_in_Prior = -1
    for i1, list in enumerate( angularData ) :
        LANG_e_in = int( list['L1'] )
        if( LANG_e_in != LANG_e_in ) : raise 'hell - need to implement this'
        e_in = list['C2']
        if( e_in == e_in_Prior ) : raise 'hell - need to implement this'
        data = list['data']
        if( MT == 526 ) :        # Some MT 526 have data like '0.999998+0 2.046010+5 0.999998+0 2.870580+5' which this fixes.
            j1, j2 = None, None
            for j3 in range( 0, len( data ), 2 ) :
                if( j2 is not None ) :
                    if( data[j2] == data[j3] ) : data[j2] -= 0.5 * ( data[j2] - data[j1] )
                j1, j2 = j2, j3
        xys = angularModule.XYs1d( data = list['data'], dataForm = 'list', 
                interpolation = interpolation, accuracy = ENDF_Accuracy, value = e_in )
        subform.append( xys )
        e_in_Prior = e_in
    return( subform )

def toPointwiseOrPiecewiseEnergy( MT, TAB2 ) :

    def getEpP( energyPrior, data ) :

        EpP = data[0]
        energy = float( EpP.value )
        if( len( data ) == 1 ) :
            EpP = energyModule.XYs1d( data = EpP, interpolation = EpP.interpolation, value = energy )
        else :
            EpP = energyModule.regions1d( value = energy )
            for datum in data :
                EpP.append( energyModule.XYs1d( data = datum, interpolation = datum.interpolation ) )
        return( energy, EpP )

    axes = energyModule.XYs2d.defaultAxes( )
    if( TAB2['NR'] == 1 ) :
        interpolation, data = TAB2['TAB2s'][0]
        interpolationQualifier, interpolation = endfFileToGNDMisc.ENDFInterpolationToGND2plusd( interpolation )
        subform = energyModule.regions2d( axes = axes )
        subformRegion = energyModule.XYs2d( axes = axes, interpolation = interpolation, interpolationQualifier = interpolationQualifier )
        energyPrior = -1
        for EpP in data :
            EIn = float( EpP[0].value )
            if( EIn == energyPrior ) :
                if( len( subformRegion ) > 1 ) : subform.append( subformRegion )
                subformRegion = energyModule.XYs2d( axes = axes, interpolation = interpolation, interpolationQualifier = interpolationQualifier )
                energyPrior = -1
            energyPrior, EpPp = getEpP( energyPrior, EpP )
            subformRegion.append( EpPp, makeCopy = False )
        if( len( subform ) == 0 ) :
            subform = subformRegion
        else :
            subform.append( subformRegion )
    else :
        subform = energyModule.regions2d( axes = axes )
        TAB2s = TAB2['TAB2s']
        for i1, ( interpolation, TAB1s ) in enumerate( TAB2s ) :
            interpolationQualifier, interpolation = endfFileToGNDMisc.ENDFInterpolationToGND2plusd( interpolation )
            region = energyModule.XYs2d( interpolation = interpolation, interpolationQualifier = interpolationQualifier )
            energyPrior = -1
            for EpP in TAB1s :
                energyPrior, EpPp = getEpP( energyPrior, EpP )
                region.append( EpPp, makeCopy = False )
            subform.append( region )

    return( subform )

def readMF2( info, MF2, warningList ) :
    """
    parse MF2 into resonances class (and sub-classes)
    """
    from xData import table as tableModule

    # store MT #s for all reactions that need to include resonance data:
    resonanceMTs = set()
    # check for inconsistent spins
    targetSpins = set()

    def readResonanceSection( LRU, LRF, NRO, NAPS ):
        """Helper function, read in resonance info for one energy range."""

        if NRO!=0:  # energy-dependent scattering radius
            if NAPS==2: raise BadResonances("NAPS=2 option not yet supported!")
            line1 = mf2.next()
            dum, dum, dum, dum, NR, NP = funkyFI( line1, logFile = info.logs )
            nLines = NR//3 + bool(NR%3)  +  NP//3 + bool(NP%3)
            data = [line1] + [mf2.next() for i in range(nLines)]
            axes = axesModule.axes( labelsUnits = { 1 : ( 'energy_in', 'eV' ), 0 : ( 'radius', '10*fm' ) } )
            dataLine, TAB1, regions = endfFileToGNDMisc.getTAB1Regions( 0, data, axes = axes, logFile = info.logs )
            if TAB1['NR']!=1:
                raise NotImplementedError("multi-region scattering radius")
            data = regions[0]
            data = data.convertAxisToUnit(0,'fm')
            scatRadius = resonancesModule.scatteringRadius( data )

        if LRU==0:  # scattering radius only. Note AP given in 10*fm
            SPI, AP, dum, dum, NLS, dum = funkyFI( mf2.next(), logFile = info.logs )
            targetSpins.add( SPI )
            info.particleSpins['target'] = ( xParticleModule.spin(SPI), 0 ) # no parity information
            value = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( AP*10 ), "fm" )
            scatRadius = resonancesModule.scatteringRadius(
                        resonancesModule.constantScatteringRadius( value, bounds = ( EL, EH ), label = info.style ) )
            return scatRadius

        elif LRU==1 and LRF in (1,2):   #SLBW or MLBW form
            SPI, AP, dum, dum, NLS, dum = funkyFI( mf2.next(), logFile = info.logs )
            targetSpins.add( SPI )
            info.particleSpins['target'] = ( xParticleModule.spin(SPI), 0 )
            if NRO==0:
                scatRadius = resonancesModule.scatteringRadius(
                        resonancesModule.constantScatteringRadius( PQUModule.PQU(
                            PQUModule.pqu_float.surmiseSignificantDigits( AP*10 ),"fm") ) )
            resList = []
            negativeJs = False
            for lidx in range(NLS):
                AWRI, QX, L, LRX, tmp, NRS = funkyFI( mf2.next(), logFile = info.logs )
                if tmp!=6*NRS:
                    raise BadResonances( "incorrect values in resonance section line %d" % mf2.index )
                for line in range(NRS):
                    e,j,gtot,gn,gg,gf = funkyF( mf2.next(), logFile = info.logs )
                    if j<0: negativeJs = True
                    resList.append( [e,L,j,gtot,gn,gg,gf] )
            if negativeJs: raise BadResonances("Encountered negative J-values for SLBW/MLBW")

            table = tableModule.table(
                columns = [
                tableModule.columnHeader( 0, name="energy", units="eV" ),
                tableModule.columnHeader( 1, name="L", units="" ),
                tableModule.columnHeader( 2, name="J", units="" ),
                tableModule.columnHeader( 3, name="totalWidth", units="eV" ),
                tableModule.columnHeader( 4, name="neutronWidth", units="eV" ),
                tableModule.columnHeader( 5, name="captureWidth", units="eV" ),
                tableModule.columnHeader( 6, name="fissionWidthA", units="eV" ), ],
                data = sorted(resList, key=lambda(res): res[0]) )   # sort by energy only
            for column in ("totalWidth","fissionWidthA"):
                if not any( table.getColumn(column) ):
                    table.removeColumn(column)

            if LRF==1:
                return resonancesModule.SLBW( resonancesModule.resonanceParameters(table), scatteringRadius=scatRadius,
                        calculateChannelRadius=not(NAPS) )
            else:
                return resonancesModule.MLBW( resonancesModule.resonanceParameters(table), scatteringRadius=scatRadius,
                        calculateChannelRadius=not(NAPS) )

        elif LRU==1 and LRF==3:     # Reich-Moore form
            SPI, AP, LAD, dum, NLS, NLSC = funkyFI( mf2.next(), logFile = info.logs )
            targetSpins.add( SPI )
            info.particleSpins['target'] = ( xParticleModule.spin(SPI), 0 ) # store spin in GND particle list
            assert NRO==0
            defaultScatteringRadius = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( AP*10 ),"fm" )
            resList = []
            LdependentAP = {}
            negativeJs = False  # in Reich-Moore, j<0 indicates channel spin
            for lidx in range(NLS):
                AWRI, APL, L, dum, tmp, NRS = funkyFI( mf2.next(), logFile = info.logs )
                if tmp!=6*NRS:
                    raise BadResonances("incorrect values in resonance section line %d" % mf2.index)
                if APL:
                    LdependentAP[L] = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( APL*10 ),"fm")
                for line in range(NRS):
                    e,j,gn,gg,gf1,gf2 = funkyF( mf2.next(), logFile = info.logs )
                    resList.append( [e,L,j,gn,gg,gf1,gf2] )
                    if j<0: negativeJs = True
            if LdependentAP:
                scatRadiusForm = resonancesModule.LdependentScatteringRadii( defaultScatteringRadius, LdependentAP,
                                                                           label = info.style )
            else:
                scatRadiusForm = resonancesModule.constantScatteringRadius( defaultScatteringRadius, label = info.style )
            scatRadius = resonancesModule.scatteringRadius( scatRadiusForm )

            columns = [
                tableModule.columnHeader( 0, name="energy", units="eV" ),
                tableModule.columnHeader( 1, name="L", units="" ),
                tableModule.columnHeader( 2, name="J", units="" ),
                tableModule.columnHeader( 3, name="neutronWidth", units="eV" ),
                tableModule.columnHeader( 4, name="captureWidth", units="eV" ),
                tableModule.columnHeader( 5, name="fissionWidthA", units="eV" ),
                tableModule.columnHeader( 6, name="fissionWidthB", units="eV" ), ]

            if negativeJs:
                # for incident neutrons, channel spin s = vector sum(target spin, 1/2)
                # in Reich-Moore, the channel spin is encoded in the sign of J
                e,L,j,gn,gg,gf1,gf2 = zip(*resList)
                channelSpins = [SPI + math.copysign(0.5,jval) for jval in j]
                jlist = [abs(jval) for jval in j]
                resList = [list(row) for row in zip( e,L,jlist,channelSpins,gn,gg,gf1,gf2 )]
                columns.insert(3, tableModule.columnHeader( 3, name="channelSpin", units="" ) )
                for i in range(len(columns)): columns[i].index = i

            table = tableModule.table( columns = columns,
                    data = sorted(resList, key=lambda(res): res[0]) ) # sort by resonance energy
            for column in ('fissionWidthA','fissionWidthB'):
                if not any( table.getColumn(column) ):
                    table.removeColumn(column)

            return resonancesModule.RM( resonancesModule.resonanceParameters(table), scatteringRadius=scatRadius,
                    calculateChannelRadius=not(NAPS), computeAngularDistribution=LAD,
                    LvaluesNeededForConvergence=NLSC )

        elif LRU==1 and LRF==4:     # Adler-Adler, not currently supported
            raise BadResonances( "Adler-Adler resonance formalism not yet supported!" )

        elif LRU==1 and LRF==7:     # R-Matrix Limited
            dum,dum,IFG,KRM,NJS,KRL = funkyFI( mf2.next(), logFile = info.logs )
            if KRM==3:
                approximation = 'Reich_Moore'
            elif KRM==4:
                approximation = 'Full R-Matrix'
            else:
                raise BadResonances( "R-Matrix with KRM=%d not yet implemented!\n" % KRM )

            dum,dum,NPP,dum,tmp1,tmp2 = funkyFI( mf2.next(), logFile = info.logs )
            if tmp1!=12*NPP or tmp2!=2*NPP:
                raise BadResonances( "incorrect LRF7 header!" )

            # some helper functions:
            def getOutgoingParticles( MT, targZA, projZA ):
                reacStr = endf_endl.ENDF_MTZAEquation(projZA,targZA, MT)[1]
                outgoing = reacStr.split('->')[1].strip()
                pA, pB = outgoing.split()[::2]
                return pA, pB

            def translateENDFJpi(I,P):
                # endf uses weird convention for Jpi. Translate to simpler version:
                spin = abs(I)
                if I: parity = abs(I)/I
                else: parity = P or 1   # if (I,P)=(0,0) treat as '0+'
                return xParticleModule.spin( spin ), xParticleModule.parity( parity )

            # ENDF R-Matrix starts by listing outgoing particle pairs
            # these are referred back to later on. Store in 'IPPlist':
            openChannels = resonancesModule.openChannels()
            for idx in range(NPP):
                MA, MB, ZA, ZB, IA, IB = funkyF( mf2.next(), logFile = info.logs )
                Q, PNT, SHF, MT, PA, PB = funkyF( mf2.next(), logFile = info.logs )
                MT = int(MT)
                resonanceMTs.add(MT)

                Qvalue = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( Q ),'eV')
                calculatePenetrability = False
                if PNT==1 or (PNT==0 and MT in (19,102)):
                    calculatePenetrability = True

                # identify the channel using ZA and MT:
                projectileZA = info.reactionSuite.projectile.getZ_A_SuffixAndZA()[-1]
                pA,pB = getOutgoingParticles(MT,int(ZAM),projectileZA)
                # get target spin. In future, this should already be present in particle list
                info.particleSpins[pA] = translateENDFJpi(IA,PA)
                info.particleSpins[pB] = translateENDFJpi(IB,PB)

                channelName = "%s + %s" % (pA,pB)
                gndChannel, = [chan for chan in info.reactionSuite.reactions if chan.ENDF_MT == MT]
                reactionLink = linkModule.link(gndChannel)
                openChannels.add(
                    resonancesModule.openChannel( label=str(idx), reactionLink=reactionLink, name=channelName, ENDF_MT=MT,
                    Q=Qvalue, calculatePenetrability=calculatePenetrability, calculateShift=bool(SHF) ) )


            # next we have NJS spin groups, each containing channels and resonances
            spinGroups = []
            radii = {}
            for spinGroupIndex in range(NJS):
                # read scattering radius, binding, etc:
                AJ, PJ, KBK, KPS, tmp, NCH = funkyFI( mf2.next(), logFile = info.logs )
                if tmp!=6*NCH:
                    raise BadResonances("incorrect LRF7 header, line %d" % mf2.index)
                idx = 0
                channelData = [ tableModule.columnHeader(idx, name="energy", units="eV") ]
                channelNames = []
                for chidx in range(NCH):
                    idx += 1
                    IPP, L, SCH, BND, APE, APT = funkyF( mf2.next(), logFile = info.logs )
                    thisChannel = openChannels[str(int(IPP)-1)]
                    channelName = "%s width" % thisChannel.name
                    jdx = 2
                    while True:
                        if channelName not in channelNames:
                            channelNames.append( channelName ); break
                        channelName = '%s width_%d' % (thisChannel.name, jdx)
                        jdx += 1

                    radii.setdefault(thisChannel, []).append( (spinGroupIndex,
                        PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( APT*10 ),"fm"),
                        PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( APE*10 ),"fm") ) )

                    channelHeaders = tableModule.columnHeader(idx, name=channelName, units="eV",
                        L=int(L), channelSpin=SCH, boundaryCondition=BND)
                    channelData.append( channelHeaders )

                # resonances for this J:
                dum,dum,dum,NRS,tmp,NX = funkyFI( mf2.next(), logFile = info.logs )
                if tmp!=6*NX:
                    raise BadResonances("incorrect LRF7 header, line %d" % mf2.index)
                if NRS==0: mf2.next()   # skip empty line
                resonances = []
                for i in range(NRS):
                    nlines = int(math.ceil( (NCH+1)/6.0 )) # Extra "1" is for the Eres column
                    vals = []
                    for j in range(nlines):
                        vals += funkyF( mf2.next(), logFile = info.logs )
                    resonances.append( vals[:NCH+1] )
                table = tableModule.table( columns=channelData, data=resonances )
                # done with this spin group:
                J, pi = translateENDFJpi(AJ,PJ)
                spinGroups.append( resonancesModule.spinGroup(spinGroupIndex, J, pi, KBK, KPS,
                        resonancesModule.resonanceParameters(table) ) )

            # make one default value for true radius, effective radius and boundary condition.
            # Only include them in specific channels if we need to override the default value.

            from collections import Counter
            for channel_ in openChannels:
                sgIndices, trueRad, effRad = zip(*radii[channel_])
                if not any(trueRad) and not any(effRad): continue     # ignore radii for (n,gamma)

                APTmostCommon = Counter(trueRad).most_common(1)[0][0]
                channel_.scatteringRadius = resonancesModule.scatteringRadius(
                    resonancesModule.constantScatteringRadius( APTmostCommon ) )
                APEmostCommon = Counter(effRad).most_common(1)[0][0]
                if APEmostCommon != APTmostCommon:
                    channel_.effectiveRadius = resonancesModule.effectiveRadius(
                        resonancesModule.constantScatteringRadius( APEmostCommon ) )

                for (idx,val1,val2) in radii[channel_]:
                    override = resonancesModule.channelOverride( channel_.label )
                    if val1 != APTmostCommon:
                        override.scatteringRadius = resonancesModule.scatteringRadius(
                            resonancesModule.constantScatteringRadius( val1 ) )
                    if val2 != APEmostCommon:
                        override.effectiveRadius = resonancesModule.effectiveRadius(
                            resonancesModule.constantScatteringRadius( val2 ) )

                    if override:
                        spinGroups[idx].overrides.add( override )

            kwargs = {}
            chanList = [chan for sg in spinGroups for chan in sg.resonanceParameters.table.columns]
            for option in ('boundaryCondition',):
                vals = [chan.attributes.get(option, None) for chan in chanList]
                valSet = [v for v in set(vals) if v is not None]
                if not valSet: continue
                counts = [vals.count(v) for v in valSet]
                mostCommon = valSet[ counts.index( max(counts) ) ]
                kwargs[option] = mostCommon
                for chan in chanList:
                    opt = chan.attributes.get(option,None)
                    if opt == mostCommon or not opt: chan.attributes[ option ] = None

            # also fix-up the penetrability and phase shift: define only once at the top if possible:
            for option in ('calculatePenetrability','calculateShift'):
                vals = [getattr(chan, option, None) for chan in openChannels]
                valSet = [v for v in set(vals)]
                counts = [vals.count(v) for v in valSet]
                mostCommon = valSet[ counts.index( max(counts) ) ]
                kwargs[option] = mostCommon
                for chan in openChannels:
                    opt = getattr(chan,option,None)
                    if opt == mostCommon: setattr(chan, option, None)

            # end of spin groups. write RMatrix class:
            #spinGroups.sort()   # sort by Jpi. Disable for easier comparison with ENDF-6
            return resonancesModule.RMatrix( openChannels, spinGroups, approximation=approximation,
                    relativisticKinematics=bool(KRL), reducedWidthAmplitudes=bool(IFG),
                    calculateChannelRadius=not(NAPS), **kwargs )

        elif LRU==2: # unresolved
            L_list = []
            LRF_LFW = "LRF,LFW=%d,%d" % (LRF, LFW)

            SPI,AP,LSSF,dum,NE,NLS = funkyFI( mf2.next(), logFile = info.logs )
            targetSpins.add( SPI )
            if info.target.name not in info.particleSpins:
                info.particleSpins[ info.target.name ] = ( xParticleModule.spin(SPI), 0 )
            if NRO==0:
                scatRadius = resonancesModule.scatteringRadius(
                        resonancesModule.constantScatteringRadius( PQUModule.PQU(
                            PQUModule.pqu_float.surmiseSignificantDigits( AP*10 ),"fm"), label = info.style ) )

            if LFW==0 and LRF==1:   # 'Case A', see ENDF 2010 manual page 70
                NLS = NE
                for lidx in range(NLS):
                    AWRI, dum, L, dum, tmp, NJS = funkyFI( mf2.next(), logFile = info.logs )
                    if tmp!=6*NJS:
                        raise BadResonances("bad unresolved flag, line %d" % mf2.index)
                    J_list = []
                    for jidx in range(NJS):
                        D,AJ,AMUN,GNO,GG,dum = funkyF( mf2.next(), logFile = info.logs )
                        if AMUN.is_integer(): AMUN = int(AMUN)
                        D,GNO,GG,GF = [PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( a ),'eV') for a in (D,GNO,GG,0)]
                        J_list.append( resonancesModule.URR_Jsection( xParticleModule.spin(AJ), eDepWidths=[],
                            constantWidths={'levelSpacing':D,'neutronWidth':GNO,'captureWidth':GG, 'fissionWidthA':GF,
                                            'competitiveWidth':PQUModule.PQU(0,'eV')},
                            neutronDOF=AMUN, gammaDOF=0, competitiveDOF=0, fissionDOF=0 ) )
                    L_list.append( resonancesModule.URR_Lsection( L, J_list ) )
                return resonancesModule.unresolvedTabulatedWidths( L_list,
                        interpolation = standardsModule.interpolation.linlinToken, 
                        scatteringRadius=scatRadius, forSelfShieldingOnly=bool(LSSF),
                        ENDFconversionFlag=LRF_LFW)

            elif LFW==1 and LRF==1: # 'Case B'
                nlines = int(math.ceil(NE/6.0))
                energyList = []
                for i in range(nlines):
                    #energyList += [PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( a ),'eV') for a in funkyF(mf2.next(), logFile = info.logs )]
                    energyList += funkyF(mf2.next(), logFile = info.logs)
                for lidx in range(NLS):
                    AWRI,dum,L,dum,NJS,dum = funkyFI( mf2.next(), logFile = info.logs )
                    J_list = []
                    for jidx in range(NJS):
                        eDepWidths = []
                        dum,dum,L,MUF,tmp,dum = funkyFI( mf2.next(), logFile = info.logs )
                        if tmp!=NE+6:
                            raise BadResonances("Bad unresolved flag, line %d" % mf2.index)
                        D,AJ,AMUN,GNO,GG,dum = funkyF( mf2.next(), logFile = info.logs )
                        if AMUN.is_integer(): AMUN = int(AMUN)
                        D,GNO,GG,GX = [PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( a ),'eV') for a in (D,GNO,GG,0)]
                        widthList = []
                        for i in range(nlines):
                            #widthList += [PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( a ),'eV') for a in funkyF( mf2.next(), logFile = info.logs )]
                            widthList += funkyF( mf2.next(), logFile = info.logs )
                        table = tableModule.table( [
                            tableModule.columnHeader( 0, name="energy", units="eV" ),
                            tableModule.columnHeader( 1, name="fissionWidthA", units="eV" ) ] )
                        for e in range(NE):
                            table.addRow( [energyList[e], widthList[e]] )
                            #eDepWidths.append({'energy':energyList[e], 'fissionWidthA':widthList[e]})
                        J_list.append( resonancesModule.URR_Jsection( xParticleModule.spin(AJ), table,
                            constantWidths = { 'levelSpacing' : D, 'neutronWidth' : GNO, 'captureWidth' : GG, 'competitiveWidth' : GX },
                            neutronDOF = AMUN, fissionDOF = MUF, competitiveDOF = 0 ) )
                    L_list.append( resonancesModule.URR_Lsection( L, J_list ) )
                return resonancesModule.unresolvedTabulatedWidths( L_list,
                        interpolation = standardsModule.interpolation.linlinToken, 
                        scatteringRadius=scatRadius, forSelfShieldingOnly=bool(LSSF),
                        ENDFconversionFlag=LRF_LFW)

            elif LRF==2:            # 'Case C', most common in ENDF-VII.1
                NLS = NE
                interpolations = []
                for Lidx in range(NLS):
                    J_list = []
                    AWRI,dum,L,dum,NJS,dum = funkyFI( mf2.next(), logFile = info.logs )
                    for jidx in range(NJS):
                        resList = []
                        AJ,dum,INT,dum,tmp,NE = funkyFI( mf2.next(), logFile = info.logs )
                        interpolations.append(INT)
                        if tmp!=6*NE+6:
                            raise BadResonances("bad unresolved flag, line %d" % mf2.index)
                        dof = {}
                        dum,dum,dof['AMUX'],dof['AMUN'],dof['AMUG'],dof['AMUF'] = funkyF( mf2.next(), logFile = info.logs )
                        for key in ('AMUX','AMUN','AMUG','AMUF'):
                            if type(dof[key]) is float and dof[key].is_integer(): dof[key] = int(dof[key])
                        for i in range(NE):
                            resList.append( funkyF( mf2.next(), logFile = info.logs ) )
                        table = tableModule.table( columns= [
                            tableModule.columnHeader( 0, name="energy", units="eV" ),
                            tableModule.columnHeader( 1, name="levelSpacing", units="eV" ),
                            tableModule.columnHeader( 2, name="competitiveWidth", units="eV" ),
                            tableModule.columnHeader( 3, name="neutronWidth", units="eV" ),
                            tableModule.columnHeader( 4, name="captureWidth", units="eV" ),
                            tableModule.columnHeader( 5, name="fissionWidthA", units="eV" ), ],
                            data = resList )
                        J_list.append( resonancesModule.URR_Jsection( xParticleModule.spin(AJ), table,
                                constantWidths={}, competitiveDOF=dof['AMUX'], neutronDOF=dof['AMUN'],
                                gammaDOF=dof['AMUG'], fissionDOF=dof['AMUF'] ) )
                        # sometimes energy-dependent flag is used, but widths are constant in energy:
                        J_list[-1].eliminateRedundantInfo()
                    L_list.append( resonancesModule.URR_Lsection( L, J_list ) )
                if len(set(interpolations))>1:
                    warningList.append( 'inconsistent interpolations in unresolved region will be ignored!' )
                return resonancesModule.unresolvedTabulatedWidths( L_list,
                        interpolation = endfFileToGNDMisc.ENDFInterpolationToGND1d( interpolations[0] ),
                        scatteringRadius = scatRadius, forSelfShieldingOnly = bool( LSSF ),
                        ENDFconversionFlag = LRF_LFW )

        else:
            info.logs.write( "Unexpected LRU=%d, LRF=%d encountered\n" % ( LRU, LRF ) )

    # end of helper functions.
    # now read MF2 data:
    mf2 = myIter(MF2) # mf2.next() to get each line

    scatteringRadius = None
    resolvedList = []
    unresolvedList = []

    # store particle spins, will be stored in the particle list:
    info.particleSpins = {}

    # read MF2 header:
    ZAM, AWR, dum, dum, NIS, dum = funkyFI( mf2.next(), logFile = info.logs )
    ZAM = int( ZAM )
    printAWR_mode( info, -1, 2, ZAM, AWR )
    countZAMasses( info, ZAM, AWR )
    if NIS!=1: info.logs.write( "careful, more than one isotope in MF2!" )
    ZAI, ABN, dum, LFW, NER, dum = funkyFI( mf2.next(), logFile = info.logs )

    for erange in range(NER):
        # each energy range
        EL, EH, LRU, LRF, NRO, NAPS = funkyFI( mf2.next(), logFile = info.logs )
        EL, EH = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( EL ),"eV"), PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( EH ),"eV")
        resonanceSection = readResonanceSection( LRU, LRF, NRO, NAPS )
        if LRU==0:
            scatteringRadius = resonanceSection
        else:
            resonanceMTs.update( [1,2,3,18,19,102] )
        if LRU==1:
            resolvedList.append( (resonanceSection,EL,EH) )
        elif LRU==2:
            unresolvedList.append( (resonanceSection,EL,EH) )
    if not resolvedList: resolved = None
    elif len(resolvedList)==1:
        resolved = resonancesModule.resolved( *resolvedList[0] )
    else:
        warningList.append( "multiple resolved energy intervals are deprecated!" )
        resolved = resonancesModule.resolved( multipleRegions=True )
        idx = 0
        for resonanceSection, EL, EH in resolvedList:
            resolved.regions.append( resonancesModule.energyInterval(idx,resonanceSection,EL,EH) )
            idx += 1
    if not unresolvedList: unresolved = None
    elif len(unresolvedList)==1:
        unresolved = resonancesModule.unresolved( *unresolvedList[0] )
    else:
        raise BadResonances( "multiple unresolved regions not supported" )

#    if mf2.index != mf2.length: raise BadResonances("Not all resonance data converted!")
    if mf2.index != mf2.length: warningList.append("Not all resonance data converted!")
    if len( targetSpins ) > 1: warningList.append("multiple values for nuclear spin encountered: %s" % targetSpins)
    resonances = resonancesModule.resonances( scatteringRadius, resolved, unresolved )
    return resonances, sorted(resonanceMTs)


def readMF3( info, MT, MF3Data, warningList ) :

    axes = crossSectionModule.defaultAxes( )
    dataLine, TAB1, crossSectionRegions = endfFileToGNDMisc.getTAB1Regions( 1, MF3Data, allowInterpolation6 = True, 
            logFile = info.logs, axes = axes, cls = crossSectionModule.XYs1d )
    QM, QI, LR = TAB1['C1'], TAB1['C2'], int( TAB1['L2'] )
    breakupProducts = None
    if(   LR == 0 ) :
        pass
    elif( LR in [ 22, 23, 24, 25, 28, 29, 30, 32, 33, 34, 35, 36 ] ) :
        info.logs.write( ' : MF=3, LR=%s' % LR )
        breakupProducts, productCounts = {}, endf_endl.endfMTtoC_ProductLists[LR].productCounts
        for product in productCounts :
            if( productCounts[product] != 0 ) : breakupProducts[product] = productCounts[product]
        breakupProducts[info.projectile] -= 1
        if( breakupProducts[info.projectile] == 0 ) : del breakupProducts[info.projectile]
    elif( LR == 31 ) :
        warningList.append( 'Invalid LR = %s for MT = %s is being ignored' % ( LR, MT ) )
    elif( LR in [ 1, 39, 40 ] ) :
        if( LR == 40 ) :
            warningList.append( 'LR = %s for MT = %s is being ignored' % ( LR, MT ) )
        else :
            if ( MT != 5 ):
                warningList.append( "Breakup LR = %s is not supported: MT = %s" % ( LR, MT ) )
                #raise NotImplementedError( "Breakup LR = %s is not supported: MT = %s" % ( LR, MT ) )
    else :
        raise Exception( "Invalid breakup flag LR %s: MT = %d" % ( LR, MT ) )

    crossSection = getCrossSectionForm( info, crossSectionRegions )

    return( QM, QI, crossSection, breakupProducts )

def readMF4( info, product, MT, MF4Data, formClass, warningList ) :

    ZA, AWR, LVT, LTT, dummy, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF4Data[0], logFile = info.logs )
    ZA = int( ZA )
    printAWR_mode( info, MT, 4, ZA, AWR )
    countZAMasses( info, ZA, AWR )
    LVT = int( LVT )                # 1: transformation matrix given. Must be 0 for endf/b6 format but not older formats.
    LTT = int( LTT )                # 0: isotropic, 1: Legendre, 2: table, 3: Legendre for low E and table for high E.

    dummy, AWR_, LI, LCT, NK, NM = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF4Data[1], logFile = info.logs )
    if( AWR != AWR_ ) : printAWR_mode( info, MT, -4, ZA, AWR )
    LI = int( LI )                  # if 1, gammas isotropic
    LCT = int( LCT )                # 1 for lab frame, 2 for center of mass
    NK = int( NK )                  # number of entries in transformation matrix
    NM = int( NM )                  # maximum Legendre order
    if( ( LCT != 2 ) and ( formClass == angularModule.twoBodyForm ) ):
        raise Exception( "Discrete two-body must be in the center-of-mass frame: LCT = %d MT = %d." % ( LCT, MT ) )

    firstDataLine = 2
    if( LVT != 0 ) :
        warningList.append( 'MF = 4, MT = 2 contains obsolete matrix used to transform Legendre coefficients between frames.' )
        firstDataLine += ( NK + 5 ) / 6

    info.logs.write( ' : MF=4, LTT = %s' % LTT )
    if( LTT == 0 ) :                # Purely isotropic angular distribution.
        subform = angularModule.isotropic( )
    elif( LTT == 1 ) :              # Legendre polynomial coefficient
        nextDataLine, angularData = endfFileToGNDMisc.getTAB2_Lists( firstDataLine, MF4Data, logFile = info.logs )
        subform = angularLegendreToPointwiseOrPiecewiseLegendre( MT, angularData, warningList, 4, 'LTT = 1' )
    elif( LTT == 2 ) :              # Tabulated probability distribution
        nextDataLine, angularTable = endfFileToGNDMisc.getTAB2_TAB1s( firstDataLine, MF4Data, logFile = info.logs )
        subform = convertAngularToPointwiseOrPiecewiseFromTAB2_TAB1( MT, angularTable, warningList )
    elif( LTT == 3 ) :              # Mixed Legendre and Tabulated probability distribution
        nextDataLine, angularData = endfFileToGNDMisc.getTAB2_Lists( firstDataLine, MF4Data, logFile = info.logs )
        nextDataLine, angularTable = endfFileToGNDMisc.getTAB2_TAB1s( nextDataLine, MF4Data, logFile = info.logs )
        subformPointwise = convertAngularToPointwiseOrPiecewiseFromTAB2_TAB1( MT, angularTable, warningList )
        subform = angularLegendreToPointwiseOrPiecewiseLegendre( MT, angularData, warningList, 4, 'LTT = 3', subformPointwise )
    else:
        raise Exception("Encountered unknown LTT=%d in MF4" % LTT)

    if( formClass is None ) : return( subform )
    form = formClass( info.style, frames[LCT], subform )
    product.distribution.add( form )
    return( form )

def readMF5( info, MT, MF5Data, warningList, delayNeutrons = False, product = None ) :

    ZA, AWR, dummy, dummy, NK, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF5Data[0], logFile = info.logs )
    ZA = int( ZA )
    printAWR_mode( info, MT, 5, ZA, AWR )
    countZAMasses( info, ZA, AWR )
    NK = int( NK )                 # number of partial energy distributions.
    dataLine = 1
    energySubforms, weights = [], []
    for k in xrange( NK ) :
        dataLine, productData = endfFileToGNDMisc.getTAB1( dataLine, MF5Data, logFile = info.logs )
        if( productData['NR'] > 1 ) :
            oldInterpolations = productData['interpolationInfo']
            if( oldInterpolations == [ [ 3, 1 ], [ 4, 2 ], [ 5, 1 ] ] ) :       # This is a kludge for about 5 data sets, but does a good job.
                productData['NR'] = 1
                productData['data'].insert( 1, [ ( 1. - FUDGE_EPS ) * productData['data'][1][0], productData['data'][0][1] ] )
                productData['interpolationInfo'] = [ [ len( productData['data'] ), 2 ] ]
            else :
                raise Exception( "Currently only one interpolation flag is supported" )
        LF = int( productData['L2'] )           # breakup flag
        xPrior, addWarning = None, True
        if( productData['data'][0][0] == productData['data'][1][0] ) : del productData['data'][0]
        weight = None
        if( ( NK > 1 ) and not( delayNeutrons ) ) :
            if( productData['NR'] != 1 ) : raise Exception( "Currently only one interpolation flag is supported for weight for MF=5, MT=%s" % MT )
            x1, y1 = productData['data'][0]
            discontinuities = []
            if( productData['interpolationInfo'][0][1] != 1 ) :
                for i1, xy in enumerate( productData['data'][1:] ) :                # Check if data can be treated as flat interpolation.
                    x2, y2 = xy
                    if( x1 != x2 ) :
                        if( y1 != y2 ) : raise Exception( 'Weight data must be convertible to flat interpolation' )
                    else :
                        discontinuities.insert( 0, i1 )
                    x1, y1 = x2, y2
                for discontinuity in discontinuities : del productData['data'][discontinuity]
            interpolation = endfFileToGNDMisc.ENDFInterpolationToGND1d( 1 )    # Flat interpolation.
            axes = axesModule.axes( labelsUnits = { 1 : ( 'energy_in' , 'eV' ), 0 : ( 'weight' , '' ) } )
            weight = XYsModule.XYs1d( data = productData['data'], axes = axes, accuracy = ENDF_Accuracy, interpolation = interpolation )
        else :
            for xy in productData['data'] :
                if( xPrior is not None ) :
                    if( xy[0] < xPrior ) : raise Exception( 'xy[0] = %s < xPrior = %s for MT=%d, MF=5' % ( xy[0], xPrior, MT ) )
                    if( xy[0] == xPrior ) :
                        xy[0] *= ( 1 + FUDGE_EPS )
                        if( addWarning ) : warningList.append( 'weights have same energies, second one being incremented for MT=%d, MF=5' % MT )
                        addWarning = False
                xPrior = xy[0]
        interpolation = endfFileToGNDMisc.ENDFInterpolationToGND1d( productData['interpolationInfo'][0][1] )
        axes = multiplicityModule.XYs1d.defaultAxes( multiplicityName = '' )
        weights.append( XYsModule.XYs1d( data = productData['data'], axes = axes, accuracy = ENDF_Accuracy, interpolation = interpolation ) )   # weights is only used for delayed nu_bar data

        info.logs.write( ' : MF=5, LF=%s' % LF )
        U = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( productData['C1'] ), 'eV' )    # Upper energy limit.
        if( LF == 1 ) :
            axes2d = energyModule.XYs2d.defaultAxes( )
            axes = axesModule.axes( )
            axes[0] = axes2d[0]
            axes[1] = axes2d[1]
            dataLine, EEpETable = endfFileToGNDMisc.getTAB2_TAB1s( dataLine, MF5Data, logFile = info.logs, axes = axes )
            subform = toPointwiseOrPiecewiseEnergy( MT, EEpETable )
        elif( LF == 5 ) :
            dataLine, TAB1, thetas = endfFileToGNDMisc.toEnergyFunctionalData( info, dataLine, MF5Data, 5, 'theta', 'eV' )
            dataLine, TAB1, gs = endfFileToGNDMisc.toEnergyFunctionalData( info, dataLine, MF5Data, 5, 'g', '',
                    xLabel = 'energy_out / theta( energy_in )', xUnit = '' )
            subform = energyModule.generalEvaporationSpectrum( U, thetas, gs )
        elif( LF == 7 ) :
            dataLine, TAB1, thetas = endfFileToGNDMisc.toEnergyFunctionalData( info, dataLine, MF5Data, 7, 'theta', 'eV' )
            subform = energyModule.simpleMaxwellianFissionSpectrum( U, thetas )
        elif( LF == 9 ) :
            dataLine, TAB1, thetas = endfFileToGNDMisc.toEnergyFunctionalData( info, dataLine, MF5Data, 9, 'theta', 'eV' )
            subform = energyModule.evaporationSpectrum( U, thetas )
        elif( LF == 11 ) :
            dataLine, TAB1, a = endfFileToGNDMisc.toEnergyFunctionalData( info, dataLine, MF5Data, 11, 'a', 'eV' )
            dataLine, TAB1, b = endfFileToGNDMisc.toEnergyFunctionalData( info, dataLine, MF5Data, 11, 'b', '1/eV' )
            subform = energyModule.WattSpectrum( U, a, b )
        elif( LF == 12 ) :
            dataLine, TAB1, Ts = endfFileToGNDMisc.toEnergyFunctionalData( info, dataLine, MF5Data, 12, 'T_M', 'eV' )
            EFL, EFH = [ PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( tmp ), 'eV' ) for tmp in ( TAB1['C1'], TAB1['C2'] ) ]
            subform = energyModule.MadlandNix( EFL, EFH, Ts )
        else :
            raise Exception( "Unsupported LF = %d" % LF )
        if( not( delayNeutrons ) and ( NK > 1 ) ) : subform.weight = weight

        energySubforms.append( subform )

    if( not( delayNeutrons ) ) :
        if( NK > 1 ) :
            info.logs.write( ' -- using energy.weightedFunctionals subform --' )
            weightedSubform = energyModule.weightedFunctionals( )
            for functional in energySubforms :
                weight = energyModule.weighted( functional.weight, functional )
                weightedSubform.addWeight( weight )
                del functional.weight
            energySubforms = weightedSubform
        else :
            energySubforms = energySubforms[0]
    return( energySubforms, weights )

def readMF6( MT, info, MF6Data, productList, warningList, undefinedLevelInfo, isTwoBody ) :

    ZA, AWR, dummy, LCT, NK, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF6Data[0], logFile = info.logs )
    ZA = int( ZA )
    printAWR_mode( info, MT, 6, ZA, AWR )
    countZAMasses( info, ZA, AWR )
    LCT = int( LCT )
    if( LCT == 0 ) :        # Happens for electro-atomic data.
        LCT = 1
        if( MT == 526 ) : LCT = 2       # Not sure this is right. Maybe it should be 1.
    LCTLight, LCTWeight = LCT, LCT
    if( LCT == 3 ) : LCTLight, LCTWeight = 2, 1
    NK = int( NK )                  # number of outgoing particle data sets

    dataLine, discreteGammas, discretePrimaryGammas = 1, {}, []
    info.logs.write( ' : MF=6' )
    for outGoing in xrange( NK ) :
        isLegendreConvertedToEnergy = False
        axes = multiplicityModule.XYs1d.defaultAxes( )
        dataLine, productData, multiplicityRegions = endfFileToGNDMisc.getTAB1Regions( dataLine, MF6Data, logFile = info.logs,
                axes = axes )
            # ZAP is ZA for outgoing particle; AWP is its atomic mass, LIP: 0 for residual in ground state, 1 for first excited state, etc
        ZAP, AWP, LIP, LAW, NP = int( productData['C1'] ), productData['C2'], productData['L1'], productData['L2'], productData['NR']
        if( ZAP < 2005 ) : LIP = 0 # LIP has multiple meanings. For light products, signifies different products with the same ZAP.
        if( ZAP == 11 ) :          # 11 is the ZAP for and electron, however, ENDF/B-VII.1 mislabels the photo as 11 also.
            ZAP = 0
            if( ( LAW == 8 ) or ( MT in [ 526, 528 ] ) or ( MT >= 534 ) ) : ZAP = 9
            if( ZAP == 0 ) : warningList.append( 'photon most likely labelled as an electron (ZAP = 11), converting to ZAP = 0' )
        printAWR_mode( info, MT, 6, ZAP, AWP )
        countZAMasses( info, ZAP, AWP )
        if( ZAP not in info.ZAMasses ) :
            info.ZAMasses[ZAP] = AWP * info.ZAMasses[1]
        elif( info.ZAMasses[ZAP] is None ) :
            info.ZAMasses[ZAP] = -AWP * info.ZAMasses[1]
        LCTp = LCTLight
        if( ZAP % 1000 > 4 ) : LCTp = LCTWeight
        LAW = int( LAW )
        frame = frames[LCTp]

        info.logs.write( ' : ZAP=%s, LAW=%s' % ( ZAP, LAW ) )
        isTwoBodyGamma = False
        form = None
        if( LAW == 0 ) :
            if( MT == 5 ) :                 # Need to check if this should be done for all reactions???????//
                form = unknownModule.form( info.style, frame )
        elif( LAW == 1 ) :              # Continuum Energy-Angle distributions
            dummy, dummy, LANG, LEP, NR, NE = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF6Data[ dataLine ], logFile = info.logs )
            LANG = int( LANG )          # identifies the type of data
            info.logs.write( ', LANG=%s' % LANG )
            LEP = int( LEP )            # interpolation type for outgoing energy
            interpolationLEP = endfFileToGNDMisc.ENDFInterpolationToGND1d( LEP )
            if( LANG == 1 ) :               # Legendre expansion
                NR = int( NR )              # number of interpolation regions for incident energy
                if( NR != 1 ) : raise Exception( "Currently only one interpolation region is supported for MF = 6, LAW = 1, LANG = 2; MT = %s" % MT )
                NE = int( NE )              # number of incident energies

                dataLine += 1
                EinInterpolationTypes = endfFileToGNDMisc.nStringsToInts( NR, dataLine, MF6Data, dimension = 2 )
                interpolationQualifier, EinInterpolation = endfFileToGNDMisc.ENDFInterpolationToGND2plusd( EinInterpolationTypes[0][1] )
                dataLine += 1 + ( NR - 1 ) / 3    # the next data is energy-angle distributions

                massRatio = AWR / ( 1. + AWR )
                maxLegendre = 0
                EEpClsData = []
                for EinCount in xrange( NE ) :
                    dummy, Ein, ND, NA, NW, NEP = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF6Data[dataLine], logFile = info.logs )
                    ND = int( ND )          # number of discrete gammas (nonzero only for gammas)
                    NA = int( NA )          # number of angular parameters (i.e., lMax).
                    NW = int( NW )          # number of data values for this incident energy
                    NEP = int( NEP )        # number of outgoing energy values
                    maxLegendre = max( maxLegendre, NA )
                    dataLine += 1

                    if( ND != 0 ) :
                        if( NA > 0 ) : raise Exception( "Logic here currently only supports isotropic scattering." )
                        discreteGammasAtE, EoutData = endfFileToGNDMisc.readDiscreteAndLegendre( ND, NEP - ND, dataLine, MF6Data,
                                dimension = 2 + NA, logFile = info.logs )
                        discretePrimaryGammasAtE = []
                        EgPrior = -1
                        for Eg, P in discreteGammasAtE :                    # Divide gammas into primary and discrete.
                            if( Eg == EgPrior ) : Eg -= float( "%0.e" % ( 1e-7 * Eg ) )
                            EgPrior = Eg
                            if( Eg < 0 ) :                                  # Primary gamma.
                                Epg = -Eg
                                if( Ein > 1e-2 ) : Epg -= massRatio * Ein   # Only include projectile correction for Ein > 1e-2 as typically Eg ~ 1e6.
                                discretePrimaryGammasAtE.append( [ Epg, [ Ein, P ] ] )
                            else :                                          # Discrete gamma.
                                if( Eg not in discreteGammas ) : discreteGammas[Eg] = []
                                discreteGammas[Eg].append( [ Ein, P ] )
                        if( len( discretePrimaryGammasAtE ) > 0 ) :
                            if( len( discretePrimaryGammas ) == 0 ) :
                                discretePrimaryGammas = discretePrimaryGammasAtE
                            else :      # Now we need to match level energies (aka binding energies) from prior with Ein's with current Ein.
                                        # This is needed since for different Ein's, the calculation of Epg will vary slightly (hopefully less than 1e-4).
                                if( len( discretePrimaryGammas ) != len( discretePrimaryGammasAtE ) ) :
                                    raise Exception( 'number of primary gammas at Ein = %s is different then first incident energy' % Ein )
                                for index, Eg in enumerate( discretePrimaryGammas ) :
                                    if( abs( Eg[0] - discretePrimaryGammasAtE[index][0] ) > 1e-4 * Eg[0] ) : raise \
                                        Exception( 'primary energy of %s is not close to primary energy %s' % ( Eg[0], discretePrimaryGammasAtE[index][0] ) )
                                    Eg.append( discretePrimaryGammasAtE[index][1] )
                    else :
                        EoutData = endfFileToGNDMisc.nFunkyFloatStringsToFloats( NEP, dataLine, MF6Data, dimension = 2 + NA, logFile = info.logs )
                    if( len( EoutData ) > 0 ) :
                        # Test for some common problems.
                        # 1: all outgoing energies == 0 for continuum gammas.
                        # This usually only affects one incident energy, so just replace that energy with empty outgoing distribution:
                        energy_out_list = [ e[0] for e in EoutData ]
                        if( sum( energy_out_list ) == 0 ) :
                            warningList.append("At Ein=%s eV, continuum gamma outgoing energies are all 0.0 (MT=%d, ZAP=%d)!"
                                    % ( Ein, MT, ZAP ) )
                            EoutData = [ [ 0.0, 0.0 ], [ 1.0, 0.0 ] ]
                        # 2: trouble with duplicate outgoing energies:
                        elif( ( max( [ energy_out_list.count( a ) for a in energy_out_list ] ) > 2 ) or ( energy_out_list.count( 0.0 ) > 1 ) ) :
                            warningList.append( "Too many duplicate outgoing energies for Ein=%s eV (MT=%d, ZAP=%d)!"
                                    % ( Ein, MT, ZAP ) )
                            tmp = []
                            i1 = 0
                            if( energy_out_list.count( 0.0 ) > 1 ) :
                                finalZeroIndex = energy_out_list.count( 0.0 ) - 1
                                energy_out_list = energy_out_list[finalZeroIndex:]
                                EoutData = EoutData[finalZeroIndex:]
                            while( i1 < len( energy_out_list ) ) :
                                eout = energy_out_list[i1]
                                eoutCount = energy_out_list.count( eout )
                                if( eoutCount > 2 ) :
                                    tmp.extend( [ EoutData[i1], EoutData[i1+eoutCount-1] ] )
                                    i1 += eoutCount
                                else :
                                    tmp.append( EoutData[i1] )
                                    i1 += 1
                            EoutData = tmp                          # End of TENDL-specific tests.
                        EpClsDatas, EpPrior, ClsPrior = [], -1, []
                        if( LEP == 2 ) :                            # Remove some data that are not needed and complicates things.
                            if( len( EoutData ) > 2 ) :
                                if( EoutData[-1][0] == EoutData[-2][0] ) :
                                     if( EoutData[-1][1] == EoutData[-2][1] ) : EoutData.pop()
                                     elif (EoutData[-1][1] == 0):   # Discontinuity used to make final point == 0. Bump previous energy back slightly
                                         EoutData[-2][0] = nudgeValue( EoutData[-2][0], -1 )
                        regions = []
                        for i1, EpCs in enumerate( EoutData ) :
                            e_out, Cls = EpCs[0], EpCs[1:]
                            if( e_out == EpPrior ) :
                                if( Cls == ClsPrior ) :
                                    if( skippedLast ) : raise 'hell - need to fix'
                                    warningList.append( 'skipping duplicate e_out = %s, i1 = %s %s %s' % ( e_out, i1, EinCount, Ein ) )
                                    skippedLast = True
                                    continue
                                else:   # create new region:
                                    regions.append( EpClsDatas )
                                    EpClsDatas = []
                            skippedLast = False
                            EpClsData = [ e_out, Cls ]
                            EpClsDatas.append( EpClsData )
                            EpPrior, ClsPrior = e_out, Cls
                        if regions:
                            regions.append( EpClsDatas )
                            EEpClsData.append( [ Ein, regions, True ] )
                        else:
                            EEpClsData.append( [ Ein, EpClsDatas, False ] )
                    dataLine += 1 + ( NW -  1 ) / 6    # Offset for the next incident energy

                if( maxLegendre == 0 ) :    # Only have l=0 for each outgoing energy. Convert this to
                                            # uncorrelated with P(E_out|E_in) and isotropic angular distribution.
                    if( len( EEpClsData ) != 0 ) :  # length == 0 happens when there are only discrete gammas and no continuum gamma data.
                        isLegendreConvertedToEnergy = True
                        angularSubform = angularModule.isotropic( )

                        axes = energyModule.XYs2d.defaultAxes( energyUnit = "eV" )
                        EPrior = -1
                        energySubform = energyModule.regions2d( axes = axes )
                        energySubformRegion = energyModule.XYs2d( axes = axes, interpolation = EinInterpolation,
                                interpolationQualifier = interpolationQualifier )
                        axes2 = axesModule.axes( rank = 2 )
                        axes2[0] = axes[0]
                        axes2[1] = axes[1]
                        for index, ( Ein, EpCls, multiRegion ) in enumerate( EEpClsData ) :
                            if multiRegion:
                                xData = energyModule.regions1d( value = Ein, axes = axes2 )
                                for region in EpCls:
                                    EpP = [ [ Ep, P[0] ] for Ep, P in region ]
                                    xData.append( energyModule.XYs1d( data = EpP, interpolation = interpolationLEP ) )
                            else:
                                EpP = [ [ Ep, P[0] ] for Ep, P in EpCls ]
                                xData = energyModule.XYs1d( data = EpP, interpolation = interpolationLEP, 
                                        value = Ein, axes = axes2 )
                            if( Ein == EPrior ) :
                                energySubform.append( energySubformRegion )
                                energySubformRegion = energyModule.XYs2d( axes = axes, interpolation = EinInterpolation,
                                interpolationQualifier = interpolationQualifier )
                            energySubformRegion.append( xData )
                            EPrior = Ein
                        if( len( energySubform ) > 0 ) :
                            energySubform.append( energySubformRegion )
                        else :
                            energySubform = energySubformRegion

                        form = uncorrelated( info.style, frame, angularSubform, energySubform )
                else :
                    axes = axesModule.axes( rank = 4 )
                    axes[3] = axesModule.axis( 'energy_in', 3, 'eV' )
                    axes[2] = axesModule.axis( 'energy_out', 2, 'eV' )
                    axes[1] = axesModule.axis( 'l', 1, '' )
                    axes[0] = axesModule.axis( 'C_l(energy_out|energy_in)', 0, '1/eV' )
                    EPrior = -1
                    energyAngularSubform = energyAngularModule.regions3d( axes = axes )
                    energyAngularSubformRegion = energyAngularModule.XYs3d( axes = axes, interpolation = EinInterpolation,
                            interpolationQualifier = interpolationQualifier )
                    for i1, ( e_in, EpCls, multiRegion ) in enumerate( EEpClsData ) :
                        if( multiRegion ) : raise NotImplemented
                        multiD_2d = energyAngularModule.XYs2d( value = e_in, interpolation = interpolationLEP )
                        for i2, ( e_out, Cls ) in enumerate( EpCls ) :
                            multiD_2d.append( energyAngularModule.Legendre( coefficients = Cls, value = e_out ) )
                        if( e_in == EPrior ) :
                            energyAngularSubform.append( energyAngularSubformRegion )
                            energyAngularSubformRegion = energyAngularModule.XYs3d( axes = axes, interpolation = EinInterpolation,
                            interpolationQualifier = interpolationQualifier )
                        energyAngularSubformRegion.append( multiD_2d )
                        EPrior = e_in
                    if( len( energyAngularSubform ) > 0 ) :
                        energyAngularSubform.append( energyAngularSubformRegion )
                    else :
                        energyAngularSubform = energyAngularSubformRegion
                    form = energyAngularModule.form( info.style, frame, energyAngularSubformRegion )

            elif( LANG == 2 ) :             # Kalbach-Mann data
                if( LCTp != 2 ) : raise Exception( 'LCT = %s != 2 as required for Kalbach-Mann data for MF=6, MT=%s' % ( LCTp, MT ) )
                dataLine, KalbachMannData = endfFileToGNDMisc.getTAB2_Lists( dataLine, MF6Data, logFile = info.logs )
                if( KalbachMannData['NR'] != 1 ) :
                    raise Exception( "Currently only one interpolation flag is supported for MF = 6, LAW = 1, LANG = 2; MT = %s" % MT )
                interpolationQualifier, interpolation = endfFileToGNDMisc.ENDFInterpolationToGND2plusd( KalbachMannData['interpolationInfo'][0][1] )
                NA = int( KalbachMannData['Lists'][0]['L2'] )
                dataPerPoint = NA + 2

                if( interpolationQualifier == standardsModule.interpolation.noneQualifierToken ) :
                    interpolationQualifier = standardsModule.interpolation.unitBaseToken
                axes = KalbachMannModule.subform.defaultAxes( 'f' )
                fData = multiD_XYsModule.XYs2d( axes = axes, interpolation = interpolation, 
                        interpolationQualifier = interpolationQualifier )

                if( interpolationQualifier == standardsModule.interpolation.unitBaseToken ) :
                    interpolationQualifier = standardsModule.interpolation.unitBaseUnscaledToken
                else :
                    interpolationQualifier = standardsModule.interpolation.correspondingPointsUnscaledToken
                axes = KalbachMannModule.subform.defaultAxes( 'r' )
                rData = multiD_XYsModule.XYs2d( axes = axes, interpolation = interpolation, 
                        interpolationQualifier = interpolationQualifier )

                aData = None
                if( NA == 2 ) :
                    axes = KalbachMannModule.subform.defaultAxes( 'a' )
                    aData = multiD_XYsModule.XYs2d( axes = axes, label = "a", interpolation = interpolation, 
                            interpolationQualifier = interpolationQualifier )

                priorE, priorEp_f_r_ = -1, [ -1, -1, -1 ]
                for i1, data in enumerate( KalbachMannData['Lists'] ) :
                    value = data['C2']
                    Ep_ = data['data'][::dataPerPoint]
                    f_ = data['data'][1::dataPerPoint]
                    r_ = data['data'][2::dataPerPoint]
                    priorEp_f_r__ = [ Ep_, f_, r_ ]
                    if( value == priorE ) :
                        if( i1 == 1 ) :
                            fData.pop( 0 )
                            rData.pop( 0 )
                            if( aData is not None ) : aData.pop( 0 )
                        else :
                            if( priorEp_f_r_ == priorEp_f_r__ ) : continue        # For TENDL files with duplicate data.
                            print '\nMT=%d' % MT
                            print value, priorEp_f_r_
                            print priorE, priorEp_f_r__
                            raise 'hell - need to support regions'
                    priorE, priorEp_f_r_ = value, priorEp_f_r__
                    fData.append( XYsModule.XYs1d( data = ( Ep_, f_ ), dataForm = 'xsandys',
                                                 value = value, interpolation = interpolationLEP ) )
                    rData.append( XYsModule.XYs1d( data = ( Ep_, r_ ), dataForm = 'xsandys',
                                                 value = value, interpolation = interpolationLEP ) )
                    if( aData is not None ) :
                        a_ = data['data'][3::dataPerPoint]
                        aData.append( XYsModule.XYs1d( data = ( Ep_, a_ ), dataForm = 'xsandys',
                                                     value = value, interpolation = interpolationLEP ) )

                fSubform = KalbachMannModule.fSubform( fData )
                rSubform = KalbachMannModule.rSubform( rData )
                aSubform = KalbachMannModule.aSubform( aData )
                form = KalbachMannModule.form( info.style, frame, fSubform, rSubform, aSubform )
            else :
                raise Exception( "Unsupported LANG = %d for continuum energy-angle distribution MF = 6: ZAP = %d, LAW = %d: MT = %d" % \
                        ( LANG, ZAP, LAW, MT ) )
        elif( LAW == 2 ) :
            if( LCT != 2 ): raise Exception( "Discrete two-body must be in the center-of-mass frame: LCT = %d MT = %d." % ( LCT, MT ) )
            isTwoBody = True
            if( ZAP == 0 ) : isTwoBodyGamma = True
            dataLine, angularData = endfFileToGNDMisc.getTAB2_Lists( dataLine, MF6Data, logFile = info.logs )
            LANG = int( angularData['Lists'][0]['L1'] )
            info.logs.write( ', LANG=%s' % LANG )
            if( angularData['NR'] != 1 ) :
                raise Exception( "Currently only one interpolation flag is supported for MF = 6, LAW = 2, LANG = %s; MT = %s" % ( LANG, MT ) )
            interpolationQualifier, interpolationE_in = endfFileToGNDMisc.ENDFInterpolationToGND2plusd( angularData['interpolationInfo'][0][1] )

            if( LANG == 0 ) :
                angularSubform = angularLegendreToPointwiseOrPiecewiseLegendre( MT, angularData, warningList, 6, 'LAW = 2, LANG = 0' )
                form = angularModule.twoBodyForm( info.style, frame, angularSubform )
            elif( LANG in [ 12, 14 ] ) :
                angularSubform = convertAngularToPointwiseOrPiecewiseFromTAB2_List( MT, LANG, angularData, warningList )
                form = angularModule.twoBodyForm( info.style, frame, angularSubform )
            else :
                raise Exception( "LANG = %d for LAW = %d not supported: MT = %d" % ( LANG, LAW, MT ) )
            if( ( ZAP == 0 ) and ( AWP != 0 ) ) :
                energySubform = energyModule.constant( PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( AWP ), 'eV' ) )
                form = uncorrelated( info.style, frame, angularSubform, energySubform )
        elif( LAW == 3 ) :
            subform = angularModule.isotropic( )
            form = angularModule.twoBodyForm( info.style, frame, subform )
        elif( LAW == 4 ) :
            subform = angularModule.recoil( product.distribution[info.style] )
            form = angularModule.twoBodyForm( info.style, frame, subform )
        elif( LAW == 5 ) :  # charged-particle elastic scattering
            assert LCT == 2, "Charged-particle elastic must be in the center-of-mass frame: LCT = %d MT = %d." % ( LCT, MT )
            dataLine, angularData = endfFileToGNDMisc.getTAB2_Lists( dataLine, MF6Data, logFile = info.logs )
            SPI = angularData['C1']
            LIDP = angularData['L1']
            info.particleList[ info.projectile ].attributes['spin'] = xParticleModule.spin(SPI)
            LTP = int( angularData['Lists'][0]['L1'] )
            info.logs.write( ', LTP=%s' % LTP )
            interpolationQualifier, interpolationE_in = endfFileToGNDMisc.ENDFInterpolationToGND2plusd( angularData['interpolationInfo'][0][1] )
            # LTP flag changes interpretation of the data:
            if( LTP == 1 ) :
                (nuclear, real, imaginary) = convertNuclearPlusInterferenceDataToPiecewise( MT, angularData, warningList, 6, 'LAW = 5, LTP = %d'%LTP, LIDP )
                form = angularModule.CoulombExpansionForm( info.style, frame, nuclear, real, imaginary )
            elif( LTP == 2 ) :
                raise NotImplemented( "MF=6 LAW=5 LTP=2 not yet implemented (MT%d)!" % MT )
            elif( LTP in ( 12, 14, 15 ) ) :
                subform = convertAngularToPointwiseOrPiecewiseFromTAB2_List( MT, LTP, angularData, warningList )
                form = angularModule.NuclearPlusCoulombInterferenceForm( info.style, frame, subform )
            else:
                raise Exception( "unknown LTP encountered for MF=6, LAW=5, MT=%s" % MT )
        elif( LAW == 6 ) :
            APSX, dummy, dummy, dummy, dummy, NPSX = endfFileToGNDMisc.sixFunkyFloatStringsToFloats( MF6Data[ dataLine ], logFile = info.logs )
            dataLine += 1
            APSX *= info.masses.getMassFromZA( 1 )
            angularSubform = angularModule.isotropic( )
            energySubform = energyModule.NBodyPhaseSpace( int( NPSX ), PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( APSX ), 'amu' ) )
            # Some data has the wrong frame, should always be com hence frames[2].
            form = uncorrelated( info.style, frames[2], angularSubform, energySubform )

        elif( LAW == 7 ) :
            dataLine, NEData = endfFileToGNDMisc.getTAB2Header( dataLine, MF6Data, logFile = info.logs )
            NR = int( NEData['NR'] )                    # number of interpolation regions for this incident energy
            if( NR != 1 ) : raise Exception( "Currently only one interpolation flag is supported for MF = 6, LAW = 7; MT = %s" % MT )
            axes = angularEnergyModule.XYs3d.defaultAxes( )
            interpolationQualifier, interpolation = endfFileToGNDMisc.ENDFInterpolationToGND2plusd( NEData['interpolationInfo'][0][1] )
            angularEnergySubform = angularEnergyModule.XYs3d( axes = axes, interpolation = interpolation,
                    interpolationQualifier = interpolationQualifier )
            axes2d = axesModule.axes( 3 )
            axes2d[0] = angularEnergySubform.axes[0]
            axes2d[1] = angularEnergySubform.axes[1]
            axes2d[2] = angularEnergySubform.axes[2]
            axes1d = axesModule.axes( 2 )
            axes1d[0] = angularEnergySubform.axes[0]
            axes1d[1] = angularEnergySubform.axes[1]
            for i1 in xrange( int( NEData['NZ'] ) ) :   # Loop over incident energies
                dataLine, muEpPTable = endfFileToGNDMisc.getTAB2_TAB1s( dataLine, MF6Data, logFile = info.logs, axes = axes1d )
                muEpP = muEpPTable['TAB2s']
                if( len( muEpP ) != 1 ) : raise 'hell - need to fix'
                interpolation, muEpP = muEpP[0]
                interpolationQualifier, interpolation = endfFileToGNDMisc.ENDFInterpolationToGND2plusd( muEpPTable['interpolationInfo'][0][1] )
                muEpP_multiD = angularEnergyModule.XYs2d( value = muEpPTable['C2'], 
                        interpolation = interpolation, interpolationQualifier = interpolationQualifier, axes = axes2d )
                for i2, EpP in enumerate( muEpP ) :
                    if( len( EpP ) > 1 ) : raise 'hell - need to fix'
                    muEpP_multiD.append( angularEnergyModule.XYs1d.returnAsClass( EpP[0], EpP[0] ) )
                angularEnergySubform.append( muEpP_multiD )
            form = angularEnergyModule.form( info.style, frame, angularEnergySubform )
        elif( LAW == 8 ) :
            angularSubform = angularModule.forward( )

            dataLine, TAB1, energyLoss = endfFileToGNDMisc.getTAB1Regions( dataLine, MF6Data, logFile = info.logs )
            if( len( energyLoss ) != 1 ) : raise 'hell - fix me'
            axes = energyModule.energyLoss.defaultAxes( )
            energySubform = energyModule.energyLoss( data = energyLoss[0], axes = axes, interpolation = energyLoss[0].interpolation )

            form = uncorrelated( info.style, frame, angularSubform, energySubform )
        else :
            raise Exception( "LAW = %d not implemented: MT = %d." %  ( LAW, MT ) )

        multiplicityData = productData['data']
        rangeMin, integerMultiplicity = multiplicityData[0][1], False
        if( ( ZAP != 0 ) and ( int( rangeMin ) == rangeMin ) ) :
            integerMultiplicity = True
            for x, y in multiplicityData : integerMultiplicity &= ( rangeMin == y )
        if( integerMultiplicity ) :
                multiplicity = int( rangeMin )
        else :
            multiplicity = getMultiplicityPointwiseOrPieceWise( info, multiplicityRegions, warningList )
        if( isinstance( multiplicity, multiplicityModule.XYs1d ) and ( multiplicity.rangeMin( ) < 0 ) ) :
            warningList.append( "Negative multiplicity encountered for MF6, MT%d %s" %
                ( MT, toGNDMisc.getTypeNameENDF( info, ZAP, undefinedLevelInfo ) ) )

        if( ( ZAP == 0 ) and ( AWP == 0 ) ) : # Gammas. Appear to always be stored using LAW = 1 and LANG = 1.

            def addGammaProduct( form, multiplicity ) :

                if( isinstance( multiplicity, multiplicityModule.regions1d ) ) :
                    if( len( multiplicity ) == 1 ) :
                        label = multiplicity.label
                        multiplicity = multiplicity[0]
                        multiplicity.label = label
                if( isinstance( multiplicity, multiplicityModule.XYs1d ) ) :
                    if( len( multiplicity ) < 2 ) : return( None )
                product = toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, ZAP, None ), multiplicity = multiplicity )
                product.distribution.add( form )
                productList.append( product )
                return( product )       # May be required for LAW = 4 logic.

            def addGammaAdjustWeight( angularSubform, energySubform, totalGammaMultiplicity ) :
                """
                Have Legendre section with total multiplicity and only l=0 for discrete gammas.
                Convert to isotropic angular distributions, plus adjust weights.
                """

                def getPointwiseMultiplicity( angularSubform, totalGammaMultiplicity, EPrior ) :

                    data = []
                    for energy_in in angularSubform :
                        multiplicity = getMultiplicity( totalGammaMultiplicity, EPrior, energy_in.value )
                        data.append( [ energy_in.value, energy_in.coefficients[0] * multiplicity ] )
                        EPrior = energy_in.value
                    multiplicity = multiplicityModule.XYs1d( data = data, axes = axes )
                    return( EPrior, multiplicity )

                if( len( angularSubform ) < 2 ) : return( None )
                axes = multiplicityModule.XYs1d.defaultAxes( )
                if( isinstance( angularSubform, angularModule.XYs2d ) ) :
                    EPrior, multiplicity = getPointwiseMultiplicity( angularSubform, totalGammaMultiplicity, angularSubform[0].value )
                    multiplicity.label = info.style
                elif( isinstance( angularSubform, angularModule.regions2d ) ) :
                    multiplicity = multiplicityModule.regions1d( label = info.style, axes = axes )
                    EPrior = -1
                    for region in angularSubform :
                        EPrior, multiplicityRegion = getPointwiseMultiplicity( region, totalGammaMultiplicity, EPrior )
                        multiplicity.append( multiplicityRegion )
                else :
                    raise Exception( 'Unsupport angular subform = "%s"' % angularSubform.moniker )

                product = toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, ZAP, None ), multiplicity = multiplicity )

                angularSubform = angularModule.isotropic( )
                form = uncorrelated( info.style, frame, angularSubform, energySubform )
                product.distribution.add( form )
                productList.append( product )
                return( product )       # May be required for LAW = 4 logic.

            def addPrimaryOrDiscreteGamma( distinctType, Eg, ELegenres, totalGammaMultiplicity, frame ) :

                if( distinctType == 'discrete' ) :
                    energySubform = energyModule.constant( PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( Eg ), 'eV' ) )
                else :
                    energySubform = energyModule.primaryGamma( PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( Eg ), 'eV' ) )

                axes = angularModule.XYs2d.defaultAxes( asLegendre = True )
                angularSubform = angularModule.XYs2d( axes = axes )
                EPrior = -1
                maxLegendre = 0
                angularRegion = angularModule.regions2d( axes = axes )
                for i1, ( energy, coefficients ) in enumerate( ELegenres ) :
                    maxLegendre = max( maxLegendre, len( coefficients ) )
                    if( energy == EPrior ) :
                        angularRegion.append( angularSubform )
                        angularSubform = angularModule.XYs2d( axes = axes )
                    angularSubform.append( angularModule.Legendre( coefficients = coefficients, value = energy ) )
                    EPrior = energy
                if( len( angularRegion ) > 0 ) :
                    angularRegion.append( angularSubform )
                    angularSubform = angularRegion

                if( maxLegendre == 1 ) :    # Only have l=0 for each incident energy. Distribution needs to be normalized.
                    return( addGammaAdjustWeight( angularSubform, energySubform, totalGammaMultiplicity ) )

                if( isTwoBodyGamma ) :
                    raise 'hell - when does this happen'
                    form = angularModule.twoBodyForm( info.style, frame, angularSubform )
                else :
                    form = uncorrelated( info.style, frame, angularSubform, energySubform )

                return( addGammaProduct( form, totalGammaMultiplicity ) )

            totalGammaMultiplicity = multiplicity               # Makes the name more explicit.
            if( ( len( discreteGammas ) + len( discretePrimaryGammas ) ) > 0 ) :
                info.totalMF6_12_13Gammas[MT] = [ 6, multiplicity ]

            if( form is not None ) :                       # Continuum gamma present

                if( discreteGammas or discretePrimaryGammas ) :

                    def fixContinuumGammaSpectrum( energyDist ) :

                        EpPrior = 2e-5
                        energyDist[0] = ( 0, 1e-16 )            # Make distribution with small probability near 0.
                        energyDist[1] = ( EpPrior, 0 )
                        for i1 in range( 2, len( energyDist ) ) :
                            Ep, P1 = energyDist[i1]
                            if( Ep < EpPrior ) : Ep = 2 * EpPrior
                            EpPrior = Ep
                            energyDist[i1] = ( Ep, 0 )

                    if(   isinstance( totalGammaMultiplicity, XYsModule.XYs1d ) ) :
                        newMultiplicity = multiplicityModule.XYs1d( axes = totalGammaMultiplicity.axes,
                                interpolation = totalGammaMultiplicity.interpolation, label = info.style )
                    elif( isinstance( totalGammaMultiplicity, regionsModule.regions ) ) :
                        newMultiplicity = multiplicityModule.regions1d( axes = totalGammaMultiplicity.axes, label = info.style )
                    else :
                        raise Exception( 'Unsupported multiplicity form "%s"' % totalGammaMultiplicity.moniker )
                    if( isinstance( energySubform, energyModule.XYs2d ) ) : energySubform = [ energySubform ]
                    EPrior = energySubform[0][0].value
                    for i1, region in enumerate( energySubform ) :
                        for i2, energyDist in enumerate( region ) :
                            Ein = energyDist.value
                            integral = energyDist.integrate( )
                            if( integral == 0 ) :
                                if( not( info.continuumSpectraFix ) ) :
                                    raise Exception( "Zero norm continuum gamma spectrum energies = %s (MT=%d). Try option '--continuumSpectraFix'" %
                                            ( energyDist.value, MT ) )
                                fixContinuumGammaSpectrum( energyDist )
                                integral = energyDist.integrate( )
                            if(   isinstance( newMultiplicity, XYsModule.XYs1d ) ) :
                                newMultiplicity.setValue( Ein, getMultiplicity( totalGammaMultiplicity, EPrior, Ein ) * integral )
                            elif( isinstance( newMultiplicity, regionsModule.regions ) ) :
                                if( Ein == EPrior ) :
                                    if( ( i1 + i2 ) > 0 ) : newMultiplicity.append( newRegionMultiplicity )
                                    newRegionMultiplicity = multiplicityModule.XYs1d( axes = totalGammaMultiplicity.axes )
                                newRegionMultiplicity.setValue( Ein, getMultiplicity( totalGammaMultiplicity, EPrior, Ein ) * integral )
                            region[i2] = region[i2].normalize( )
                            EPrior = Ein
                    if( isinstance( energySubform, list ) ) : energySubform = energySubform[0]
                    form = uncorrelated( info.style, form.productFrame, form.angularSubform.data, energySubform )
                    if( isinstance( newMultiplicity, multiplicityModule.regions1d ) ) :
                        newMultiplicity.append( newRegionMultiplicity )
                    product = addGammaProduct( form, newMultiplicity )
                else:
                    product = addGammaProduct( form, totalGammaMultiplicity )
            for Eg in sorted( discreteGammas ) :
                product = addPrimaryOrDiscreteGamma( 'discrete', Eg, discreteGammas[Eg], totalGammaMultiplicity, frame )
            for EgEinPs in discretePrimaryGammas :
                product = addPrimaryOrDiscreteGamma( 'primary', EgEinPs[0], EgEinPs[1:], totalGammaMultiplicity, frame )
        else :                          # Non gamma particle.
            if( info.targetLevel > 0 ) and ( MT in range( 51, 90 ) ) and ( ZAP == info.target.getZ_A_SuffixAndZA()[-1] ) :
                # for isomeric targets, need to adjust the product level: MT51 goes to ground state, etc.
                if( ( MT - 50 ) <= info.targetLevel ) : undefinedLevelInfo['levelIndex'] -= 1
            thisParticle = toGNDMisc.getTypeNameENDF( info, ZAP, undefinedLevelInfo )
            if( LIP > 0 ) :
                thisParticleBaseName = thisParticle.name
                metaStableName = aliasModule.aliases.nuclearMetaStableName( thisParticleBaseName, LIP )
                if( info.reactionSuite.hasAlias( metaStableName ) ) :
                    levelName = info.reactionSuite.getAlias( metaStableName ).getValue( )
                    thisParticle = info.reactionSuite.getParticle( levelName )
                else :
                    try :
                        levelIndex, level = metaStableData[metaStableName]
                        thisParticle = toGNDMisc.getTypeNameGamma( info, ZAP, level = level, levelIndex = levelIndex )
                        info.reactionSuite.addNuclearMetaStableAlias( thisParticleBaseName, thisParticle.name, LIP )
                    except :
                        raise KeyError( 'Meta stable data not available for %s' % metaStableName )
            product = toGNDMisc.newGNDParticle( info, thisParticle, multiplicity = multiplicity )
            if( form is not None ) : product.distribution.add( form )

            productList.append( product )
            if( isLegendreConvertedToEnergy ) : product.addAttribute( 'ENDFconversionFlag', 'MF6' )

        if( ( len( discreteGammas ) > 0 ) and ( productData['interpolationInfo'][0][1] != 2 ) ) :
            info.logs.write( 'interpolation %s is not linear for gamma multiplicity' % productData['interpolationInfo'][0][1] )

    for product in productList :
        if 'ENDFconversionFlag' not in product.attributes:
            product.addAttribute( 'ENDFconversionFlag', 'MF6' )
    return( isTwoBody )

def readMF8_454_459( info, MT, MTData, warningList ) :
    "Parse fission yield data."

    if( MT not in [ 454, 459 ] ) : raise Exception( 'MT most be 454 or 459, not %s' % MT )
    MF8Data = MTData[8]
    ZA, AWR, LE1 = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF8Data[0], intIndices = [ 0, 2 ], logFile = info.logs )[0:3]
    iLine = 1
    allData = {}
    for iLE in range( LE1 ):
        iLine, thisData = endfFileToGNDMisc.getList( iLine, MF8Data, logFile = info.logs )
        NFP = thisData[ 'N2' ]
        E = thisData[ 'C1' ]
        for iFP in range( NFP ):
            ZAFP = int( thisData[ 'data' ][ iFP*4 ] )
            FPS =  int( thisData[ 'data' ][ iFP*4 + 1 ] )
            Y =    thisData[ 'data' ][ iFP*4 + 2 ]
            DY =   thisData[ 'data' ][ iFP*4 + 3 ]
            key = ( ZAFP, FPS )
            if( key not in allData ) : allData[key] = []
            allData[key].append( ( E, Y, DY ) )
    return( allData )

def readMF8( info, MT, MTData, warningList ) :
    "Regular decay data."

    firstLMF, radioactiveDatas = None, []
    if( 8 in MTData ) :
        dataLine, MF8Data = 1, MTData[8]
        ZA, AWR, LIS, LISO, NS, NO = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF8Data[0], intIndices = [ 0, 2, 3, 4, 5 ], logFile = info.logs )
        countZAMasses( info, ZA, AWR )
        printAWR_mode( info, MT, 8, ZA, AWR )
        MF9Data = readMF9or10( info, MT, MTData, 9, LIS, warningList )
        MF10Data = readMF9or10( info, MT, MTData, 10, LIS, warningList )
        metastables = {}
        for idx in xrange( NS ) :
            ZAP, ELFS, LMF, LFS, ND6, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF8Data[dataLine], intIndices = [ 0, 2, 3, 4 ], logFile = info.logs )

            if( ZAP not in metastables ) : metastables[ZAP] = 0
            if( LFS == 0 and ELFS != 0 ) :
                warningList.append( "MF8 claims non-zero ELFS = %s for the ground state, MT = %s. Converting ELFS to 0" % (ELFS, MT ) )
                info.doRaise.append( warningList[-1] )
                ELFS = 0
            if( ELFS != 0 ) : metastables[ZAP] += 1

            if( firstLMF is None ) : firstLMF = LMF
            if( LMF != firstLMF ) : raise Exception( 'LMF changing from %s to %s is not allowed' % ( firstLMF, LMF ) )
            dataLine += 1

            crossSection, multiplicity = None, 1
            QM, QI = None, None
            if( LMF == 3 ) :
                pass
            elif( LMF == 6 ) :  # multiplicity in MF6, which hasn't been read. Set multiplicity to 'unknown', to be overridden later.
                multiplicity = multiplicityModule.unknown( 'unknown' )
            elif( LMF in [ 9, 10 ] ) :
                if( ( ( LMF == 9 ) and not( MF9Data ) ) or ( ( LMF == 10 ) and not( MF10Data ) ) ) :    # BRB ????? Why are we checking for bad data.
                    LMF1 = ( 10 if( LMF == 9 ) else 9 )
                    warningList.append( 'LMF = %d, but no MF%d found for MT%d. Trying LMF=%d instead' % ( LMF, LMF, MT, LMF1 ) )
                    info.doRaise.append( warningList[-1] )
                    if( ( ( LMF1 == 9 ) and not( MF9Data ) ) or ( ( LMF1 == 10 ) and not( MF10Data ) ) ) :
                        # neither MF9 or MF10 exist.
                        continue
                    LMF = LMF1
                if( LMF == 9 ) :
                    TAB1, multiplicity = MF9Data[idx]
                    multiplicity = getMultiplicityPointwiseOrPieceWise( info, multiplicity, warningList )
                else :
                    TAB1, crossSection = MF10Data[idx]
                QM, QI, IZAP, LFS9or10 = TAB1['C1'], TAB1['C2'], int( TAB1['L1'] ), int( TAB1['L2'] )
                ELFS9or10 = QM - QI
                if( abs( ELFS - ELFS9or10 ) > 2e-4 * abs( ELFS ) ) :
                    warningList.append( "MF8 residual level energy = %s for level %s of ZA = %d not close to MF%s's value = %s for MT = %s"
                            % ( ELFS, LIS, ZAP, LMF, ELFS9or10, MT ) )
                    info.doRaise.append( warningList[-1] )
                if( LFS != LFS9or10 ):
                    warningList.append("MF8 claims level index = %d but MF9/10 claim level index = %d" % (LFS, LFS9or10))
                    info.doRaise.append( warningList[-1] )

            try :
                residual = toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameGamma( info, ZAP, level = ELFS, levelIndex = LFS ), multiplicity = multiplicity )
            except :
                info.logs.write( '\nMT = %s\n' % MT )
                raise
            radioactiveDatas.append( [ ZAP, residual, crossSection, LFS, QI ] )

            if( metastables[ZAP] and ( ZAP != 0 ) ) :
                aliasName = aliasModule.aliases.nuclearMetaStableName( residual.particle.groundState.name, metastables[ZAP] )
                if( not( info.reactionSuite.hasAlias( aliasName ) ) ) :
                    info.reactionSuite.addNuclearMetaStableAlias( residual.particle.groundState.name, residual.particle.name, metastables[ZAP] )

            if( NO == 0 ) : dataLine += ( ND6 + 5 ) / 6
    return( firstLMF, radioactiveDatas )

def readMF9or10( info, MT, MTData, MF, targetLIS, warningList ) :

    if( MF not in MTData.keys( ) ) : return( None )
    dataLine, MFData, MF9or10 = 1, MTData[MF], []
    ZA, AWR, LIS, dummy, NS, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MFData[0], intIndices = [ 0, 2, 4 ], logFile = info.logs )
    printAWR_mode( info, MT, MF, ZA, AWR )
    countZAMasses( info, ZA, AWR )
    if( LIS != targetLIS ) :
        warningList.append( "residual's LIS = %s does not match target's LIS = %s: MT=%d, MF=%d" % ( LIS, targetLIS, MT, MF ) )
        info.doRaise.append( warningList[-1] )
        return
    if( MF == 10 ) :
        axes = crossSectionModule.defaultAxes( crossSectionUnit = 'b' )
    else:
        axes = axesModule.axes( labelsUnits = { 1 : ( 'energy_in' , 'eV' ), 0 : ( 'weight' , '' ) } )
    cls = crossSectionModule.XYs1d
    if( MF == 9 ) : cls = multiplicityModule.XYs1d
    for i1 in xrange( NS ) :
        QM, QI, IZAP, LFS, NR, NP = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MFData[0], intIndices = [ 0, 2, 4 ], logFile = info.logs )
        dataLine, TAB1, regions = endfFileToGNDMisc.getTAB1Regions( dataLine, MFData, axes = axes, logFile = info.logs )
        if( MF == 10 ) :
            XSec = getCrossSectionForm( info, regions )
        else :
            XSec = regions
        MF9or10.append( [ TAB1, XSec ] )
    return( MF9or10  )

def readMF12_13( info, MT, MTData, productList, warningList, gammaBRTolerance = 1e-6 ) :

    def addMF12_13GammaToList( gList, EGk, ESk, LP, LF, regions ) :
        """EGk is the gamma's energy and ESk is the gamma's origination level energy."""

        if( EGk in gList ) : raise Exception( 'Gammas with the same energy (%s) are not supported: MT=%s' % ( EGk, MT ) )
        gList.append( { 'EGk' : EGk, 'ESk' : ESk, 'LP' : LP, 'LF' : LF, 'yield' : regions, 'angularSubform' : None, 'energySubform' : None } )

    def checkAngularData( gamma ) :

        angularSubform = gamma['angularSubform']
        if( angularSubform is None ) : angularSubform = angularModule.isotropic( )
        gamma['angularSubform'] = angularSubform

    def addGammaProduct( info, MF, gamma, productList, warningList, ESk ) :

        yields = gamma['yield']
        multiplicity = getMultiplicityPointwiseOrPieceWise( info, yields, warningList )
        if( MF == 13 ) : multiplicity._temp_divideByCrossSection = True

        attributes = {}
        if( gamma['ESk'] != 0 ) : attributes['originationLevel'] = "%s" % PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( ESk ), "eV" ).toString( keepPeriod = False )
        product = toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, 0, None ), multiplicity = multiplicity, attributes = attributes )

        angularSubform = gamma['angularSubform']
        energySubform = gamma['energySubform']
        form = uncorrelated( info.style, frames[1], angularSubform, energySubform )
        product.distribution.add( form )
        if( MF == 13 ) : product.attributes['ENDFconversionFlag'] = 'MF13'
        productList.append( product )

    if( 12 in MTData ) :
        if( 13 in MTData ) :
            raise Exception( 'MF = 12 and 13 present for MT=%s, this is not supported' % MT )
        MF = 12
    elif( 13 in MTData ) :
        MF = 13
    elif( ( 14 in MTData ) or ( 15 in MTData ) ) :
        warningList.append('MF 14 and/or 15 data and no MF 12 or 13 data: MT=%s MFs=%s' % ( MT, MTData.keys( ) ) )
        info.doRaise.append( warningList[-1] )
        return
    else :
        return

    MF12_13Data = MTData[MF]
    ZA, AWR, LO, LG, NK, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF12_13Data[0], intIndices = [ 0, 2, 3, 4 ], logFile = info.logs )
    printAWR_mode( info, MT, MF, ZA, AWR )
    countZAMasses( info, ZA, AWR )
    info.logs.write( ' : MF=%s LO=%s : ZAP=0 ' % ( MF, LO ) )
    dataLine, continuousGamma, discreteGammas, primaryGammas, branchingGammas = 1, [], [], [], []
    if( ( ( MF == 12 ) and ( LO == 1 ) ) or ( ( MF == 13 ) and ( LO == 0 ) ) ) :
        if( MF == 12 ) :
            axes = multiplicityModule.XYs1d.defaultAxes( )
        else :
            axes = crossSectionModule.defaultAxes( crossSectionUnit = 'b' )
        if( NK > 1 ) :
            dataLine, TAB1, regions = endfFileToGNDMisc.getTAB1Regions( dataLine, MF12_13Data, axes = axes, logFile = info.logs )
            info.totalMF6_12_13Gammas[MT] = [ MF, getMultiplicityPointwiseOrPieceWise( info, regions, warningList ) ]
        for i in xrange( NK ) :
            dataLine, TAB1, regions = endfFileToGNDMisc.getTAB1Regions( dataLine, MF12_13Data, axes = axes, logFile = info.logs )
            EGk, ESk, LP, LF = TAB1['C1'], TAB1['C2'], int( TAB1['L1'] ), int( TAB1['L2'] )
            if( EGk == 0. ) :
                if( LP not in [ 0, 1 ] ) : raise Exception( 'LP = %s != 0 for continuous gamma for MT = %s' % ( LP, MT ) )
                if( len( continuousGamma ) == 0 ) :
                    addMF12_13GammaToList( continuousGamma, EGk, ESk, LP, LF, regions )
                else :
                    raise Exception( 'continuous gamma information for MF=%s, MT = %s already exist' % ( MF, MT ) )
            elif( LP == 2 ) :
                addMF12_13GammaToList( primaryGammas, EGk, ESk, LP, LF, regions )
            else :
                addMF12_13GammaToList( discreteGammas, EGk, ESk, LP, LF, regions )
    elif( ( MF == 12 ) and ( LO == 2 ) ) :
        dataLine, LO2 = endfFileToGNDMisc.getList( dataLine, MF12_13Data, logFile = info.logs )
        LP = int(LO2['L1'])
        if( ( LP == 2 ) and ( MT not in ( 91, 649, 699, 749, 799, 849 ) ) ) :
            warningList.append("Incorrect 'primary gamma' flag for MF12 MT%d" % MT)
        NT, LGp = LO2['N2'], LG + 1
        NK = NT
        for idx in xrange( NT ) :
            parentEnergy, finalEnergy = LO2['C1'], LO2['data'][idx*LGp]
            if parentEnergy==0:
                raise Exception("Gamma decay from ground state in MF12 MT%d" % MT)
            if abs((parentEnergy-finalEnergy)/parentEnergy)<0.001:
                raise Exception("Zero-energy gamma from %f eV to %f eV in MF12 MT%d" % (parentEnergy,finalEnergy,MT))
            branchingGammas.append( {'ES' : LO2['C1'], 'EGk' : 0, 'ESk' : LO2['data'][idx*LGp], 'angularSubform' : None,
                                     'branching' : LO2['data'][idx*LGp+1:idx*LGp+LGp]} )
        gammaBRList = [ g['branching'][0] for g in branchingGammas ]
        sumGammaBRList = sum( gammaBRList )
        if abs( sumGammaBRList - 1.0 ) > gammaBRTolerance:
            warningList.append( "sum of gamma BR's for MT="+str(MT)+" MF=12 is " + str(sumGammaBRList)+' != 1.0' )
            # leave re-normalization for checker/fixer codes
            #for i in xrange( len( branchingGammas ) ): branchingGammas[ i ][ 'branching' ][ 0 ] /= sumGammaBRList
        info.MF12_LO2[MT] = branchingGammas
    else :
        raise Exception( 'LO=%s is not valid for MF=%s, MT=%s' % ( LO, MF, MT ) )

    readMF14( info, MT, MTData, MF, NK, warningList, discreteGammas, primaryGammas, continuousGamma, branchingGammas )
    readMF15( info, MT, MTData, continuousGamma, warningList )

    for gamma in branchingGammas :
        if( not( gamma['angularSubform'] is None ) ) :
            info.logs.write( 'NON-isotropic gamma' )
            break
    if( 14 in MTData ) : info.logs.write( ': MF=14 ' )
    if( 15 in MTData ) : info.logs.write( ': MF=15 ' )

    if( len( continuousGamma ) ) :
        gamma = continuousGamma[0]
        checkAngularData( gamma )
        addGammaProduct( info, MF, gamma, productList, warningList, gamma['ESk'] )

    for gamma in discreteGammas :
        checkAngularData( gamma )
        gamma['energySubform'] = energyModule.constant( PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( gamma['EGk'] ), 'eV' ) )
        addGammaProduct( info, MF, gamma, productList, warningList, gamma['ESk'] )
    for gamma in primaryGammas :
        checkAngularData( gamma )
        gamma['energySubform'] = energyModule.primaryGamma( PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( gamma['EGk'] ), 'eV' ) )
        addGammaProduct( info, MF, gamma, productList, warningList, gamma['ESk'] )
    for gamma in branchingGammas : checkAngularData( gamma )

def readMF14( info, MT, MTData, MF, NK, warningList, discreteGammas, primaryGammas, continuousGamma, branchingGammas ) :

    allGammas = discreteGammas + primaryGammas + continuousGamma + branchingGammas

    def addAngularData( EGk, ESk, subform ) :

        possibleGammas = []
        for gamma in allGammas :
            if( ( abs( EGk - gamma['EGk'] ) < 1e-4 ) and ( abs( ESk - gamma['ESk'] ) < 1e-4 ) ) : possibleGammas.append( gamma )
        if( len( possibleGammas ) == 0 ) :
            sys.stderr.write( "MF = 14 gamma data: EGk = %s, ESk = %s\n" % ( EGk, ESk ) )
            sys.stderr.write( 'MF = %s gamma list is,\n' % MF )
            for gamma in allGammas : sys.stderr.write( "    EGk = %s, ESk = %s\n" % ( gamma['EGk'], gamma['ESk'] ) )
            raise Exception( 'No matching gamma from MF = %s found for MF = 14 gamma: MT=%s: see above' % ( MF, MT ) )
        if( len( possibleGammas ) != 1 ) : raise Exception( 'Multiple possible gammas found for MF=14 data: EGk=%s, ESk=%s, MT=%s' % ( EGk, ESk, MT ) )
        gamma = possibleGammas[0]
        if( gamma['angularSubform'] is None ) :
            gamma['angularSubform'] = subform
        else :
            raise Exception( 'Gamma already has MF=14 angular data: EGk=%s, ESk=%s, MT=%s' % ( EGk, ESk, MT ) )

    if( not( 14 in MTData ) ) : return
    MF14Data = MTData[14]
    ZA, AWR, LI, LTT, NK14, NI = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF14Data[0], intIndices = [ 0, 2, 3, 4, 5 ], logFile = info.logs )
    printAWR_mode( info, MT, 14, ZA, AWR )
    countZAMasses( info, ZA, AWR )
    if( ( NK14 != NK ) and info.printBadNK14 ) :
        warningList.append( 'MF14 NK = %s != MF12/13 NK = %s for MT = %s' % ( NK14, NK, MT ) )
    dataLine, frame = NI + 1, standardsModule.frames.labToken
    if( LTT == 0 ) :                                     # All distributions are isotropic
        pass
    elif( LTT == 1 ) :
        for i in xrange( NK14 - NI ) :
            dataLine, angularData = endfFileToGNDMisc.getTAB2_Lists( dataLine, MF14Data, logFile = info.logs )
            if( angularData['NR'] != 1 ) : raise Exception( 'Currently only one interpolation flag is supported: NR=%s, MT=%s' % ( angularData['NR'], MT ) )
            EGk, ESk = angularData['C1'], angularData['C2']
            subform = angularLegendreToPointwiseOrPiecewiseLegendre( MT, angularData, warningList, 4, 'LTT = 1' )
            addAngularData( EGk, ESk, subform )
    else :
        raise Exception( 'MF=14, LI=%s not implemented' % LI )

def readMF15( info, MT, MTData, continuousGamma, warningList ) :

    if( 15 not in MTData ) :
        if( len( continuousGamma ) ) :
            warningList.append( 'Continous gamma with no MF=15 data: MT=%s' % MT )
            info.doRaise.append( warningList[-1] )
        return
    if( len( continuousGamma ) == 0 ) :
        warningList.append( 'MF=15 data and no continous gamma MF=12,13 data: MT=%s' % MT )
        info.doRaise.append( warningList[-1] )
    MF15Data = MTData[15]
    ZA, AWR, dummy, dummy, NC, dummy = endfFileToGNDMisc.sixFunkyFloatStringsToIntsAndFloats( MF15Data[0], intIndices = [ 0, 4 ], logFile = info.logs )
    printAWR_mode( info, MT, 15, ZA, AWR )
    countZAMasses( info, ZA, AWR )
    if( NC != 1 ) : raise Exception( 'NC = %s > 1 is not supported: MT=%s' % ( NC, MT ) )

    dataLine, weights = endfFileToGNDMisc.getTAB1( 1, MF15Data, logFile = info.logs )
    for Ein, weight in weights['data'] :
        if( weight != 1 ) : raise Exception( 'For MF15 data weight other than 1 is not currently supportd: MT%d' % MT )

    LF = int( weights['L2'] )
    if( LF != 1 ) : raise Exception( 'For MF=15 data, only LF=1 is currently supported, not LF=%s : MT%s' % ( LF, MT ) )
    axes2d = energyModule.XYs2d.defaultAxes( )
    axes = axesModule.axes( )
    axes[0] = axes2d[0]
    axes[1] = axes2d[1]
    dataLine, EEpETable = endfFileToGNDMisc.getTAB2_TAB1s( dataLine, MF15Data, logFile = info.logs, axes = axes )
    continuousGamma[0]['energySubform'] = toPointwiseOrPiecewiseEnergy( MT, EEpETable )

def genID( cov_info, MT, MF, MT2=None, MF2=None, MAT2=None, QI=None, QI2=None, linkType='rowColumn' ):
    """
    For covariances we need a unique id for each section, and also need link info.
    This is messy: lots of special cases.
    """

    evalStyle = cov_info['style']
    MTdict = cov_info['MTdict']
    def getReaction( MT, MF, QI=None ):
        # find the <reaction> corresponding to this covariance info, if any
        if MT in (452,455,456):
            MT = 18 # nubar actually refers to fission
        if MT in MTdict:
            chThisMT = MTdict[MT]
            if len(chThisMT)==1:
                return chThisMT[0]
            elif MF==40:   # MF40 section may contain covariances for multiple excited states of residual
                thisChannel = [ch for ch in chThisMT if isinstance(ch, productionModule.production) and ch.getQ('eV')==QI ]
                if len(thisChannel)==1: return thisChannel[0]
                elif MT in cov_info:
                    # Q-values disagree between MF10/MF40. Assume 1st covariance is for lowest-energy residual, etc
                    index = cov_info[MT]
                    cov_info[MT] += 1
                else:
                    index = 0; cov_info[MT] = 1
                return thisChannel[ index ]
            elif MF==33:
                # residual must be in ground state unless MT in 51-90, etc
                thisChannel = [ch for ch in chThisMT if (isinstance(ch, reactionModule.reaction) and
                        sum( [p.getLevelAsFloat('eV') for p in ch.outputChannel] )==0) or
                        isinstance(ch, sumsModule.crossSectionSum)]
                if len(thisChannel) != 1: raise BadCovariance("MF33 MT%d covariance doesn't correspond to any channel!"
                        % MT )
                return thisChannel[0]
            else :
                raise BadCovariance( "Can't determine which reaction this covariance (MF%d MT%d) corresponds to!" % ( MF, MT ) )
        return

    rowReaction = getReaction( MT, MF, QI )

    def makeID(MT, reaction):
        if reaction is not None:
            #Id = str( reaction.outputChannel )
            Id = str( reaction )
        elif MT in (452,455,456):
            Id = { 452 : 'total', 455 : tokensModule.delayedToken, 456 : tokensModule.promptToken }[MT]
        elif MT in (1,4,103,104,105,106,107):
            Id = {1:'total', 4:'sum(z,n)', 103:'sum(z,p)', 104:'sum(z,d)',
                    105:'sum(z,t)', 106:'sum(z,He3)', 107:'sum(z,a)'}[MT]
        elif 850 < MT < 871:
            Id = "lump%d" % (MT-851)
        else: Id = "MF%d_MT%d" % (MF,MT)
        return Id

    rowId = makeID( MT, rowReaction )

    if MAT2:
        # cross-material covariance
        colReaction = None  # must create an 'externalReaction' and link to it below
        versus = ' vs. '
        colId = "MAT%d_MF%d_MT%d" % (MAT2,MF2,MT2)
    elif (MT2 and MT2!=MT):
        MF2 = MF2 or MF
        colReaction = getReaction( MT2, MF2, QI2 )
        versus = ' vs. '
        colId = makeID( MT2, colReaction )
    else:
        colId = versus = ''
    qualifier = {31:' [nubar]', 33:'', 34:' [angular distribution]', 35:' [spectrum]', 40:''}[MF]
    ID = rowId + versus + colId + qualifier

    """ also, create links from covariances to data. If we're trying to point to nonexistant data
    (ie for lumped channels), make a 'fake channel' """
    filename = ''   # should hold file for covariances
    def makeLink( MT, reaction, Id, MAT2, linkClass ):
        link_ = linkClass(ENDF_MFMT="%d,%d" % (MF,MT))
        if reaction is not None:
            if MF in (33,40):
                link_.link = reaction.crossSection[evalStyle]
                if isinstance(link_.link, crossSectionModule.resonancesWithBackground):
                    link_.link = link_.link.tabulatedData
            elif MF==31: link_.link = reaction.outputChannel[0].multiplicity[evalStyle]
            elif MF==34:
                link_.link = reaction.outputChannel[0].distribution[evalStyle].subforms[0]
            elif MF==35:
                distribution = reaction.outputChannel[0].distribution[evalStyle]
                if( isinstance( distribution, uncorrelatedModule.form ) ) :
                    link_.link = distribution.energySubform
                else :
                    link_.link = distribution
        elif MAT2:
            # cross-material covariance: get other isotope from the MAT number
            try: otherTarget = endf_endl.getParticleNameFromMAT( MAT2 )
            except KeyError:
                raise BadCovariance( "Encountered illegal MAT number %d in covariances!" % MAT2 )
            quant = covarianceSectionModule.externalReaction(id=Id, target=otherTarget, ENDF_MFMT=(MF,MT))
            cov_info['externalReactions'][(MAT2,MT,MF)] = quant
            link_.link = quant
        else:
            if 850 < MT < 871: pass
            elif MT==4: cov_info['MTL_2'][(MT,MF)] = zip(range(50,92),[33]*41)
            elif MT==103: cov_info['MTL_2'][(MT,MF)] = zip(range(600,650),[33]*49)
            elif MT==104: cov_info['MTL_2'][(MT,MF)] = zip(range(651,700),[33]*49)
            elif MT==105: cov_info['MTL_2'][(MT,MF)] = zip(range(701,750),[33]*49)
            elif MT==106: cov_info['MTL_2'][(MT,MF)] = zip(range(751,800),[33]*49)
            elif MT==107: cov_info['MTL_2'][(MT,MF)] = zip(range(801,850),[33]*49)
            elif MT==1: cov_info['MTL_2'][(MT,MF)] = zip(range(2,1000),[33]*998)
            if (MT,MF) not in cov_info['lumpedChannels']:
                quant = covarianceSectionModule.reactionSum(id=Id, ENDF_MFMT=(MF,MT))
                cov_info['lumpedChannels'][(MT,MF)] = quant
            quant = cov_info['lumpedChannels'][(MT,MF)]

            link_.link = quant

        return link_

    linkClass = covarianceSummedModule.summand
    if linkType == 'rowColumn':
        linkClass = covarianceSectionModule.rowData
    rowData = makeLink( MT, rowReaction, rowId, MAT2, linkClass )
    colData = None
    if (MT2 and MT2!=MT) or MAT2:
        if linkType == 'rowColumn':
            linkClass = covarianceSectionModule.columnData
        colData = makeLink( MT2, colReaction, colId, MAT2, linkClass )
#        colData.label = 'columnData'           # CALEB, what is this for?
    return ID, [rowData, colData]

def readMatrix( info, LS,LB, NT,NP, dat ):
    """ matrix format is very similar for MF31, MF33 and MF35, so we can
    generalize parsing the matrix """

    if NP in (0,1):   # matrix has at most one energy boundary, no data
        return None
    nlines, remainder = divmod(NT,6)
    subsec = []
    for line in range(nlines): subsec += funkyF(dat.next(), logFile = info.logs)
    if remainder: subsec += funkyF(dat.next(), logFile = info.logs)[:remainder]
    # LB flag tells how to interpret this data:
    if LB in (0,1,2,8): # diagonal
        subsec = subsec[:NT]
        energyBounds = [subsec[::2]]
        data = subsec[1::2]
        matrix = arrayModule.diagonal( (len(data),len(data)), data )
    elif LB==5:
        if LS==1: #symmetric upper-diagonal
            energyBounds = [subsec[:NP],]
            data = linearAlgebraModule.switchSymmetry( subsec[NP:NT] )
            matrix = arrayModule.full( (NP-1,NP-1), data, symmetry=arrayModule.symmetryLowerToken )
        else: # asymmetric, but same energy grids for both axes:
            energyBounds = [subsec[:NP]]
            matrix = subsec[NP:]
            matrix = arrayModule.full( (NP-1,NP-1), matrix )
    elif LB==6: # asymmetric
        NER = NP
        NEC = (NT-1)//NER
        energyBounds = [subsec[:NER], subsec[NER:NER+NEC]]
        matrix = subsec[NER+NEC:NT]
        NEC-=1; NER-=1  # matrix dimensions
        matrix = arrayModule.full( (NER,NEC), matrix )
    else:
        return None

    # FIXME: covariances need 'flat' interpolation along both independent axes, but gridded doesn't support that yet
    axes = axesModule.axes( labelsUnits={0:('matrix_elements',''), 1:('column_energy_bounds','eV'), 2:('row_energy_bounds','eV')} )
    axes[2] = axesModule.grid( axes[2].label, axes[2].index, axes[2].unit,
                style = axesModule.boundariesGridToken, values = valuesModule.values( energyBounds[0] ) )
    if len(energyBounds)==2:
        axes[1] = axesModule.grid( axes[1].label, axes[1].index, axes[1].unit,
                style = axesModule.boundariesGridToken, values = valuesModule.values( energyBounds[1] ) )
    else:
        axes[1] = axesModule.grid( axes[1].label, axes[1].index, axes[1].unit,
                style = axesModule.linkGridToken, values = linkModule.link( link = axes[2].values, relative = True ) )
    return griddedModule.gridded( axes, matrix )

def readMF31_33( info, dat, mf, mt, cov_info, warningList ):
    """ nubar and cross section covariances have basically the same form
    in ENDF, so we can treat them the same way: """

    dat = myIter(dat)
    # dat contains one MFMT section, but in gnd each cross-correlation gets its own section:
    sectionList, linkData = [], []

    ZA,AWR,dum,MTL,dum,NL = funkyFI(dat.next(), logFile = info.logs)
    if MTL!=0:
        # MTL!=0 implies section is a placeholder pointing to a lumped channel
        if (MTL,mf) in cov_info['MTL']:
            cov_info['MTL'][(MTL,mf)].append((mt,mf))
        else: cov_info['MTL'][(MTL,mf)] = [(mt,mf)]
    for subsection in range(NL):
        XMF1,XLFS1,MAT1,MT1,NC,NI = funkyFI(dat.next(), logFile = info.logs)
        if MAT1 == info.MAT:
            # Some files explicitly give MAT1==MAT for internal covariance.
            # for simplicity, just set MAT1=0 unless it is actually a different material
            MAT1 = 0
        XMF1,XLFS1 = int(XMF1),int(XLFS1)

        covarsThisSection = []

        if XMF1!=0 or XLFS1!=0:
            #info.logs.write( "non-zero XMF1/XLFS1!!!" ) # may not be properly dealt with
            raise BadCovariance( "non-zero XMF1/XLFS1 in covariances not currently handled!" )

        for NCdx in range(NC):
            dum,dum,dum,LTY,dum,dum = funkyFI(dat.next(), logFile = info.logs)
            if LTY==0:
                E1,E2,dum,dum,NCI2,NCI = funkyFI(dat.next(), logFile = info.logs)
                subsec = []
                nlines = int(math.ceil(NCI2/6.0))
                for line in range(nlines): subsec += funkyF(dat.next(), logFile = info.logs)
                #coefs = subsec[:NCI2][::2]
                pointerList = []
                for coef,mtnum in zip(subsec[:NCI2][::2],subsec[:NCI2][1::2]):
                    Id, pointers = genID( cov_info, int(mtnum), mf, linkType='summand' )
                    link = pointers[0]
                    link.attributes['coefficient'] = coef
                    pointerList.append( link )
                covarsThisSection.append( covarianceSummedModule.summedCovariance(
                        label = info.style,
                        lowerBound = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( E1 ),'eV'),
                        upperBound = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( E2 ),'eV'),
                        pointerList = pointerList) )
                cov_info['NC_data'].append( covarsThisSection[-1] )
            else:
                warningList.append( 'non-zero LTY in MF33' )

        for NIdx in range(NI):
            dum,dum,LS,LB,NT,NP = funkyFI(dat.next(), logFile = info.logs)
            matrix = readMatrix( info, LS,LB,NT,NP, dat )
            if LB not in (0,1,2,5,6,8):
                warningList.append( 'skipping LB%d section for MF%d MT%d' % ( LB, mf, mt ) )
                continue
            if matrix is None:
                warningList.append( 'skipping empty matrix for MF%d MT%d' % ( mf, mt ) )
                info.doRaise.append( warningList[-1] )
                continue
            Type='relative'
            if LB==0:
                Type='absolute'
                matrix.axes[0].unit = 'b**2'

            ENDFconversionFlag = None
            if LB in (2,3,4,8,): ENDFconversionFlag = "LB=%d" % LB

            covarsThisSection.append( covarianceBaseModule.covarianceMatrix( label = info.style, type=Type, matrix=matrix,
                ENDFconversionFlag=ENDFconversionFlag) )

        # create unique id for each section:
        idNow, pointers = genID( cov_info, mt, mf, MT2=MT1, MF2=(XMF1 or mf), MAT2=MAT1 )
        rowdat, coldat = pointers
        section = covarianceSectionModule.section( id=idNow, rowData=rowdat, columnData=coldat )

        if len(covarsThisSection)>1:
            form = covarianceMixedModule.mixedForm( label = info.style, components=covarsThisSection )
            for idx in range( len( form ) ) :
                form[idx].label = str(idx)
        elif len(covarsThisSection)==1:
            form = covarsThisSection[0]
        else:
            #raise Exception("Encountered empty covariance section!!!")
            info.logs.write("Missing covariance data from section!")
            continue
        section.add( form )

        sectionList.append( section )
        linkData.append( (mt,mf,MT1,XMF1, idNow) )
        # end loop over NL subsections

    if dat.index != dat.length: raise BadCovariance("Not all covariance data converted, MF%d MT%d" % (mf,mt))
    return sectionList, linkData

def readMF32( info, dat, mf, mt, cov_info, warningList ) :
    # MF=32 resonance parameter covariances. Must be synchronized with MF=2 resonance parameters

    try: import numpy
    except ImportError:
        warningList.append("Skipping MF32 since NumPy is unavailable")
        return [],[]
    resonances = cov_info['resonances']

    def swaprows( matrix, i1, i2, nrows ):
        # matrix rows may be out-of-order and need resorting
        rows = matrix[i1:i1+nrows].copy()
        matrix[i1:i1+nrows] = matrix[i2:i2+nrows]; matrix[i2:i2+nrows] = rows
        cols = matrix[:,i1:i1+nrows].copy()
        matrix[:,i1:i1+nrows] = matrix[:,i2:i2+nrows]; matrix[:,i2:i2+nrows] = cols

    dat = myIter(dat)
    ZA,AWR,dum,dum,NIS,dum = funkyFI(dat.next(), logFile = info.logs)
    if (NIS!=1): raise BadCovariance( "Can't handle multi-isotope file 32!" )
    ZAI,ABN,dum,LFW,NER,dum = funkyFI(dat.next(), logFile = info.logs)

    sections = []
    for subsection in range(NER):
        EL,EH,LRU,LRF,NRO,NAPS = funkyFI(dat.next(), logFile = info.logs)
        if (NRO!=0): raise BadCovariance( "Can't handle non-zero NRO in MF32!" )
        # format is determined mainly by combination of LCOMP and LRU/LRF
        if LRU==1:  # resolved resonance covariance section
            attributes = {'endfConversionFlags':''}
            if LRF==7:
                warningList.append("skipping LRF=7 resonance covariances!")
                continue
            if LRF in (1,2,3):  # Breit-Wigner and simplified Reich-Moore formats are similar
                mf2_elist = zip( resonances.resolved.evaluated.resonanceParameters.table.getColumn('energy'),
                        resonances.resolved.evaluated.resonanceParameters.table.getColumn('neutronWidth'),
                        resonances.resolved.evaluated.resonanceParameters.table.getColumn('captureWidth') )
                SPI,AP,dum,LCOMP,NLS,ISR = funkyFI(dat.next(), logFile = info.logs)
                DAP = []
                if ISR>0:   # scattering radius uncertainty
                    if LRF in (1,2):
                        dum,dap,dum,dum,dum,dum = funkyFI(dat.next(), logFile = info.logs)
                        DAP = [10*dap]
                    elif LRF==3:
                        dum,dum,dum,dum,MLS,one = funkyFI(dat.next(), logFile = info.logs)
                        for idx in range( int( math.ceil(MLS/6.0) ) ):
                            DAP.extend( funkyF(dat.next(), logFile = info.logs) )
                        DAP = [10*dap for dap in DAP[:MLS]] # convert from 10*fm to fm
                    else:
                        raise BadCovariance("ISR>0 not yet supported for LRF=%d!" % LRF)
                if LCOMP==0:
                    # internal correlations given for each resonance, no cross-resonance terms
                    attributes['endfConversionFlags'] = 'LCOMP=0'
                    mf32_resonances, mf32_covars = [],[]
                    NRS = 0
                    for Lval in range(NLS):
                        AWRI, dum, L, dum, tmp, nrs_ = funkyFI(dat.next(), logFile = info.logs)
                        NRS += nrs_
                        for i in range(nrs_):
                            mf32_resonances.append( funkyF(dat.next(), logFile = info.logs) )
                            mf32_covars.append( funkyF(dat.next(), logFile = info.logs) + funkyF(dat.next(), logFile = info.logs) )

                    dEsq, dNsq, dNdG, dGsq, dNdF, dGdF, dFsq, dJdN, dJdG, dJdF, dJsq, dum = zip(*mf32_covars)
                    MPAR = 3
                    if any(dFsq): MPAR = 4
                    if any(dJsq): raise BadCovariance("Encountered uncertainty on J in MF32!")
                    matrix = numpy.zeros((MPAR*NRS,MPAR*NRS))
                    for ridx in range(NRS):
                        matrix[ridx*MPAR,ridx*MPAR] = dEsq[ridx]
                        matrix[ridx*MPAR+1,ridx*MPAR+1] = dNsq[ridx]
                        matrix[ridx*MPAR+2,ridx*MPAR+1:ridx*MPAR+3] = [dNdG[ridx], dGsq[ridx]]
                        if MPAR==4:
                            matrix[ridx*MPAR+3,ridx*MPAR+1:ridx*MPAR+4] = [dNdF[ridx], dGdF[ridx], dFsq[ridx]]
                    # symmetrize:
                    for ridx in range( MPAR*NRS ):
                        matrix[ridx,ridx:] = matrix[ridx:,ridx]

                    mf32_elist = [(lis[0],lis[3],lis[4]) for lis in mf32_resonances]
                    Type="absolute"
                    GNDmatrix = gndMatrix.matrix( matrix, form="sparse_symmetric" )

                if LCOMP==1:
                    AWRI,dum,dum,dum,NSRS,NLRS = funkyFI(dat.next(), logFile = info.logs)
                    dum,dum,MPAR,dum,tmp,NRB = funkyFI(dat.next(), logFile = info.logs)
                    dim = NRB * MPAR  # num. of resonances * num. parameters per resonance
                    matrixSize = dim * (dim+1) // 2
                    if matrixSize + 6*NRB != tmp:
                        raise BadCovariance("Incorrect dimension for the matrix!")

                    # resonanances are listed again (redundant!):
                    mf32_resonances = [ funkyF(dat.next(), logFile = info.logs) for i in range(NRB) ]
                    if LRF in (1,2):
                        mf32_elist = [(lis[0],lis[3],lis[4]) for lis in mf32_resonances]
                    else:
                        mf32_elist = [(lis[0],lis[2],lis[3]) for lis in mf32_resonances]
                    # then the matrix:
                    data = []
                    nLines, rem = divmod( matrixSize, 6 )
                    for i in range(nLines): data += funkyF(dat.next(), logFile = info.logs)
                    if rem: data += funkyF(dat.next(), logFile = info.logs)[:rem]
                    matrix = numpy.zeros((dim,dim))
                    start, length = 0, dim
                    for i in range(dim):
                        # data stores upper-diagonal matrix. Symmetrize:
                        matrix[i,i:] = matrix[i:,i] = data[start:start+length]
                        start = start+length; length = length-1
                    Type="absolute"
                    GNDmatrix = gndMatrix.matrix( matrix, form="symmetric" )

                elif LCOMP==2:
                    attributes['endfConversionFlags'] = 'LCOMP=2'
                    AWRI, QX, dum, LRX, tmp, NRSA = funkyFI(dat.next(), logFile = info.logs)
                    # resonance parameters + uncertainties:
                    mf32_resonances = [funkyF(dat.next(), logFile = info.logs) for i in range(NRSA*2)]
                    if LRF in (1,2):
                        mf32_elist = [(lis[0],lis[3],lis[4]) for lis in mf32_resonances[::2]]
                    else:
                        mf32_elist = [(lis[0],lis[2],lis[3]) for lis in mf32_resonances[::2]]
                    # for LCOMP==2, off-diagonal terms are given as correlation matrix:
                    dum,dum,NDIGIT,NNN,NM,dum = funkyFI(dat.next(), logFile = info.logs)
                    MPAR = NNN/NRSA
                    matrix = numpy.eye(NNN)
                    diagonal = []
                    for idx in range(NRSA):
                        if LRF in (1,2):
                            dE,dum,dum,dGammaN,dGammaG,dGammaF = mf32_resonances[2*idx+1]
                            diagonal.extend( [dE,dGammaN,dGammaG] )
                            if MPAR==4: diagonal.append( dGammaF )
                        elif LRF==3:
                            dE,dum,dGammaN,dGammaG,dGammaF1,dGammaF2 = mf32_resonances[2*idx+1]
                            diagonal.extend( [dE,dGammaN,dGammaG] )
                            if MPAR==4: diagonal.extend( [dGammaF1] )
                            elif MPAR==5: diagonal.extend( [dGammaF1,dGammaF2] )
                    if len(diagonal)!=NNN: raise BadCovariance( "Incorrect dimensions for LCOMP=2 matrix!" )

                    # off-diagonal parts of matrix are stored as sparse correlation matrix:
                    attributes['endfConversionFlags'] += ',NDIGIT=%d' % NDIGIT
                    matrix = numpy.eye(NNN, dtype="float") * 10**NDIGIT
                    for idx in range(NM):
                        row,col,vals = endfFileToGNDMisc.readEndfINTG( dat.next(), NDIGIT )
                        vals = vals[:row-col]
                        # go to 0-based index:
                        row -= 1
                        col -= 1
                        if row>=NNN or col>=row:
                            raise BadCovariance("Matrix indices out of range for MF32 LCOMP=2 matrix")
                        matrix[row,col:col+len(vals)] = vals
                    for idx in range(NNN):    # symmetrize
                        matrix[idx,idx:] = matrix[idx:,idx]
                    matrix /= float(10**NDIGIT)

                    # convert correlation -> relative covariance matrix
                    rsd = numpy.sqrt( numpy.array( diagonal ) )
                    matrix = matrix * numpy.outer(rsd,rsd)
                    Type="absolute"
                    GNDmatrix = gndMatrix.matrix( matrix, form="sparse_symmetric" )

                # mf32 may cover only the low-energy portion of the resonance region:
                if len(mf2_elist) != len(mf32_elist):
                    mf2_elist = mf2_elist[:len(mf32_elist)]

                # check whether resonances are sorted by L rather than by energy. If so, rearrange:
                if mf32_elist != mf2_elist or LCOMP==0:
                    if attributes['endfConversionFlags'] == '': attributes['endfConversionFlags'] = 'sortByL'
                    else: attributes['endfConversionFlags'] += ',sortByL'

                mf32_elist_sorted = sorted(mf32_elist, key=lambda res: res[0])
                if mf32_elist_sorted!=mf2_elist:
                    ndiffs = sum( [tmpa != tmpb for (tmpa,tmpb) in zip(mf2_elist, mf32_elist_sorted)] )
                    warningList.append( "MF2/MF32 resonance parameters differ for %d resonances. For example:" % ndiffs )
                    for idx, (mf2res,mf32res) in enumerate(zip(mf2_elist, mf32_elist_sorted)):
                        if mf2res != mf32res:
                            warningList.append( "    resonance #%d: MF2 = %s, MF32 = %s"
                                    % (idx, mf2res, mf32res) )
                            if not info.verboseWarnings: break

                    raise BadCovariance("MF32 resonance parameters don't match MF2 parameters!")

                for i1 in range(len(mf2_elist)):
                    i2 = mf32_elist.index( mf2_elist[i1] )
                    if i2 != i1:
                        swaprows( GNDmatrix.data, MPAR*i1, MPAR*mf32_elist.index( mf2_elist[i1] ), MPAR )
                        # also need to swap values in elist2:
                        val = mf32_elist[i1]
                        mf32_elist[i1] = mf32_elist[i2]; mf32_elist[i2] = val

                if DAP: # scattering radius uncertainty was specified. Expand matrix to include it:
                    if len(DAP) > 1:
                        raise BadCovariance("Energy-dependent scat. radius uncertainty not yet handled!")
                    dim = len(matrix) + len(DAP)
                    new_matrix = numpy.zeros( (dim,dim) )
                    for i in range(len(DAP)): new_matrix[i,i] = DAP[i]
                    new_matrix[ len(DAP):, len(DAP): ] = matrix
                    GNDmatrix = gndMatrix.matrix( new_matrix, form=GNDmatrix.form,
                            precision=GNDmatrix.precision )

                # switch to diagonal matrix if possible (much more compact):
                if numpy.all( matrix==( numpy.identity(len(matrix)) * matrix.diagonal() ) ):
                    GNDmatrix.form = "diagonal"

                # store into GND:
                resData = resonances.resolved.evaluated.resonanceParameters.table
                nResonances = len( mf32_elist )
                parametersPerResonance = "energy,neutronWidth,captureWidth"
                if MPAR>=4: parametersPerResonance = parametersPerResonance + ',fissionWidthA'
                if MPAR==5: parametersPerResonance = parametersPerResonance + ',fissionWidthB'
                inputParameters = [ covarianceModelParametersModule.loopOverResonanceParameters( link=resData,
                            nResonances=nResonances, parametersPerResonance=parametersPerResonance ) ]
                if DAP:
                    inputParameters.insert( 0, covarianceModelParametersModule.inputParameter("scatteringRadius",
                        resonances.resolved.evaluated.scatteringRadius, unit="fm") )
                sections.append( covarianceModelParametersModule.resonanceParameterCovariance( None, inputParameters,
                    GNDmatrix, type=Type, **attributes ) )
        else:
            # unresolved resonance parameters
            SPI,AP,dum,dum,NLS,dum = funkyFI( dat.next(), logFile = info.logs )
            for lval in range(NLS):
                AWRI,dum,L,dum,tmp,NJS = funkyFI( dat.next(), logFile = info.logs )
                if tmp!=6*NJS: raise BadCovariance( "Incorrect header in MF32 unresolved section!" )
                for jval in range(NJS):
                    D,AJ,GNO,GG,GF,GX = funkyFI( dat.next(), logFile = info.logs )

            # matrix:
            dum,dum,MPAR,dum,tmp,NPAR = funkyFI( dat.next(), logFile = info.logs )
            if tmp != (NPAR*(NPAR+1))/2: raise BadCovariance( "Incorrect header in MF32 unresolved section!" )
            nlines = int(math.ceil(tmp/6.0))
            data = []
            for line in range(nlines): data += funkyF(dat.next(), logFile = info.logs)
            matrix = numpy.zeros((NPAR,NPAR))
            start, length = 0, NPAR
            for i1 in range(NPAR):
                matrix[i1,i1:] = matrix[i1:,i1] = data[start:start+length]
                start = start+length; length = length-1
            if numpy.all( matrix==0 ):
                warningList.append("ignoring empty unresolved covariance matrix!")
            else:
                warningList.append("covariances for unresolved resonance region not yet supported!")

    return sections, []

def readMF34( info, dat, mf, mt, cov_info, warningList ):
    """ angular distribution covariances: """

    # dat contains one MFMT section
    dat = myIter(dat)
    sectionList, linkData = [], []

    ZA,AWR,dum,LTT,dum,NMT = funkyFI(dat.next(), logFile = info.logs)

    form = covarianceDistributionsModule.LegendreOrderCovarianceForm( label = info.style )
    groupIndex = 0

    for subsection in range(NMT):
        dum,dum,MAT1,MT1,NL,NL1 = funkyFI(dat.next(), logFile = info.logs)
        if MT1 == mt:
            NSS = NL*(NL+1)/2
        else:
            NSS=NL*NL1
        if MAT1 != 0: raise BadCovariance( "Cross material angular distribution covariance is not allowed in ENDF format.  Found MAT1 =",str(MAT1) )
        if MT1 != mt: raise NotImplementedError( "Cross reaction covariances in angular distribution covariance data not supported" )

        for iNSS in range(NSS):
            dum, dum, L, L1, LCT, NI = funkyFI(dat.next(), logFile = info.logs)
            frame = [ "frameOfMF4", standardsModule.frames.labToken, standardsModule.frames.centerOfMassToken][LCT]

            Lsection = covarianceDistributionsModule.LegendreLValue( L, L1, frame )
            covarsThisL = []

            for NIdx in range(NI):
                dum,dum,LS,LB,NT,NE = funkyFI(dat.next(), logFile = info.logs)
                if( ( LS not in [ 1, 0 ] ) and ( LB < 7 ) ) :
                    #info.logs.write( "Unexpected LS%d LB%d in MF35" % LB )
                    raise BadCovariance( "Unexpected LS%d LB%d in MF34" % LB )
                # pretend it's MF33 file:
                matrix = readMatrix( info, LS,LB,NT,NE, dat )
                if LB==0: raise Exception("LB=0 in Legendre covariances")
                covarsThisL.append( covarianceBaseModule.covarianceMatrix(
                    label = info.style, type='relative', matrix=matrix ) )

            if len(covarsThisL)>1:
                sectionForm = covarianceMixedModule.mixedForm( label = info.style, components=covarsThisL )
                for idx in range(len(sectionForm)): sectionForm[idx].label = str(idx)
            elif len(covarsThisL)==1:
                sectionForm = covarsThisL[0]

            Lsection.add( sectionForm )

            # end loop over subsubsections
            form.lvalues.append( Lsection )

            if frame == 'frameOfMF4':
                form.endfConversionFlag = "LCT=0"
                cov_info.setdefault( 'MF34_missingFrames', {} ).setdefault( mt, [] ).append( Lsection )

        if dat.index != dat.length: raise BadCovariance("Not all covariance data converted, MF%d MT%d" % (mf,mt))

        # add unique id to the section:
        idNow, pointers = genID( cov_info, mt, mf )
        rowdat, coldat = pointers
        section = covarianceSectionModule.section( id=idNow, rowData=rowdat, columnData=coldat)
        section.add( form )

        sectionList.append( section )
        linkData.append( (mt,mf,mt,mf, idNow) )
    return sectionList, linkData

def readMF35( info, dat, mf, mt, cov_info, warningList ):
    """ spectra covariances are fairly simple: """

    # dat contains one MFMT section
    dat = myIter(dat)
    sectionList, linkData = [], []

    ZA,AWR,dum,dum,NK,dum = funkyFI(dat.next(), logFile = info.logs)

    for subsection in range(NK):
        # each subsection contains matrix for a different incident energy interval.
        E1,E2,LS,LB,NT,NE = funkyFI(dat.next(), logFile = info.logs)
        if not (LS==1 and LB==7):
            #info.logs.write( "Unexpected LS%d LB%d in MF35" % LB )
            raise BadCovariance( "Unexpected LS%d LB%d in MF35" % LB )
        # pretend it's MF33 file:
        LS=1; LB=5
        E1,E2 = [PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( a ),'eV') for a in (E1,E2)]

        matrix = readMatrix( info, LS, LB, NT, NE, dat )
        form = covarianceBaseModule.covarianceMatrix( label = info.style, type='relative', matrix=matrix )

        # add unique id to the section:
        idNow, pointers = genID( cov_info, mt, mf )
        rowdat, coldat = pointers
        rowdat.attributes['incidentEnergyLowerBound'] = E1
        rowdat.attributes['incidentEnergyUpperBound'] = E2
        section = covarianceSectionModule.section( id=idNow, rowData=rowdat, columnData=coldat)
        section.add( form )

        sectionList.append( section )
        linkData.append( (mt,mf,mt,mf, idNow) )

        # end loop over NK subsections

    if dat.index != dat.length: raise BadCovariance("Not all covariance data converted, MF%d MT%d" % (mf,mt))
    return sectionList, linkData

def readMF40(info,dat,mf,mt,cov_info,warningList):
    """ production of radioactive isotopes. Also very similar to MF33 """

    dat = myIter(dat)
    sectionList, linkData = [],[]
    ZA,AWR,LIS,dum,NS,dum = funkyFI(dat.next(), logFile = info.logs)
    # each subsection represents different excited state of residual
    for subsection in range(NS):
        try:
            QM,QI,dum,LFS,dum,NL = funkyFI(dat.next(), logFile = info.logs)
        except StopIteration:
            warningList.append('MF40 MT%d lists incorrect number of subsections!' % mt )
            info.doRaise.append( warningList[-1] )
            break
        for subsubsection in range(NL):
            # each subsubsection is a single matrix
            XMF1,XLFS1,MAT1,MT1,NC,NI = funkyFI(dat.next(), logFile = info.logs)
            XMF1,XLFS1 = int(XMF1),int(XLFS1) # XLFS1: level index

            covarsThisSection = []
            if XMF1 not in (0,10):
                raise BadCovariance( "non-zero XMF1/XLFS1 in covariances not currently handled!" )
            if MAT1!=0:
                warningList.append( "cross-material covariance with MAT=%d" % MAT1 )

            for NCdx in range(NC):
                dum,dum,dum,LTY,dum,dum = funkyFI(dat.next(), logFile = info.logs)
                if LTY==0:
                    E1,E2,dum,dum,NCI2,NCI = funkyFI(dat.next(), logFile = info.logs)
                    subsec = []
                    nlines = int(math.ceil(NCI2/6.0))
                    for line in range(nlines): subsec += funkyF(dat.next(), logFile = info.logs)
                    #coefs = subsec[:NCI2][::2]
                    pointerList = [
                            linkModule.link( 'summand', genID( cov_info, int(mtnum), mf )[0],
                            attributes={'ENDF_MFMT':"%d,%d"%(mf,mtnum), 'coefficient':coef})
                            for coef,mtnum in zip(subsec[:NCI2][::2],subsec[:NCI2][1::2])
                            ]
                    covarsThisSection.append( covarianceSummedModule.summedCovariance(
                            label = info.style,
                            lowerBound = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( E1 ),'eV'),
                            upperBound = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( E2 ),'eV'),
                            pointerList=pointerList) )
                else:
                    warningList.append( 'non-zero LTY in MF40' )

            for NIdx in range(NI):
                dum,dum,LS,LB,NT,NP = funkyFI(dat.next(), logFile = info.logs)
                matrix = readMatrix( info, LS,LB,NT,NP, dat )
                if LB not in (0,1,5,6):
                    warningList.append( 'skipping LB%d section for MF%d MT%d' % ( LB, mf, mt ) )
                    continue
                Type='relative'
                if LB==0:
                    Type='absolute'
                    matrix.axes[0].unit = 'b**2'

                covarsThisSection.append( covarianceBaseModule.covarianceMatrix( label = info.style, type=Type, matrix=matrix) )

            if MAT1!=0:
               continue

            # create unique id for each section:
            idNow, pointers = genID( cov_info, mt, mf, MT2=MT1, MF2=(XMF1 or mf), MAT2=MAT1, QI=QI )
            rowdat, coldat = pointers
            section = covarianceSectionModule.section( id=idNow, rowData=rowdat, columnData=coldat )

            if len(covarsThisSection)>1:
                form = covarianceMixedModule.mixedForm( label = info.style, components=covarsThisSection )
                for idx in range(len(form)): form[idx].label = str(idx)
            elif len(covarsThisSection)==1:
                form = covarsThisSection[0]
            else:
                #raise Exception("Encountered empty covariance section!!!")
                info.logs.write("Missing covariance data from section!")
                continue
            section.add( form )

            sectionList.append( section )
            linkData.append( (mt,mf,MT1,XMF1, idNow) )
            # end loop over NL subsections

    if dat.index != dat.length:
        warningList.append( "Not all covariance data converted for MF%d MT%d" % (mf,mt) )
        info.doRaise.append( warningList[-1] )
    return sectionList, linkData

def fillRemainingProductsResidualForBreakup( info, decayChannel, lightIsotopeNames, breakupProducts, residualZA ) :

    residualZA2 = residualZA
    for productName in lightIsotopeNames :
        if( productName in breakupProducts ) :
            multiplicity = breakupProducts[productName]
            product = toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameGamma( info, productNameToZA[productName] ), multiplicity = multiplicity )
            decayChannel.products.add( decayChannel.products.uniqueLabel( product ) )
            if( ( residualZA % 1000 ) > 0 ) :
                residualZA2 -= multiplicity * productNameToZA[productName]
            else :
                residualZA2 -= multiplicity * ( 1000 * ( productNameToZA[productName] / 1000 ) )
    if( residualZA2 != 0 ) : decayChannel.products.add( decayChannel.products.uniqueLabel( 
            toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameGamma( info, residualZA2 ) ) ) )

def parseReaction( info, target, projectileZA, targetZA, MT, MTData, warningList, parseCrossSectionOnly = False ) :
    """
    Translate all available data for the reaction with this MT.
    :return: tuple(crossSection, outputChannel, MFKeys)  where MFKeys contains MF numbers that remain untranslated (should be empty)
    """

    productList = []
    MFKeys = MTData.keys( )
    info.logs.write( '    %3d %s' % ( MT, sorted( MFKeys ) ) )

    for MF in [ 8, 9, 10, 31, 32, 33, 34, 35, 40, 45 ] :
        if( MF in MFKeys ) : MFKeys.remove( MF )

    if( ( MT == 3 ) and ( 3 not in MFKeys ) ) :   # Kludge, for ENDF files that for MT 3 have MF 12 and 13 but not MF 3 data.
        QM, QI, crossSection, breakupProducts = 0, 0, None, None
    else :
        QM, QI, crossSection, breakupProducts = readMF3( info, MT, MTData[3], warningList )
        MFKeys.remove( 3 )
    if( parseCrossSectionOnly ) :
        channel = channelsModule.NBodyOutputChannel( )
        channel.Q.add( toGNDMisc.returnConstantQ( info.style, QM ) )
        return( crossSection, channel, MFKeys )

    fissionGenre = { 18 : channelsModule.fissionGenreTotal,
                     19 : channelsModule.fissionGenreFirstChance,
                     20 : channelsModule.fissionGenreSecondChance,
                     21 : channelsModule.fissionGenreThirdChance,
                     38 : channelsModule.fissionGenreFourthChance }.get(MT, None)

    neutronMFs = []
    for MF in [ 4, 5 ] :
        if( MF in MFKeys ) : neutronMFs.append( MF )
    if( ( neutronMFs != [] ) and ( 6 in MFKeys ) ) : raise Exception( "MF 6 present and MF 4 and/or 5 present: not allowed" ) # This should never happen!

    endfMTProductList = endf_endl.endfMTtoC_ProductLists[MT]
    compoundZA = calculateZA( targetZA, projectileZA, minus = False )
    lightIsotopeZAs = sorted( [ productNameToZA[product] for product in lightIsotopeNames ] )
    lightIsotopeZAsMultiplicity = {}
    for product in lightIsotopeNames : lightIsotopeZAsMultiplicity[productNameToZA[product]] = endfMTProductList.productCounts[product]

    if( ( 4 in neutronMFs ) or ( ( MT == 18 ) and ( neutronMFs == [ 5 ] ) ) ) : # MT == 18 and neutronMFs == [ 5 ] is a special case for bad data (g + Am241).
        ZAP = 1
        if( MT not in [ 2, 5, 18, 19, 20, 21, 38 ] )  :                # Not elastic, fission or sumOfRemainingReactions.
            for product in lightIsotopeNames :
                if( endfMTProductList.productCounts[product] > 0 ) : break
            ZAP = productNameToZA[product]

    levelIndex, decayChannel, twoBodyResidualZA = None, None, None
    if(    50 <= MT < 91 ) :
        levelIndex = MT - 50
    elif( 600 <= MT < 649 ) :
        levelIndex = MT - 600
    elif( 650 <= MT < 699 ) :
        levelIndex = MT - 650
    elif( 700 <= MT < 749 ) :
        levelIndex = MT - 700
    elif( 750 <= MT < 799 ) :
        levelIndex = MT - 750
    elif( 800 <= MT < 849 ) :
        levelIndex = MT - 800
    elif( 875 <= MT < 891 ) :
        levelIndex = MT - 875
    if( MT in [ 91, 649, 699, 749, 799, 849, 891 ] ) : levelIndex = 'c'
    if( MT in [ 103, 104, 105, 106, 107 ] ) : levelIndex = 's'
    level = QM - QI                                                 # If level > 0., residual is in an excited state.
    if( breakupProducts is not None ) :
        if( 50 <= MT < 91 ) :
            level = -QI
        elif( MT == 91 ) :
            pass
        else :
            print breakupProducts
            raise NotImplementedError( 'breakup for MT %s is not supported' % MT )
    isUndefinedTwoBody = ( MT == 91 ) or ( 102 < MT <= 107 ) or ( MT in [ 649, 699, 749, 799, 849 ] )
    isTwoBody = ( MT == 2 ) or ( 50 <= MT < 91 ) or ( ( 600 <= MT < 849 ) and not( isUndefinedTwoBody ) )

    if( isTwoBody or isUndefinedTwoBody ) :
        if( MT == 2 ) :
            ZAP = projectileZA
        else :
            for productName in endfMTProductList.productCounts :
                if( endfMTProductList.productCounts[productName] != 0 ) : break
            if( productName == 'gamma' ) :
                ZAP = 0
            else :
                ZAP = productNameToZA[productName]
        twoBodyResidualZA = calculateZA( compoundZA, ZAP )
    undefinedLevelInfo = { 'ZA' : twoBodyResidualZA, 'level' : level, 'levelIndex' : levelIndex, 'count' : 0 }
    if( neutronMFs == [ 4 ] ) :                     # This is a two-body reaction with only angular data.
        if( not( isTwoBody ) ) : raise ValueError( 'With only MF = 4 data present, reaction is assumed to be two-body and it is not for MT = %s' % MT )
        product = toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameGamma( info, ZAP ) )
        form = readMF4( info, product, MT, MTData[4], angularModule.twoBodyForm, warningList )
        MFKeys.remove( 4 )
        productList.append( product )
    elif( ( neutronMFs == [ 4, 5 ] ) or ( ( neutronMFs == [ 5 ] ) and ZAP == 1 ) ) :
            # Don't check ZAP if MT=5. Currently this combination, MT=5, MF=4/5 appears only for incident gammas
        if( MT != 5 and ZAP != 1 ) : raise ValueError( 'ZAP = %d != 1 for MFs = [ 4, 5 ] for MT = %d' % ( ZAP, MT ) )
        multiplicity = 1
        if( MT not in [ 2, 5, 18, 19, 20, 21, 38 ] )  :                # Not elastic or fission.
            for product in lightIsotopeNames :
                if( endfMTProductList.productCounts[product] > 0 ) : break
            ZAP = productNameToZA[product]
            multiplicity = endfMTProductList.productCounts[product]
        else :
            if( MT not in [ 2, 5 ] ) : multiplicity = -1            # MT 5 is special case where (g,n') put into MT 5 instead of one of 50-91.
        product = toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, ZAP, undefinedLevelInfo ), multiplicity = multiplicity )

        if( neutronMFs == [ 5 ] ) :
            angularSubform = angularModule.isotropic( )             # MF = 5 data is always in lab frame.
        else :
            angularSubform = readMF4( info, product, MT, MTData[4], None, warningList )
            MFKeys.remove( 4 )

        energySubform , weights = readMF5( info, MT, MTData[5], warningList, product = product )
        MFKeys.remove( 5 )

        form = uncorrelated( info.style, frames[1], angularSubform, energySubform )  # BRB: is frame right.
        product.distribution.add( form )
        productList.append( product )
    elif( 6 in MFKeys ) :
        isTwoBody = readMF6( MT, info, MTData[6], productList, warningList, undefinedLevelInfo, isTwoBody )
        MFKeys.remove( 6 )
    elif( neutronMFs == [] ) :
        if( isTwoBody and False ) :                 # ????????? Why False
            raise Exception( 'How did we get here.' )
            product = toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameGamma( info, ZAP ) )
            residualZA = calculateZA( compoundZA, ZAP )
            if( levelIndex is not None ) :
                if( ( levelIndex <= info.targetLevel ) and ( info.target.getZ_A_SuffixAndZA( )[3] == residualZA ) ) : levelIndex -= 1
            if( QI != QM ) :     # Residual is in an excited state.
                decayChannel = channelsModule.NBodyOutputChannel( )
                decayChannel.products.add( decayChannel.products.uniqueLabel( 
                        toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameGamma( info, residualZA ) ) ) )
            residual = toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameGamma( info, residualZA, level = level ), outputChannel = decayChannel )
            productList.append( product )
            productList.append( residual )
    else :
        pass

    readMF12_13( info, MT, MTData, productList, warningList )

    for MF in [ 12, 13, 14, 15 ] :
        if( MF in MFKeys ) : MFKeys.remove( MF )

    """ # FIXME: doesn't appear to be used anymore
    specialBe9n2nExcited, specialBe9n2nExcitedLevel = False, 0
    if( 875 <= MT < 892 ) :                                         # Special case for (z,2n[?])
        if( targetZA == 4009 ) :                                    # Special case for Be9(n,2n[?]He4)He4
            specialBe9n2nExcited = True
            specialBe9n2nExcitedLevel = ( QM - QI ) / 1e6
    """

    if( MT == 5 ) :
        if( QM != QI ) : info.logs.write( '    --QM %s != QI = %s\n' % ( QM, QI ) )
        outputChannel = channelsModule.sumOfRemainingOutputChannels( )
        outputChannel.Q.add( toGNDMisc.returnConstantQ( info.style, QM ) )
    elif( ( MT == 102 ) and not( isTwoBody ) ) :
        residualIndex, gammaMissing = -1, False
        for index, product in enumerate( productList ) :
            if( product.name != 'gamma' ) : residualIndex = index
            gammaMissing = ( product.name == 'gamma' ) or gammaMissing
        if( residualIndex == -1 ) :
            productList.insert( 0, toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, calculateZA( compoundZA, 0 ), undefinedLevelInfo ) ) )
        if( residualIndex > 0 ) : productList.insert( 0, productList.pop( residualIndex ) )
        if( not( gammaMissing ) ) : productList.append( toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, 0, undefinedLevelInfo ) ) )
        outputChannel = channelsModule.NBodyOutputChannel( )
        outputChannel.Q.add( toGNDMisc.returnConstantQ( info.style, QM ) )  # Q????? What about QI?
    elif( isTwoBody ) :
        if( ( QI == 0 ) and ( QM != 0 ) ) : raise Exception("QI = 0, QM = %f for MT=%d" % (QM,MT))
        outputChannel = channelsModule.twoBodyOutputChannel( )
        outputChannel.Q.add( toGNDMisc.returnConstantQ( info.style, QI ) )
        if( len( productList ) == 0 ) :
            for ZA in lightIsotopeZAs :
                if( lightIsotopeZAsMultiplicity[ZA] != 0 ) :
                    productList.append( toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, ZA, undefinedLevelInfo ) ) )
                    break
            if( len( productList ) == 0 ) :
                if( MT != 2 ) : raise Exception( "product data for reaction MT = %s needs to be implemented" % MT )
                productList.append( toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, projectileZA, undefinedLevelInfo ) ) )
        decayProductList = productList[1:]
        productList = productList[:1]                               # Assume first product is "b" in "a + A -> b + B" where B is the larger product.
        ZA = productList[0].particle.getZ_A_SuffixAndZA( )[3]
        residualZA = calculateZA( compoundZA, ZA )
        levelIndex = undefinedLevelInfo['levelIndex']
        if( levelIndex is not None ) :
            if( ( levelIndex <= info.targetLevel ) and ( info.target.getZ_A_SuffixAndZA( )[3] == residualZA ) ) : levelIndex -= 1
        undefinedLevelInfo['levelIndex'] = levelIndex
        for index, product in enumerate( decayProductList ) :
            ZA = product.particle.getZ_A_SuffixAndZA( )[3]
            if( residualZA == ZA ) :
                productList.append( decayProductList.pop( index ) )
                break
        if( len( productList ) < 2 ) :
            if( MT == 2 ) :
                productList.append( toGNDMisc.newGNDParticle( info, target ) )
            else :
                if( ZA == undefinedLevelInfo['ZA'] ) : undefinedLevelInfo['ZA'] = None
                productList.append( toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, residualZA, undefinedLevelInfo ) ) )
            if( info.style in productList[0].distribution ) :
                recoilForm = angularModule.twoBodyForm( info.style, standardsModule.frames.centerOfMassToken,
                        angularSubform = angularModule.recoil( productList[0].distribution[info.style] ) )
                productList[-1].distribution.add( recoilForm )

        decayZAs, decayGammaList, decayNonGammaList = 0, [], []
        for decayProduct in decayProductList :
            decayZAs += decayProduct.particle.getZ_A_SuffixAndZA( )[-1]
            if( decayProduct.name == 'gamma' ) :
                decayGammaList.append( decayProduct )
            else :
                decayNonGammaList.append( decayProduct )
        if( decayZAs == 0 ) :
            if( len( decayGammaList ) != 0 ) :
                if( len( decayNonGammaList ) == 0 ) :
                    decayNonGammaList.append( toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, residualZA, None ) ) )
            elif( len( decayNonGammaList ) != 0 ) :
                if( len( decayGammaList ) == 0 ) : decayGammaList.append( toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, 0, None ) ) )
            decayProductList = decayNonGammaList + decayGammaList
        else :
            raise Exception( "decayZAs = %d != 0" % decayZAs )

        if( breakupProducts is not None ) :
            if( decayChannel is not None ) : raise Exception( 'breakupProducts and decayChannel both not None' )
            decayChannel = channelsModule.NBodyOutputChannel( )
            decayChannel.Q.add( toGNDMisc.returnConstantQ( info.style, QM - QI ) )
            fillRemainingProductsResidualForBreakup( info, decayChannel, lightIsotopeNames, breakupProducts,
                productList[1].particle.getZ_A_SuffixAndZA( )[3] )
            productList[1].addOutputChannel( decayChannel )
        elif( len( decayProductList ) > 0 ) :                         # At this point, both two bodies are in productList and second one is redisual.
            if( QI > QM ) : raise Exception( "Negative decay Q-value for MT%d, QI = %s, QM = %s" % ( MT, QI, QM ) )
            decayChannel = channelsModule.NBodyOutputChannel( )
            decayChannel.Q.add( toGNDMisc.returnConstantQ( info.style, QM - QI ) )  # Q????? Not right?
            for decayProduct in decayProductList : decayChannel.products.add( decayChannel.products.uniqueLabel( decayProduct ) )
            productList[1].addOutputChannel( decayChannel )

    elif( endfMTProductList.isFission ) :
        if hasattr( info, 'fissionEnergies' ):
            ER = info.fissionEnergies[info.style].data[ 'nonNeutrinoEnergy' ] # gets whole polynomial + uncertianty
            useThisQM = ER[ 0 ][ 0 ] # grab constant term and take the value, not the uncertainty
            if abs( QM - useThisQM )/QM > 1e-7: warningList.append( "Fission QM inconsistent with energy release data for MT = " + str( MT ) )
        else:
            useThisQM = QM
        outputChannel = channelsModule.fissionChannel( fissionGenre = fissionGenre )
        outputChannel.Q.add( toGNDMisc.returnConstantQ( info.style, useThisQM ) )
        if( MT == 18 ) :
            if( hasattr( info, 'fissionEnergies' ) ) : outputChannel.addFissionEnergyReleased( info.fissionEnergies )
            if( len( productList ) > 0 ) : outputChannel.products.add( outputChannel.products.uniqueLabel( productList.pop( 0 ) ) )
            if( len( outputChannel ) == 0 ) :
                multiplicity = multiplicityModule.unknown( info.style )
                product = toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameGamma( info, 1 ), multiplicity = multiplicity )
                outputChannel.products.add( outputChannel.products.uniqueLabel( product ) )
            else :
                for product in outputChannel :
                    if( product.name == 'n' ) : break
                info.firstFissionNeutron = product
                if( tokensModule.promptToken in info.totalOrPromptFissionNeutrons ) :
                    product.multiplicity.remove( info.style )
                    product.multiplicity.add( info.totalOrPromptFissionNeutrons[tokensModule.promptToken] )
                    product.multiplicity.remove( 'constant' )
                    product.addAttribute( 'emissionMode', tokensModule.promptToken )
# BRB uncomment                   if( 'total' in info.totalOrPromptFissionNeutrons ) :
#   ditto                     warningList.append( 'have prompt fission nu_bar so not including total' )
                elif( 'total' in info.totalOrPromptFissionNeutrons ) :
                    product.multiplicity.remove( info.style )
                    product.multiplicity.add( info.totalOrPromptFissionNeutrons['total'] )
                    product.multiplicity.remove( 'constant' )
                    product.addAttribute( 'emissionMode', 'total' )
                if( hasattr( info, 'delayedFissionDecayChannel' ) ) :
                    for delayedNeutron in info.delayedFissionDecayChannel : outputChannel.products.add( outputChannel.products.uniqueLabel( delayedNeutron ) )
        else :
            if( neutronMFs == [] ) :
                pass    # we used to add a reference to total fission nubar and PFNS, but that's not physically correct
                """
                if( hasattr( info, 'firstFissionNeutron' ) ) :
                    multiplicity = multiplicityModule.reference( link=info.firstFissionNeutron.multiplicity, label = info.style )
                else :                                              # When singleMTOnly is fission MT != 18.
                    multiplicity = multiplicityModule.unknown( label = info.style )
                product = toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameGamma( info, 1 ), multiplicity = multiplicity )
                if( hasattr( info, 'firstFissionNeutron' ) ) :
                    form = referenceModule.form( link = info.firstFissionNeutron.distribution, label = info.style )
                    product.distribution.add( form )
                outputChannel.products.add( outputChannel.products.uniqueLabel( product ) )
                """

        # July 2011: some files have distributions for 1stChanceFission etc, but should still link to total nubar:
        for product in productList:
            multiplicity = product.multiplicity[info.style]
            if( ( isinstance( multiplicity, multiplicityModule.constant ) ) and ( product.multiplicity[info.style].getValue( 0 ) == -1 ) ) :
                if hasattr( info, 'firstFissionNeutron' ):
                    product.multiplicity.remove( info.style )
                    multiplicity = multiplicityModule.reference( info.firstFissionNeutron.multiplicity, label = info.style )
                    product.multiplicity.add( multiplicity )

        while( len( productList ) > 0 ) : outputChannel.products.add( outputChannel.products.uniqueLabel( productList.pop( 0 ) ) )
    else :
        outputChannel = channelsModule.NBodyOutputChannel( )
        outputChannel.Q.add( toGNDMisc.returnConstantQ( info.style, QI ) )      # Q?????
        if( MT not in [ 18, 19, 20, 21, 38 ] ) :
            residualZA, ZAsMultiplicities, productAsResidual, biggestProduct = compoundZA, {}, None, 0
            for index, product in enumerate( productList ) :
                if( product.name == 'gamma' ) : continue
                Z, A, suffix, ZA = product.particle.getZ_A_SuffixAndZA( )
                multiplicity = product.multiplicity[info.style]
                if( isinstance( multiplicity, multiplicityModule.constant ) ) :
                    m = multiplicity.value
                else :
                    info.logs.write( '\n\nIncorrect multiplicity in ENDF file! MT = %s\n' % MT )
                    info.logs.write( 'Multiplicity should be constant but is (%s).\n' % product.multiplicity.nativeData )
                    info.logs.write( '%s\n' % product.multiplicity[product.multiplicity.nativeData] )
                    raise 'multiplicity should be a constant and it is not.'
                if( ZA in lightIsotopeZAs ) :
                    residualZA = calculateZA( residualZA, m * ZA, minus = True )
                        # If we have different distributions for both neutrons in (n,2n), n shows up twice in the productList.
                    if ZA in ZAsMultiplicities: ZAsMultiplicities[ZA] += m
                    else: ZAsMultiplicities[ZA] = m
                else :
                    if( productAsResidual is not None ) :
                        raise Exception( 'multiple residuals for MT = %, %s %s' % ( MT, productAsResidual.name, product.name ) )
                    productAsResidual = product

            _residualZA = compoundZA
            for ZA in lightIsotopeZAsMultiplicity : _residualZA = calculateZA( _residualZA, lightIsotopeZAsMultiplicity[ZA] * ZA, minus = True )
            if( residualZA != 0 ) :
                for ZA in lightIsotopeZAs :
                    if( ZA not in ZAsMultiplicities ) : ZAsMultiplicities[ZA] = 0
                    if( ZAsMultiplicities[ZA] == lightIsotopeZAsMultiplicity[ZA] ) : continue       # All this ZA accounted for.
                    if( ZAsMultiplicities[ZA] > lightIsotopeZAsMultiplicity[ZA] ) :
                        raise Exception( 'negative multiplicity for ZA = %s for MT = %s' % ( ZA, MT ) )
                    multiplicity = lightIsotopeZAsMultiplicity[ZA] - ZAsMultiplicities[ZA]
                    productList.append( toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, ZA, None ), multiplicity = multiplicity ) )
                    residualZA = calculateZA( residualZA, multiplicity * ZA, minus = True )
                if( productAsResidual is None ) :
                    if( residualZA > 0 ) : productList.append( toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, residualZA, undefinedLevelInfo ) ) )

            if( MT in [ 103, 104, 105, 106, 107, 91, 649, 699, 749, 799, 849, 891 ] ) :
                gammaIndices = []
                for index, product in enumerate( productList ) :
                    if( product.name == 'gamma' ) : gammaIndices.append( index )
                if( len( gammaIndices ) > 0 ) :
                    decayChannel = channelsModule.NBodyOutputChannel( )
                    decayChannel.Q.add( toGNDMisc.returnConstantQ( info.style, QM - QI ) )
                    finalResidual = toGNDMisc.newGNDParticle( info, toGNDMisc.getTypeNameENDF( info, residualZA, None ) )
                    decayChannel.products.add( decayChannel.products.uniqueLabel( finalResidual ) )
                    if( productList[-1].getAttribute( 'ENDFconversionFlag' ) ) : finalResidual.addAttribute( 'ENDFconversionFlag', productList[-1].getAttribute( 'ENDFconversionFlag' ) )
                    for index in gammaIndices : decayChannel.products.add( decayChannel.products.uniqueLabel( productList[index] ) )
                    gammaIndices.reverse( )
                    for index in gammaIndices : del productList[index]
                    productList[-1].addOutputChannel( decayChannel )
                    if( productList[-1].distribution.hasData() ) :  # One must be careful here as this assumes that distributions
                        u_distribution = productList[-1].distribution             # can be simply traded.
                        productList[-1].distribution = finalResidual.distribution
                        finalResidual.distribution = u_distribution
            if( breakupProducts is not None ) :
                if( MT == 91 ) :
                    if( decayChannel is not None ) : raise Exception( 'breakupProducts and decayChannel both not None' )
                    decayChannel = channelsModule.NBodyOutputChannel( )
                    decayChannel.Q.add( toGNDMisc.returnConstantQ( info.style, QM - QI ) )
                    fillRemainingProductsResidualForBreakup( info, decayChannel, lightIsotopeNames, breakupProducts,
                        productList[1].particle.getZ_A_SuffixAndZA( )[3] )
                    productList[1].addOutputChannel( decayChannel )
                else :
                    raise Exception( 'breakup not supported for MT %d' % MT )

    for product in productList :
        if( len( product.distribution ) == 0 ) :
            frame = standardsModule.frames.labToken
            if( isTwoBody ) : frame = standardsModule.frames.centerOfMassToken
            form = unspecifiedModule.form( info.style, productFrame = frame )
            product.distribution.add( form )
        outputChannel.products.add( outputChannel.products.uniqueLabel( product ) )

    return( crossSection, outputChannel, MFKeys )

def parseCovariances( info, MTDatas, MTdict, singleMTOnly=None, resonances=None ):

    covarianceSuite = covarianceSuiteModule.covarianceSuite( )
    linkData = []    # which mf/mts need to be linked for each covariance?
    if( singleMTOnly ) : return( covarianceSuite, linkData )

    # make list of available covariance information:
    warningList = []
    cov_info = {'MTL':{}, 'MTL_2':{}, 'lumpedChannels':{}, 'externalReactions':{}, 'mfmts':[], 'MTdict':MTdict,
            'resonances':resonances, 'NC_data':[], 'style' : info.style }
    for mt in MTDatas.keys():
        if (singleMTOnly is not None) and (mt!=singleMTOnly): continue
        for mf in MTDatas[mt].keys():
            if mf >= 30:    # all covariance-type data
                cov_info['mfmts'].append((mf,mt))
    cov_info['mfmts'].sort()  # sorting first by MF, then by MT

    for mf,mt in cov_info['mfmts']:
        try:
            if mf in (31,33):
                covars, tmp = readMF31_33( info, MTDatas[mt][mf], mf, mt, cov_info, warningList )
            elif mf == 32:
                covars, tmp = readMF32( info, MTDatas[151][32], mf, mt, cov_info, warningList )
            elif mf == 34:
                covars, tmp = readMF34( info, MTDatas[mt][mf], mf, mt, cov_info, warningList )
            elif mf == 35:
                covars, tmp = readMF35( info, MTDatas[mt][mf], mf, mt, cov_info, warningList )
            elif mf == 40:
                covars, tmp = readMF40( info, MTDatas[mt][mf], mf, mt, cov_info, warningList )
            else:
                warningList.append( 'MF%d not yet supported' % mf)
                continue
            for cov in covars:
                if mf == 32: covarianceSuite.addModelParameterCovariance(cov)
                else:
                    covarianceSuite.addSection(cov)
                    if cov.columnData is None:
                        centralValue = cov.rowData.link
                        if( mf == 35 ) : centralValue = centralValue.data
                        centralValue.uncertainties = uncertaintiesModule.uncertainties(
                                [ uncertaintiesModule.uncertainty( type='covariance', functional=linkModule.link(link=cov) ) ] )
            linkData += tmp
        except BadCovariance, e:
            warningList.append('MF%d MT%d covariance conversion failed with message "%s"' % (mf,mt,e) )
            info.doRaise.append( warningList[-1] )

    # fix links for summed matrices:
    for summedMatrix in cov_info['NC_data']:
        for pointer in summedMatrix.pointerList:
            pointed_to = [sec for sec in covarianceSuite.sections if sec.columnData is None
                    and sec.rowData['ENDF_MFMT'] == pointer['ENDF_MFMT']]
            if len(pointed_to) != 1:
                thisMFMT = summedMatrix.ancestor.rowData['ENDF_MFMT']
                warningList.append( "Covariance for MF,MT=%s attempts to sum over non-existant covariance %s"
                        % (thisMFMT, pointer['ENDF_MFMT']) )
                info.doRaise.append( warningList[-1] )
                continue
            pointer.link = pointed_to[0]

    # fix lumped channel covariances (MT851-871) and summed channels (MT1,4,103-107)
    summedReactions = cov_info['MTL'].copy();  summedReactions.update( cov_info['MTL_2'] )
    for (mt,mf) in summedReactions:
        try:
            lumpedChannels = cov_info['lumpedChannels'][(mt,mf)]
        except KeyError as e:
            warningList.append("Cannot find lumped channel %s" %str(e) )
            continue
        for (mt2,mf2) in summedReactions[(mt,mf)]:
            if mt not in range(851,872) and mt2 not in cov_info['MTdict']: continue
            label, pointers = genID(cov_info, mt2, mf2)
            lumpedChannels.reactions.append( pointers[0] )
    for lc in sorted(cov_info['lumpedChannels']):
        covarianceSuite.addReactionSum( cov_info['lumpedChannels'][lc] )
    for exReac in sorted(cov_info['externalReactions']):
        covarianceSuite.addExternalReaction( cov_info['externalReactions'][exReac] )

    if cov_info.get('MF34_missingFrames'):
        for (mt, LegendreLVals) in cov_info['MF34_missingFrames'].items():
            reaction, = [tmp for tmp in info.reactionSuite.reactions if tmp.ENDF_MT==mt]
            # MF=34 is only for neutrons, so only need to look at neutron product:
            frame = reaction.outputChannel.getProductWithName('n').distribution[ info.style ].productFrame
            for lval in LegendreLVals:
                lval.frame = frame
                # FIXME: add an 'endfConversionFlag' to write back as LCT=0

    # add labels
    for label,section in enumerate(covarianceSuite.sections + covarianceSuite.modelParameterCovariances):
        section.label = str(label)

    sys.stdout.flush( )
    for warning in warningList : info.logs.write( "       WARNING: %s\n" % warning, stderrWriting = True )
    return covarianceSuite, linkData
