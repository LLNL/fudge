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

"""
This module contains the method toGND for the class endlZA. This module must only be import by the
toGND method of the endlZA class.

Nuclear particle name: Name depends on whether natural or specific isotope as follows,
    isotope: An isotope's name is composed of SA[_L]
        S is the element's symbol (e.g., 'He' for helium).
        A is the nucleon number (e.g., 'He4' for He-4).
        L is an optional level designation (e.g., 'O16_e3' for the third excition state of O16).
            The level name can be an excitation level designation (e.g., 'e#') or an unknown 
            designation (i.e., 'c') where '#' is an integer (e.g., 'e2' for the second excited state).
    natural: A natural's name is composed of SA[_L] (see above for meanings of S, A and L except A is
            always 'natural'.

gamma: name is simply given as 'gamma'.
"""

import os
import copy

from pqu import PQU

import xData.standards as standardsModule
import xData.axes as axesModule
import xData.XYs as XYsModule
import xData.multiD_XYs as multiD_XYsModule

from fudge import gnd
from fudge.legacy.endl import endl2, endlmisc, endl_C, endl_y, endlIClasses
import endf_endl
from . import toGNDMisc as toGNDMiscModule

from fudge.gnd.reactions import reaction as reactionModule
from fudge.gnd.reactions import partialGammaProduction as partialGammaProductionModule
from fudge.gnd.productData.distributions import angular as angularModule
from fudge.gnd.productData.distributions import energy as energyModule
from fudge.gnd.productData.distributions import Legendre as LegendreModule
from fudge.gnd.productData.distributions import uncorrelated as uncorrelatedModule
from fudge.gnd.productData.distributions import angularEnergy as angularEnergyModule
from fudge.gnd.productData.distributions import reference as referenceModule

ENDL_Accuracy = 1e-3
FUDGE_EPS = 1e-8

def getXYInterpolation( data ) :

    if( data.columns != 2 ) : raise Exception( 'Only 2 column data supported, # of columns = %s for %s' % ( data.columns, str( data ) ) )
    if( data.interpolation == 0 ) : return( standardsModule.interpolation.linlinToken )
    if( data.interpolation == 1 ) : return( standardsModule.interpolation.linlogToken )
    if( data.interpolation == 2 ) : return( standardsModule.interpolation.loglinToken )
    if( data.interpolation == 3 ) : return( standardsModule.interpolation.loglogToken )
    raise Exception( 'Unsupported interpolation = "%s" for %s' % ( data.interpolation, str( data ) ) )

def returnConstantQ( label, Q_MeV ) :

    return( toGNDMiscModule.returnConstantQ( label, Q_MeV, unit = 'MeV' ) )

def uncorrelated( style, frame, angularSubform, energySubform ) :

    _angularSubform = uncorrelatedModule.angularSubform( angularSubform )
    _energySubform = uncorrelatedModule.energySubform( energySubform )
    return( uncorrelatedModule.form( style, frame, _angularSubform, _energySubform ) )

def toGND( self, evaluationLibrary, evaluationVersion, xenslIsotopes = None, verbose = 0 ) :
    """Returns an gnd.reactionSuite.reactionSuite for self where self is an endlZA class."""
#
#   A discrete N-body reaction is one where all the yo's are known and they have energy independent multiplicities.
#
    def getDistributionIs( yo, yosDatas ) :
        """Gets the list of I-values matching the list [ 1, 3, 4, 7, 9 ] for yo."""

        Is = []
        for data in yosDatas[yo] :
            if( data.I in ( 1, 3, 4, 7, 9 ) ) : Is.append( data.I )
        return( Is )

    def getIDistributionData( yo, yosDatas, I, S = None ) :
        """Return the data for the requested yo/I/S or None if not data present."""

        for data in yosDatas[yo] :
            if( S is None ) : S = data.S
            if( ( data.I == I ) and ( data.S == S ) ) : return( data )
        return( None )

    def getAllDistributionIs( yoDatas, IMax = 9, IExtras = [] ) :
        """For a yo get all distribution I-values (i.e., I in [1,9] and in IExtras)."""

        Is = []
        for yoData in yoDatas :
            I = yoData.I
            if( ( ( I in IExtras ) or ( I <= IMax ) ) and ( I not in Is ) ) : Is.append( I )
        Is.sort( )
        return( Is )

    def isYoAndResidaulDataTheSame( yo, yosDatas ) :
        """Test if yo and the residual have the same distribution data (i.e., I = 1, 3, 4, 7 and 9)."""

        yoData = yosDatas[yo]
        resData = yosDatas[yo+10]
        if( ( len( yoData ) == 0 ) or ( len( resData ) == 0 ) ) : return( False )
        yoIs = getDistributionIs( yo, yosDatas )
        resIs = getDistributionIs( yo + 10, yosDatas )
        if( yoIs != resIs ) : return( False )
        for I in yoIs :
            yoIData = getIDistributionData( yo, yosDatas, I )
            resIData = getIDistributionData( yo + 10, yosDatas, I )
            if( yoIData.data != resIData.data ) : return( False )
        return( True )

    def getMultiplicityYos( self, yos, residuals, yosDatas ) :
        """Gets the (multiplicty, yo) pairs and adjust residuals if needed."""

        if( len( residuals ) > 1 ) : raise Exception( 'len( residuals ) > 1 not supported: residual = %s' % `residuals` )
        yoOld = None
        mYos = []
        for yo in yos :
            if( yo == yoOld ) :
                n += 1
            else :
                if( yoOld is not None ) : mYos.append( [ n, yoOld ] )
                yoOld = yo
                n = 1
        if( len( yos ) > 0 ) :
            if( residuals[0] == yo ) :
                rYo = endl2.ZAToYo( residuals[0] )
                if( isYoAndResidaulDataTheSame( rYo, yosDatas ) ) :
                    residuals = residuals[1:]           # residuals should now be empty.
                    n += 1
                    yosDatas[rYo + 10] = []
                    for yosData in yosDatas[rYo] :
                        if( yosData.I == 10 ) : yosData.set( yosData * ( n / float( n - 1 ) ) )

            mYos.append( [ n, yo ] )
        if( len( residuals ) > 0 ) : mYos.append( [ 1, residuals[0] ] )
        return( mYos )

    def addDistributionDataAndRemove( particle, yo, yosDatas, promptNeutronParticle = None ) :

        def angular( data ) :

            axes = angularModule.XYs2d.defaultAxes( energyUnit = 'MeV' )
            subform = angularModule.XYs2d( axes = axes )
            for i1, EMuP in enumerate( data.data ) :
                E1, muP = EMuP
                subform[i1] = angularModule.XYs1d( data = muP, accuracy = ENDL_Accuracy, value = E1 )
            return( subform )

        def energy( data ) :

            axes = energyModule.XYs2d.defaultAxes( energyUnit = 'MeV' )
            subform = energyModule.XYs2d( axes = axes, interpolationQualifier = standardsModule.interpolation.unitBaseToken )
            for E1, EpP in data.data[0][1] : subform.append( energyModule.XYs1d( data = EpP, accuracy = ENDL_Accuracy, value = E1 ) )
            return( subform )

        def LLNLLegendrePointwise( data ) :

            axes = LegendreModule.LLNLPointwise.defaultAxes( )
            subform = LegendreModule.LLNLPointwise( axes )
            axes2D = axesModule.axes( 3 )
            axes2D[0] = axes[0]
            axes2D[1] = axes[1]
            axes2D[2] = axes[2]
            axes1D = axesModule.axes( 2 )
            axes1D[0] = axes[0]
            axes1D[1] = axes[1]
            for l, EEpPs in data.data :
                w_xys = multiD_XYsModule.XYs2d( axes = axes2D, value = l )
                for index, E_EpPs in enumerate( EEpPs ) :
                    xPrior = -1
                    E1, EpPs = E_EpPs
                    for xy in EpPs :
                        if( xy[0] == xPrior ) : xy[0] *= ( 1 + FUDGE_EPS )
                        xPrior = xy[0]
                    w_xys.append( XYsModule.XYs1d( data = EpPs, axes = axes1D, accuracy = ENDL_Accuracy, value = E1 ) )
                subform.append( w_xys )
            return( subform )

        def LLNLLegendrePointwise_L0_only( info, data ) :
            """Has only L=0 form, can be converted to uncorrelated angular (isotropic) and energy distributions."""

            angularSubform = angularModule.isotropic( )
            axes = energyModule.XYs2d.defaultAxes( energyUnit = 'MeV' )
            energySubform = energyModule.XYs2d( axes = axes, interpolationQualifier = standardsModule.interpolation.unitBaseToken )
            for index, ( energy, probability_list ) in enumerate( data.data[0][1] ) :
                energySubform.append( energyModule.XYs1d( data = probability_list, accuracy = ENDL_Accuracy, value = energy ) )
            return( uncorrelated( info.style, standardsModule.frames.labToken, angularSubform, energySubform ) )

        def LLNLAngularEnergy( data ) :

            axes = angularEnergyModule.XYs3d.defaultAxes( energyUnit = "MeV", energy_outUnit = 'MeV', probabilityUnit = '1/MeV',
                    probabilityLabel = 'P(energy_out|energy_in,mu)' )
            subform = angularEnergyModule.XYs3d( axes = axes, interpolationQualifier = standardsModule.interpolation.unitBaseToken  )
            for i, E_MuEpPs in enumerate( data.data ) :
                w_xys = angularEnergyModule.XYs2d( value = E_MuEpPs[0] )
                for j, Mu_EpPs in enumerate( E_MuEpPs[1] ) :
                    w_xys.append( angularEnergyModule.XYs1d( data = Mu_EpPs[1], accuracy = ENDL_Accuracy, value = Mu_EpPs[0] ) )
                subform.append( w_xys )
            return( subform )

        def getSubform( I, datas, func ) :

            for i1, data in enumerate( datas ) :
                if( data.I == I ) :
                    del datas[i1]
                    return( func( data ) )
            raise Exception( 'I = %d not in datas' % I )

        if( yo not in yosDatas ) : return
        datas = yosDatas[yo]
        Is = getAllDistributionIs( datas, IMax = 4, IExtras = [ 941, 942 ] )
        form = None
        if( Is == [] ) :
            if( len( datas ) ) :
                if( ( len( datas ) == 1 ) and ( datas[0].S == 7 ) and ( datas[0].I == 7 ) ) :   # Special case for endl I, S = 7, 7 with no distribution data.
                    form = referenceModule.form( link = promptNeutronParticle.distribution[info.style], label = info.style )
                else :
                    print
                    for data in datas : print data
                    raise Exception( "Unsupported data for particle %s: I's = %s." % ( particle.name, `Is` ) )
        elif( ( 941 in Is ) or ( 942 in Is ) ) :
            raise Exception( 'I = 941 or 942 not supported' )
        elif( Is == [ 1 ] ) : 
            if( datas[0].S == 3 ) :
                energySubform = energyModule.constant( PQU.PQU( PQU.pqu_float.surmiseSignificantDigits( datas[0].getX1( ) ), 'MeV' ) )
                angularSubform = getSubform( 1, datas, angular )
                form = uncorrelated( info.style, standardsModule.frames.labToken, angularSubform, energySubform )
            else :
                for data in datas :
                    if( data.I == 1 ) : productFrame = standardsModule.frames.centerOfMassToken
                subform = getSubform( 1, datas, angular )
                form = angularModule.twoBodyForm( info.style, productFrame, subform )
        elif( Is == [ 4 ] ) :
            for i1 in range( len( datas ) ) :
                if( datas[i1].I == 4 ) : break
            if( ( len( datas[i1].data ) == 1 ) and ( datas[i1].data[0][0] == 0 ) ) :
                form = LLNLLegendrePointwise_L0_only( info, datas[i1] )
                del datas[i1]
            else:
                subform = getSubform( 4, datas, LLNLLegendrePointwise )
                form = LegendreModule.form( info.style, standardsModule.frames.labToken, subform )
        elif( Is == [ 1, 3 ] ) :
            angularSubform = getSubform( 1, datas, angular )
            angularSubform = angularEnergyModule.LLNLAngularOfAngularEnergySubform( angularSubform )
            angularEnergySubform = getSubform( 3, datas, LLNLAngularEnergy )
            angularEnergySubform = angularEnergyModule.LLNLAngularEnergyOfAngularEnergySubform( angularEnergySubform )
            form = angularEnergyModule.LLNLAngularEnergyForm( info.style, standardsModule.frames.labToken,
                    angularSubform, angularEnergySubform )
        elif( Is == [ 1, 4 ] ) :
            angularSubform = getSubform( 1, datas, angular )
            energySubform = getSubform( 4, datas, energy )
            form = uncorrelated( info.style, standardsModule.frames.labToken, 
                    angularSubform, energySubform )
        else :
            raise Exception( "I's = %s not supported for particle %s." % ( `Is`, particle.name ) )
        for idx in xrange( len( datas ) - 1, -1, -1 ) :
            data = datas[idx]
            if( data.I in [ 7, 9, 10, 13, 941, 942 ] ) :
                if( data.I == 7 ) :
                    if( data.S == 7 ) : 
                        particle.addAttribute( 'emissionMode', gnd.tokens.delayedToken )
                        particle.addAttribute( 'decayRate', PQU.PQU( PQU.pqu_float.surmiseSignificantDigits( data.getX1( ) ), "1/s" ) )
                    else :
                        particle.addAttribute( 'emissionMode', gnd.tokens.promptToken )
                interpolation = getXYInterpolation( data )
                if( data.I in [ 7, 9 ] ) : 
                    axes = gnd.productData.multiplicity.XYs1d.defaultAxes( energyUnit = 'MeV')
                    particle.multiplicity.remove( info.style )
                    particle.multiplicity.add( gnd.productData.multiplicity.XYs1d( data = data.data, label = info.style, 
                            axes = axes, accuracy = ENDL_Accuracy, interpolation = interpolation ) )
                else :          # for I = 10 and 13
                    if( data.I == 10 ) :
                        axes = gnd.productData.energyDeposition.XYs1d.defaultAxes( energyUnit = 'MeV' )
                        dataForm = gnd.productData.energyDeposition.XYs1d( label = info.style, axes = axes, 
                                data = data.data, accuracy = ENDL_Accuracy, interpolation = interpolation )
                        particle.energyDeposition.add( dataForm )
                    elif( data.I == 13 ) :
                        axes = gnd.productData.momentumDeposition.XYs1d.defaultAxes( energyUnit = 'MeV', 
                                momentumDepositionUnit = 'MeV/c' )
                        dataForm = gnd.productData.momentumDeposition.XYs1d( label = info.style, axes = axes, 
                                data = data.data, accuracy = ENDL_Accuracy, interpolation = interpolation )
                        particle.momentumDeposition.add( dataForm )
                del datas[idx]
        if( not form is None ) : particle.distribution.add( form )

    def makeNBodyChannelFrom_mYos( channelClass, mYos, ZAsYos, yosDatas, Q_MeV, specialCase = None, decayChannelInfo = [], levelIndex = None ) :
        """This routine does not support yo's (including the residual) in an excited state. This is consistent with endl99.
        This routine does not work if decayChannelInfo has data (see note Note_Be_9) due to added gndPath logic (gndPath has been
        removed and may be so should this last sentence)."""

        if( ( len( yosDatas[7] ) > 0 ) and ( hasattr( yosDatas[7][0], 'specialCase' ) ) ) :
            specialCase = yosDatas[7][0].specialCase
            if( specialCase == 'Am_242_m1' ) : level = yosDatas[7][0].level
            if( specialCase == 'Sc_45_m1' ) : level = yosDatas[7][0].level
        s = ' '
        channel = channelClass( )
        channel.Q.add( returnConstantQ( info.style, Q_MeV ) )
        i = 0
        n = len( mYos ) - 1
        for m, ZA in mYos :
            decayChannel = None
            if( len( decayChannelInfo ) > 0 ) :
                if( decayChannelInfo[0] == i ) : decayChannel = decayChannelInfo[1]
            if( m > 1 ) :
                particle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, ZA ), multiplicity = m, outputChannel = decayChannel )
            elif( ( specialCase == 'unknown level' ) and ( i == 1 ) ) :
                particle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, ZA, level = 0.0, levelIndex = levelIndex ), outputChannel = decayChannel )
            elif( ( specialCase == 'Am_242_m1' ) and ( ZA == 95242 ) ) :
                particle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, ZA, level = level, levelIndex = 2 ), outputChannel = decayChannel )
            elif( ( specialCase == 'Sc_45_m1' ) and ( ZA == 21046 ) ) :
                particle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, ZA, level = level, levelIndex = 2 ), outputChannel = decayChannel )
            else :
                particle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, ZA ), outputChannel = decayChannel )
            s = ' + '
            if( ZA in ZAsYos ) :
                if( ( specialCase == 'S8' ) and ( ZA == 1 ) ) :
                    yo = 11
                else :
                    yo = ZAsYos[ZA]
                    if( ( i == n ) and ( yosDatas[yo] == [] ) ) : yo += 10
                addDistributionDataAndRemove( particle, yo, yosDatas )
            channel.products.add( channel.products.uniqueLabel( particle ) )
            i += 1
        return( channel )

    isThermalNeutronScatteringLaw = self.ZA in [ 1801, 1901, 1902, 4809, 4909, 6912, 8916 ]
    targetZA = self.ZA
    targetName = None
    if( isThermalNeutronScatteringLaw ) :
        targetZA -= 800
        if( self.ZA in [ 1901, 1902, 4909, 6912, 8916 ] ) : targetZA -= 100
        targetName = { 1801 : "H1_inH2O_TNSL",       1901 : "H1_inCH2_TNSL",    1902 : "H2_inD2O_TNSL",
                       4809 : "Be9_inBeMetal_TNSL",  4909 : "Be9_inBeO_TNSL",   6912 : "C12_inGraphite_TNSL", 8916 : "O16_inBeO_TNSL" }[self.ZA]

        I0 = self.findData( C = 10, I = 0 )
        for energy, xSec in I0.data :
            if( xSec > 0 ) : break
            TNSL_EMax = energy
    if( targetZA in [ 95241, 21045 ] ) :
        I0s = self.findDatas( C = 46, S = 0, I = 0 )
        if( len( I0s ) == 2 ) :
            I0m1, I0g = I0s
            if( I0m1.getQ( ) > I0g.getQ( ) ) : I0m1, I0g = I0g, I0m1
            level = I0g.getQ( ) - I0m1.getQ( )
            dataM1s = self.findDatas( C = 46, S = 0, Q = I0m1.getQ( ) )
            for dataM1 in dataM1s :
                if( targetZA == 95241 ) :
                    dataM1.specialCase = 'Am_242_m1'
                else :
                    dataM1.specialCase = 'Sc_45_m1'
                dataM1.level = level

    styleName = 'eval'
    yosDatas = {  0 : [],  1 : [],  2 : [],  3 : [],  4 : [],  5 : [],  6 : [],  7 : [],  8 : [],  9 : [],
                     10 : [], 11 : [], 12 : [], 13 : [], 14 : [], 15 : [], 16 : [], 17 : [], 18 : [], 19 : [] }
    yosZAs = {}
    for yo in xrange( 1, 17 ) : yosZAs[yo] = endl2.yoToZA( yo )
    ZAsYos = {}
    for yo in xrange( 1, 7 ) : ZAsYos[endl2.yoToZA( yo )] = yo
    ZAsYos[7] = 7
    ENDLCS_To_ENDFMT = endf_endl.ENDLCS_To_ENDFMT( self.projectileName )

    level = None
    I0 = self.findDatas( I = 0 )[0]
    if( I0.getELevel( ) > 0 ) : level = I0.getELevel( )
    info = toGNDMiscModule.infos( styleName, xenslIsotopes, transportables = [ 'n', 'H1', 'H2', 'H3', 'He3', 'He4', 'gamma' ] )
    projectile = toGNDMiscModule.getTypeName( info, endl2.yoToZA( self.yi ) )

    CsAndSs, residualExcitationIndexLevels = {}, {}
    if( level is not None ) : residualExcitationIndexLevels[targetZA] = [ [ None, level ] ]
    reactionsDatas = self.findReactionsDatas( )
    if( 1 in self.CList( ) ) : reactionsDatas += [ self.findDatas( C = 1 ) ]
    compoundZA = endl2.yoToZA( self.yi ) + targetZA
    for reactionDatas in reactionsDatas :               # Determine a unique set of levels for excited residuals.
        C, S = reactionDatas[0].C, reactionDatas[0].S
        if( C not in CsAndSs ) : CsAndSs[C] = []
        if( S not in CsAndSs[C] ) : CsAndSs[C].append( S )
        if( S == 1 ) :
            yos = endl_C.endl_C_yoInfo( C )
            residualZA = compoundZA - yos[0]
            if( residualZA not in residualExcitationIndexLevels ) : residualExcitationIndexLevels[residualZA] = []
            residualExcitationIndexLevels[residualZA].append( [ None, reactionDatas[0].X1 ] )
    for residualZA in residualExcitationIndexLevels :
        residualExcitationIndexLevels[residualZA].sort( )
        indexLevels, indexOffset = {}, 1
        if( residualExcitationIndexLevels[residualZA][0][1] == 0 ) : indexOffset = 0       # Should only happend for meta-stable targets.
        for index, indexAndLevel in enumerate( residualExcitationIndexLevels[residualZA] ) : 
            indexLevel = index + indexOffset
            indexLevels[indexAndLevel[1]] = indexLevel
        residualExcitationIndexLevels[residualZA] = indexLevels

    levelIndex = None
    if( level is not None ) : levelIndex = residualExcitationIndexLevels[targetZA][level]
    target = toGNDMiscModule.getTypeName( info, targetZA, name = targetName, level = level, levelIndex = levelIndex )

    maxDate = 19590101
    allData = self.findDatas( )
    for oneDatum in allData :
        if( 0 <= oneDatum.yo <= 16 ) :
            if( oneDatum.I in [ 80, 81, 84, 89, 90, 91, 92 ] ) : continue
            if( oneDatum.I not in [ 10, 11, 13, 20 ] ) :
                date = oneDatum.getDate( ) + 19 * 1000000
                if( date < 600000 ) : date += 1000000   # Must treat dates after 2000 different then those before (600000 is arbitrary).
                maxDate = max( date, maxDate )
    maxDate = str( maxDate )
    maxDate = '%s-%s-%s' % ( maxDate[:4], maxDate[4:6], maxDate[6:] )
    evaluatedStyle = gnd.styles.evaluated( info.style, 
            PQU.PQU( PQU.pqu_float.surmiseSignificantDigits( I0.getTemperature( ) ), 'MeV/k' ),
            evaluationLibrary, evaluationVersion, date = maxDate )

    reactionSuite = gnd.reactionSuite.reactionSuite( projectile, target, 
            style = evaluatedStyle, particleList = info.particleList )
    info.setReactionSuite( reactionSuite )

    if( self.filelist[0].levels[0].ELevel > 0 ) :
        reactionSuite.addNuclearMetaStableAlias( target.name.split( '_' )[0], target.name, 1 )

    I20s = []                                           # Some I20 data has bad Q values, so let's ignore all I20 data.
    for i, reactionDatas in enumerate( reactionsDatas ) :
        for j, reactionData in enumerate( reactionDatas ) :
            if( reactionData.I == 20 ) : I20s.append( [ i, j ] )
    I20s.reverse( )
    for i, j in I20s :
        del reactionsDatas[i][j]
        if( len( reactionsDatas[i] ) == 0 ) : del reactionsDatas[i]

    channelList, delayedNeutrons = [], {}
    for iChannel, reactionDatas in enumerate( reactionsDatas ) :
        C, S = reactionDatas[0].C, reactionDatas[0].S   # Phase 1
        if( isThermalNeutronScatteringLaw ) :
            if( ( C, S ) != ( 11, 1 ) ) : continue
            I0 = None
            for reactionData in reactionDatas :
                if( reactionData.I == 0 ) : I0 = reactionData
            if( 0.9e-12 < I0.getX1( ) > 2.1e-12 ) : continue
        c_Or_s_Level = 'c'
        if( ( S == 0 ) and ( 1 not in CsAndSs[C] ) ) : c_Or_s_Level = 's'
        maxReactionDate = 0
        I0, I12, I20 = None, None, None
        specialCase = None
        for yo in yosDatas : yosDatas[yo] = []
        for data in reactionDatas :                         # yosDatas should only contain I = 1, 3, 4, 7, 9, 10 and 13 type data.
            data.frame = ( endlmisc.getNumberOfColumns_( data.I, '' ) * ( '%s,' % standardsModule.frames.labToken ) )[:-1]
            if( data.S == 7 ) :                             # Special treatment for delayed neutron data.
                decayRate = data.getX1( )
                if( decayRate not in delayedNeutrons ) : delayedNeutrons[decayRate] = []
                delayedNeutrons[decayRate].append( data )
                continue
            if( 0 <= data.yo <= 16 ) :
                if( data.I in [ 80, 81, 84, 89, 90, 91, 92 ] ) : continue
                if( data.I not in [ 10, 11, 13, 20 ] ) :
                    date = data.getDate( ) + 19 * 1000000
                    if( date < 600000 ) : date += 1000000   # Must treat dates after 2000 different then those before (600000 is arbitrary).
                    maxReactionDate = max( date, maxReactionDate )
                    maxDate = max( date, maxDate )
                if( data.I == 0 ) :
                    I0 = data
                elif( data.I == 11 ) :
                    continue
                elif( data.I == 12 ) :
                    I12 = data
                elif( data.I == 20 ) :
                    I20 = data
                elif( ( data.yo == 0 ) and ( data.I == 13 ) ) :
                    continue
                else :
                    if( data.I not in [ 1, 3, 4, 7, 9, 10, 13, 20, 941, 942 ] ) :
                        print data
                        raise Exception( 'I = %d data not currently supported' % data.I )
                    if( ( data.S == 8 ) and ( data.yo == 1 ) and ( data.I == 4 ) and ( data.C == 12 ) ) :   # This should only be for endl n+d->2n+p
                        if( data.ZA != 1002 ) :
                            print data
                            raise Exception( 'S = 8, yo = 1, I = 4 only supported for ZA = 1002' )
                        yosDatas[data.yo + 10].append( data )
                        specialCase = 'S8'
                    else :
                        yosDatas[data.yo].append( data )
            else :
                print data
                raise Exception( 'Unsupported yo = %d' % data.yo )
        if( specialCase == 'S8' ) :                         # May need to add I = 10 data to yo = 11 and readjust yo = 1's I = 10 data.
            if( len( yosDatas[11] ) ) :
                I10Present = False
                for data in yosDatas[11] :
                    if( data.I == 10 ) : I10Present = True
                if( not I10Present ) :
                    yo1I10 = None
                    for data in yosDatas[1] :
                        if( data.I == 10 ) : yo1I10 = data
                    yo12I10 = None
                    for data in yosDatas[12] :
                        if( data.I == 10 ) : yo12I10 = data
                    if( ( yo1I10 is not None ) and ( yo12I10 is not None ) ) :
                        yosDatas[11].append( yo12I10 )
                        yo1I10.set( yo1I10 - yo12I10 )
        for yo in yosDatas :
            if( ( len( yosDatas[yo] ) == 1 ) and ( yosDatas[yo][0].I == 10 ) ) : yosDatas[yo] = []
        if( I0 is None ) : raise Exception( 'Missing cross section data for reaction C = %d S = %d.' % ( C, S ) )   # End of Phase 1

        if( C == 1 ) :
            yos = ( 0, )
        elif( isThermalNeutronScatteringLaw ) :
            residualZA, yos = targetZA, ( 1, targetZA )

            points = I0.data[:]
            while( points[-1][1] == 0 ) : del points[-1]
            priorEnergy = -1
            for energyXSec in points :
                if( energyXSec[0] == priorEnergy ) : energyXSec[0] += 1e-6 * energyXSec[0]
                priorEnergy = energyXSec[0]
            I0 = endlIClasses.endlI0( None, I0.yo, I0.C, I0.I, I0.S, I0.h, points, I0.bdflsFile )
            if( TNSL_EMax != I0.xMax( ) ) : raise Exception( 'TNSL_EMax = %s != I0.xMax( ) = %s: %s' % ( TNSL_EMax, I0.xMax( ), `I0` ) )
        elif( ( C == 12 ) and ( targetZA == 4009 ) and ( self.yi == 1 ) ) :  # Special case for Be_9(n,2n) which should really be Be_9(n,2n a)
            specialCase = 'Be_9(n, 2n)'
            yos = ( 1, 1, 2004 )
        elif( ( C == 13 ) and ( targetZA == 4010 ) and ( self.yi == 1 ) ) :  # Special case for Be10(n,3n)
            specialCase = 'Be_10(n, 3n)'
            yos = ( 1, 1, 1, 2004 )
        elif( ( C == 14 ) and ( targetZA == 4011 ) and ( self.yi == 1 ) ) :  # Special case for Be_11(n,4n)
            specialCase = 'Be_11(n, 4n)'
            yos = ( 1, 1, 1, 1, 2004 )
        else :
            yos = endl_C.endl_C_yoInfo( C )                 # Phase 2. Determine products (i.e., yos and residuals).
        for yo in yos :
            if( yo in [ 8, 9, 10 ] ) :
                if( C in [ 73, 74 ] ) : continue
                raise Exception( 'Electron, positron or electron capture (i.e., yo = 8, 9 or 10) not supported: C = %d' % C )
        residuals = []
        if( len( yos ) == 0 ) :
            if( C != 1 ) : print 'C = %d with no yos is currently not handled by toGND' % C
            continue
        elif( yos[0] == -1 ) :                              # Unknown yos.
            if( C != 1 ) : print 'C = %d has unknown products which is currently not handled by toGND' % C
            continue
        elif( yos[0] == -2 ) :                              # yo = yi
            if( len( yos ) != 1 ) :
                print 'product list = ', yos
                raise Exception( 'For C = %d, endl_C.endl_C_yoInfo returned -2 plus other products, only the -2 is supported.' % C )
            yos = ( yosZAs[ self.yi ], )
            residuals = [ targetZA ]
        elif( yos[0] == -3 ) :                                          # fission
            residuals = [ 99120, 99120 ]
        elif( yos[0] == -4 ) :                                          # nX?
            pass                                                        # Leave residuals empty as residual is undefined.
        else :
            yiZA = endl2.yoToZA( self.yi )
            if( self.A == 0 ) : yiZA = 1000 * ( yiZA / 1000 )
            productsZA = -yiZA
            for yo in yos :
                if( yo != 7 ) : productsZA += yo
            if( self.A == 0 ) : productsZA = 1000 * ( productsZA / 1000 )
            residuals = [ targetZA - productsZA ]        # End of phase 2

        #FIXME do we still need the next 4 lines? 4008 isn't in residuals for Be9(n,2n), Be10(n,3n) or Be11(n,4n) due to specialCase logic above
        if( 4008 in residuals ) :                                       # Special case for Be_8 breakup into two alphas
            if( ( 2004 not in yos ) and ( yosDatas[6] != [] ) and ( yosDatas[16] == [] ) ) :
                yosDatas[16] = yosDatas[6]
                yosDatas[6] = []

        if( 99000 < targetZA < 99200 ) :
            residuals = [ targetZA ]
        elif( ( C == 30 ) and ( endl2.yoToZA( self.yi ) in [ 1002, 1003 ] ) and ( targetZA in [ 1002, 1003 ] ) ) :   # Special treatment for C = 30 as 
            if( yosDatas[6] != [] ) : raise Exception( 'C = 30 as yo = 6 data' )                    # yos should be [ g, n ] and not [ g, n, a ].
            yos = ( 7, 1 )
            residuals = ( 2004, )

        if( verbose > 0 ) : print '   ', C, S,
        if( isThermalNeutronScatteringLaw ) :
            if( reactionDatas[0].getX1( ) < 1.1e-12 ) :
                TNSL_MT = 4
                channelProcess = "Incoherent inelastic scattering"
            else :
                TNSL_MT = 2
                channelProcess = "Coherent elastic scattering"
                if( self.ZA == 1901 ) : channelProcess = "Inc" + channelProcess[1:]
            outputChannel = gnd.channels.NBodyOutputChannel( process = channelProcess )
            outputChannel.Q.add( returnConstantQ( info.style, I0.getQ( ) ) )

            firstParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, 1 ) )
            addDistributionDataAndRemove( firstParticle, 1, yosDatas )
            outputChannel.products.add( outputChannel.products.uniqueLabel( firstParticle ) )

            for yo in yosDatas : yosDatas[yo] = []
        elif( ( yos[0] > 0 ) and ( C not in [ 71, 72, 73, 74 ] ) ) :  # Reactions that are two-body, two-body + break-up or discrete N-body except C = 71 - 74.
            outputChannel = None
            Is = getAllDistributionIs( yosDatas[ZAsYos[yos[0]]] )
            if( Is == [ 1 ] ) :     # Two-body initial product state.
                if( verbose > 0 ) : print '2Body',
                firstParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, yos[0] ) )    # Assume first particle is in the ground state.
                addDistributionDataAndRemove( firstParticle, ZAsYos[yos[0]], yosDatas )
                Q = returnConstantQ( info.style, I0.getQ( ) )
                if( len( yos ) == 1 ) :                                 # Simple two-body, yos[0] + Residual [ + gamma].
                    resYo = residuals[0]
                    if( resYo in ZAsYos ) : resYo = ZAsYos[resYo] + 10
                    if( S == 1 ) :
                        Q = returnConstantQ( info.style, I0.getQ( ) - I0.getX1( ) )
                        resParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, residuals[0] ) )
                        addDistributionDataAndRemove( resParticle, resYo, yosDatas )

                        decayChannel = None
                        if( 7 in yos ) : print '7 in yos'

                        if( yosDatas[7] != [] ) :
                            decayChannel = gnd.channels.NBodyOutputChannel( )
                            decayChannel.Q.add( returnConstantQ( info.style, I0.getX1( ) ) )  # Assume residual returns to ground.
                            decayChannel.products.add( decayChannel.products.uniqueLabel( resParticle ) )  # state, hence Q = getX1( )
                            gammaParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, 7 ) )
                            Is_ = sorted( [ data.I for data in yosDatas[7] ] )
                            if( Is_ == [ 10, 13 ] ) :
                                print '  gamma has I = 10 and 13 but no others',
                                yosDatas[7] = []
                            else :
                                addDistributionDataAndRemove( gammaParticle, 7, yosDatas )
                            decayChannel.products.add( decayChannel.products.uniqueLabel( gammaParticle ) )

                        level = I0.getX1( )
                        levelIndex = residualExcitationIndexLevels[residuals[0]][level]
                        if( levelIndex is None ) : level = None
                        secondParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, residuals[0], level = level, levelIndex = levelIndex ), 
                            outputChannel = decayChannel )
                    else :
                        decayChannel, level, levelIndex = None, None, None
                        if( S == 0 ) :
                            level = I0.getELevel( )
                            if( C not in [ 8, 9, 10 ] ) :
                                levelIndex = c_Or_s_Level
                                if( ( yosDatas[7] != [] ) and ( residuals[0] != 7 ) ) :
                                    decayChannel = gnd.channels.NBodyOutputChannel( )
                                    decayChannel.Q.add( returnConstantQ( info.style, I0.getX1( ) ) )
                                    resParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, residuals[0] ) )       # Residual in ground state.
                                    addDistributionDataAndRemove( resParticle, resYo, yosDatas )
                                    decayChannel.products.add( decayChannel.products.uniqueLabel( resParticle ) )
                                    gammaParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, 7 ) )
                                    addDistributionDataAndRemove( gammaParticle, 7, yosDatas )
                                    decayChannel.products.add( decayChannel.products.uniqueLabel( gammaParticle ) )
                            elif( I0.getELevel( ) != 0 ) :            # Elastic scattering off an isomer
                                levelIndex = residualExcitationIndexLevels[residuals[0]][level]
                            if( ( level == 0.0 ) and ( levelIndex is None ) ) : level = None
                        if( residuals[0] <= 2004 ) : level, levelIndex = None, None
                        secondParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, residuals[0], level = level, levelIndex = levelIndex ), 
                            outputChannel = decayChannel )
                        addDistributionDataAndRemove( secondParticle, resYo, yosDatas )
                else :                                          # Two-body with residual breaking up.
                    if( verbose > 0 ) : print 'w/breakup',
                    if( self.A == 0 ) : raise Exception( 'For C = %d, breakup of product not supported for two-body reaction for natural target.' )
                    mYos = getMultiplicityYos( self, yos[1:], residuals, yosDatas )
                    residualZA = targetZA + endl2.yoToZA( self.yi ) - yos[0]
                    levelIndex, level = c_Or_s_Level, 0.0
                    if( I0.S == 1 ) :
                        level = I0.getX1( )

                        levelIndex = residualExcitationIndexLevels[residualZA][level]
                    elif( I0.S == 8 ) : 
                        levelIndex, level = 1, I0.getX1( )
                    Q_MeV = endl2.reactionQByZAs( [ endl2.yoToZA( self.yi ), targetZA ], [ yos[0], residualZA ] ) - level
                    residualName = toGNDMiscModule.getTypeName( info, residualZA, level = level, levelIndex = levelIndex )
                    decayChannel = makeNBodyChannelFrom_mYos( gnd.channels.NBodyOutputChannel, mYos, ZAsYos, yosDatas, I0.getQ( ) - Q_MeV, 
                        specialCase = specialCase )
                    secondParticle = toGNDMiscModule.newGNDParticle( info, residualName, outputChannel = decayChannel )
                    Q = returnConstantQ( info.style, Q_MeV )
                outputChannel = gnd.channels.twoBodyOutputChannel( )
                outputChannel.Q.add( Q )
                outputChannel.products.add( outputChannel.products.uniqueLabel( firstParticle ) )
                outputChannel.products.add( outputChannel.products.uniqueLabel( secondParticle ) )
            else :                                              # Mutli-particle breakup (i.e., not two body).
                if( verbose > 0 ) : print 'NBody',
                mYos = getMultiplicityYos( self, yos, residuals, yosDatas )
                gammaPresent = False
                for m, yo in mYos : gammaPresent = gammaPresent or ( yo == 7 )
                if( ( not gammaPresent ) and ( yosDatas[7] != [] ) ) : mYos.append( [ 1, 7 ] )
                decayChannelInfo = []
                if( specialCase == 'Be_9(n, 2n)' ) :                                            # Special case for endl Be_9(n, 2n)
                    if( yosDatas[6] == [] ) : mYos = [[2, 1], [2, 2004] ]
                elif( specialCase == 'Be_10(n, 3n)' ) :
                    if( yosDatas[6] == [] ) : mYos = [[3, 1], [2, 2004] ]
                elif( specialCase == 'Be_11(n, 4n)' ) :
                    if( yosDatas[6] == [] ) : mYos = [[4, 1], [2, 2004] ]
                elif( ( mYos[-1] == [ 1, 4008 ] ) and ( yosDatas[16] != [] ) ) :                # Note_Be_9: Special case for breakup of Be_8 into two He_4's.
                    raise Exception( 'See note "Note_Be_9"' )
                    if( ( yosDatas[6] != [] ) and ( yosDatas[16] != [] ) ) :                    # This elif is probably not used anymore, replaced by
                        decayChannel = gnd.channels.NBodyOutputChannel( )                        # prior if.

                        He4Particle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, 2004 ) )
                        addDistributionDataAndRemove( He4Particle, 6, yosDatas )
                        decayChannel.products.add( decayChannel.products.uniqueLabel( He4Particle ) )

                        He4Particle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, 2004 ) ) 
                        addDistributionDataAndRemove( He4Particle, 16, yosDatas )
                        decayChannel.products.add( decayChannel.products.uniqueLabel( He4Particle ) )
                        decayChannelInfo = [ len( mYos ) - 1, decayChannel ]
                    else :
                        He4Particle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, 2004 ), multiplicity = 2 )
                        addDistributionDataAndRemove( He4Particle, 16, yosDatas )
                        decayChannel = gnd.channels.NBodyOutputChannel( )
                        decayChannel.products.add( decayChannel.products.uniqueLabel( He4Particle ) )
                        decayChannelInfo = [ len( mYos ) - 1, decayChannel ]
                        specialCase = None
                if( ( S == 0 ) and ( C in [ 11, 40, 41, 42, 44, 45 ] ) ) : specialCase = 'unknown level'
                outputChannel = makeNBodyChannelFrom_mYos( gnd.channels.NBodyOutputChannel, mYos, ZAsYos, yosDatas, I0.getQ( ), 
                    specialCase = specialCase, decayChannelInfo = decayChannelInfo, levelIndex = c_Or_s_Level )
        else :  # Reactions that are not two-body, two-body + break-up or discrete N-body except C = 71 to 74.
            if( I0.C == 1 ) :                                 # (n,total)
                outputChannel = None
            elif( I0.C == 15 ) :                              # Fission
                if( verbose > 0 ) : print 'Fission',
                outputChannel = gnd.channels.fissionChannel( fissionGenre = gnd.channels.fissionGenreTotal )
                outputChannel.Q.add( returnConstantQ( info.style, I0.getQ( ) ) )
                particle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, 1 ) )
                addDistributionDataAndRemove( particle, 1, yosDatas )
                outputChannel.products.add( outputChannel.products.uniqueLabel( particle ) )
                decayRates = delayedNeutrons.keys( )
                decayRates.sort( )
                decayRates.reverse( )
                for decayRate in decayRates:
                    gammasRemoved = []
                    for index, delayedDatum in enumerate( delayedNeutrons[decayRate] ) :
                        if( delayedDatum.yo == 7 ) : gammasRemoved.append( index )
                    if( len( gammasRemoved ) > 0 ) : print
                    while( len( gammasRemoved ) > 0 ) :
                        gammaRemoved = delayedNeutrons[decayRate].pop( gammasRemoved[-1] )
                        print '        remove delayed gamma data', `gammaRemoved`
                        gammasRemoved.pop( )
                    delayedData = { 1 : delayedNeutrons[decayRate] }
                    delayedParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, 1 ) )
                    Is_ = sorted( [ data.I for data in delayedData[1] ] )
                    if( Is_ == [ 7, 10, 13 ] ) :
                        if( verbose ) : print '  delayed neutron has I = 7, 10 and 13 but no others',
                    else :
                        addDistributionDataAndRemove( delayedParticle, 1, delayedData, promptNeutronParticle = particle )
                    outputChannel.products.add( outputChannel.products.uniqueLabel( delayedParticle ) )
                if( yosDatas[7] != [] ) :
                    gammaParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, 7 ) )
                    addDistributionDataAndRemove( gammaParticle, 7, yosDatas )
                    outputChannel.products.add( outputChannel.products.uniqueLabel( gammaParticle ) )
            elif( I0.C in [ 71, 72, 74 ] ) :
                if( verbose > 0 ) : print "7[124]",
                qualifier = { 71 : 'coherent', 72 : 'incoherent', 74 : 'pair production' }[I0.C]
                for idx in xrange( len( yosDatas[7] ) - 1, -1, -1 ) :
                    if( ( yosDatas[7][idx].I == 4 ) or ( ( I0.C == 74 ) and ( yosDatas[7][idx].I == 10 ) ) ) : del yosDatas[7][idx]
                firstParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, 7 ), attributes = { 'scattering' : qualifier } )
                addDistributionDataAndRemove( firstParticle, 7, yosDatas )
                if( C == 74 ) : 
                    raise Exception( 'A distribution class is needed here' )
                secondParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, residuals[0] ) )
                if( I0.C in [ 71, 72 ] ) :
                    outputChannel = gnd.channels.twoBodyOutputChannel( )
                else :
                    outputChannel = gnd.channels.NBodyOutputChannel( )
                outputChannel.Q.add( returnConstantQ( info.style, I0.getQ( ) ) )
                outputChannel.products.add( outputChannel.products.uniqueLabel( firstParticle ) )
                outputChannel.products.add( outputChannel.products.uniqueLabel( secondParticle ) )
            elif( I0.C == 73 ) :
                if( verbose > 0 ) : print "73",
                firstParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, 9 ), attributes = { 'scattering' : 'photo-electric' } )
                secondParticle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, residuals[0] ) )
                outputChannel = gnd.channels.NBodyOutputChannel( )
                outputChannel.Q.add( returnConstantQ( info.style, I0.getQ( ) ) )
                outputChannel.products.add( outputChannel.products.uniqueLabel( firstParticle ) )
                outputChannel.products.add( outputChannel.products.uniqueLabel( secondParticle ) )
            else :
                if( ( I0.C >= 50 ) and ( I0.C < 57 ) ) : 
                    if( verbose > 0 ) : print "production",
                    if( yos[0] != -4 ) : raise Exception( 'Does not appear to be production channel. Internal error. C = %d, yos = %s' % ( C, `yos` ) )
                    multiplicity = gnd.productData.multiplicity.partialProduction( info.style )
                    particle = toGNDMiscModule.newGNDParticle( info, toGNDMiscModule.getTypeName( info, yos[1] ), multiplicity = multiplicity )
                    addDistributionDataAndRemove( particle, ZAsYos[yos[1]], yosDatas )
                    outputChannel = gnd.channels.productionChannel( )
                    outputChannel.products.add( outputChannel.products.uniqueLabel( particle ) )
                else :
                    if( verbose > 0 ) : print
                    print "NONO:", I0
                    continue
        if( verbose > 0 ) : print
        s = str( outputChannel )
        if( S == 2 ) :
            s += ' S2'
        elif( C == 9 ) :
            s += ' n+i'
        elif( C == 55 ) :
            s += ' discreteEnergy = %s' % I0.getX1( )
        if( s in channelList ) :
            endlmisc.printWarning( '    WARNING: %s already in list for %d: C=%s, S=%s, X1=%s\n' % ( s, self.ZA, C, S, I0.getX1( ) ) )
        else :
            channelList.append( s )

        if( I12 is not None ) : 
            axes = gnd.channelData.Q.XYs1d.defaultAxes( energyUnit = 'MeV', QUnit = 'MeV' )
            Q = gnd.channelData.Q.XYs1d( label = info.style, data = I12.data, axes = axes, 
                    accuracy = ENDL_Accuracy, interpolation = getXYInterpolation( I12 ) )
            outputChannel.Q.add( Q )

        axes = gnd.reactionData.crossSection.defaultAxes( energyUnit = 'MeV' )
        crossSection = gnd.reactionData.crossSection.XYs1d( data = I0.data, label = info.style, 
                interpolation = getXYInterpolation( I0 ), axes = axes, accuracy = ENDL_Accuracy )

        CCounts = len( self.findDatas( C = I0.C, I = 0 ) )
        if( isThermalNeutronScatteringLaw ) :
            MT = TNSL_MT
        else :
            MT = ENDLCS_To_ENDFMT.getMTFromCS( I0.C, I0.S, CCounts = CCounts )
        if( I20 is not None ) : endlmisc.printWarning( '    WARNING: I20 URR data not currently supported.' )

        dateStr = str( maxReactionDate )
        dateStr = '%s-%s-%s' % ( dateStr[:4], dateStr[4:6], dateStr[6:] )
        if( C == 1 ) :
            summands = [ gnd.sums.add( link = r.crossSection ) for r in reactionSuite.reactions ]    # all other reactions have already been read in
            reaction = gnd.sums.crossSectionSum( name = "total", label = str( iChannel ), ENDF_MT = MT,
                    summands = gnd.sums.listOfSummands( summandList = summands ), date = dateStr )
            reaction.Q.add( returnConstantQ( info.style, I0.getQ( ) ) )
        else:
            reaction = reactionModule.reaction( outputChannel, "%s" % iChannel, ENDF_MT = MT, date = dateStr )
        reaction.crossSection.add( crossSection )

        if( S == 8 ) :
            reaction.process = "ENDL_S8"
        elif( C == 8 ) :
            reaction.process = "largeAngleCoulombScattering"
        elif( C == 9 ) :
            reaction.process = "nuclearInterferenceForLargeAngleCoulombScattering"
        elif( specialCase in ( 'Be_9(n, 2n)', 'Be_10(n, 3n)', 'Be_11(n, 4n)' ) ) :
            reaction.process = specialCase

        if( C == 1 ) :
            reactionSuite.sums.add( reaction )
        elif( C == 55 ) :           # FIXME. Do not set __class__ and moniker.
            reaction.__class__ = partialGammaProductionModule.partialGammaProduction
            reaction.moniker = partialGammaProductionModule.partialGammaProduction.moniker
            reactionSuite.addPartialGammaProduction( reaction, iChannel )
        else:
            reactionSuite.reactions.add( reaction )
        for yo in yosDatas :
            if( len( yosDatas[yo] ) ) :
                message = [ `yosData` for yosData in yosDatas[yo] ]
                endlmisc.printWarning( '    WARNING: unused data for ZA = %d, yo = %d\n        %s' % ( self.ZA, yo, '\n        '.join( message ) ) )

    if( isThermalNeutronScatteringLaw ) :
        TNSLStyle = gnd.miscellaneous.style( 'ENDL_TNSL', 
                attributes = { 'target' : targetName.split( '_' )[0], 'energyMax' : TNSL_EMax } )
        reactionSuite.styles.add( TNSLStyle )

    documentation = self.getDocumentation()
    if( documentation is None ) : documentation = "No documentation!"
    docs = gnd.documentation.documentation( "ENDL", documentation )
    reactionSuite.addDocumentation( docs )

    for reaction in reactionSuite.reactions : toGNDMiscModule.addUnspecifiedDistributions( info, reaction.outputChannel )
    for production in reactionSuite.productions : toGNDMiscModule.addUnspecifiedDistributions( info, production.outputChannel )

    return( reactionSuite )
