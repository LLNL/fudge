# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
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

from pqu import physicalQuantityWithUncertainty
from fudge.core.math.xData import axes, XYs, W_XYs
from fudge import gnd
from fudge.gnd.productData import distributions
from fudge.structure import masses as massModule
from fudge.legacy.endl import (endl2, endlmisc, endl_C, endl_y, endlIClasses,)
import endf_endl

ENDL_Accuracy = 1e-3
FUDGE_EPS = 1e-8

def getXYInterpolation( data ) :

    if( data.columns != 2 ) : raise Exception( 'Only 2 column data supported, # of columns = %s for %s' % ( data.columns, str( data ) ) )
    if( data.interpolation == 0 ) : return( axes.linearToken, axes.linearToken, )
    if( data.interpolation == 1 ) : return( axes.logToken,    axes.linearToken )
    if( data.interpolation == 2 ) : return( axes.linearToken, axes.logToken )
    if( data.interpolation == 3 ) : return( axes.logToken,    axes.logToken )
    raise Exception( 'Unsupported interpolation = "%s" for %s' % ( data.interpolation, str( data ) ) )

class infos :

    def __init__( self, xenslIsotopes = None, transportables = [] ) :

        self.masses = massModule
        self.xenslIsotopes = xenslIsotopes
        if( xenslIsotopes is None ) : self.xenslIsotopes = {}
        self.transportables = transportables
        self.particleList = gnd.xParticleList.xParticleList( ) 
        self.newIndices( )
        self.ZAMasses = {}

    def newIndices( self ) :

        self.indices = tokenAndIndices( )

    def setReactionSuite( self, reactionSuite )  :

        self.reactionSuite = reactionSuite

def getTypeName( info, ZA, name = None, levelIndex = None, level = 0., levelUnit = 'MeV' ) :
    """Returns the name for this ZA and level if present. Returned name is of the form Am_242, Am_242_m1, Am_242_e2.
    levelIndex must be None, 'c', 's' or an integer > 0."""

    if( levelIndex not in [ 'c', 's', None ] ) :
        levelIndex = int( levelIndex )
        if( levelIndex < 0 ) : raise Exception( 'levelIndex = %d must be >= 0' % levelIndex )
    if( ZA == 17 ) : ZA = 7
    particleType, particleQualifiers = gnd.xParticle.particleType_Nuclear, {} # Default values.
    if( not ( name is None ) ) : 
        pass
    elif( ZA == 7 ) : 
        particleType, name = gnd.xParticle.particleType_Photon, 'gamma'
    elif( ZA == 9 ) :
        particleType, name = gnd.xParticle.particleType_Lepton, 'electron'
    else :
        name = endl2.endlToGNDName( ZA )
        if( levelIndex not in [ None, 0 ] ) :
            if( level < 0. ) : raise Exception( 'Negative value = %s for continuum is not allowed' % level )
            if( levelIndex in [ 'c', 's' ] ) :        # to continuum or sum of all levels.
                particleQualifiers['level'] = gnd.miscellaneous.undefinedLevel( physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( level, levelUnit ) )
                name += "_%s" % levelIndex
            else :
                particleQualifiers['level'] = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( level, levelUnit )
                if( levelIndex != 0 ) : name += '_e%d' % levelIndex
    if( not( info.particleList.hasParticle( name ) ) ) :
        if( ZA in info.ZAMasses ) :
            if( not( info.ZAMasses[ZA] is None ) and ( info.ZAMasses[ZA] < 0 ) ) : info.ZAMasses[ZA] *= -1
            mass = info.ZAMasses[ZA]
        else :
            mass = info.masses.getMassFromZA( ZA )
            if( mass == None ) :                                        # If mass is not found, let's attempt to make something up that is close.
                A = ZA % 1000
                for i in xrange( 1, A ) :
                    massLower = info.masses.getMassFromZA( ZA - i )
                    massHigher = info.masses.getMassFromZA( ZA + i )
                    if( not( massLower is None ) ) :
                        mass = massLower + i * info.masses.getMassFromZA( 1 )
                        break
                    elif( not( massHigher is None ) ) :
                        mass = massHigher - i * info.masses.getMassFromZA( 1 )
                        break
                if( mass == None ) : raise Exception( 'Could not obtain mass for ZA = %s' % ZA )
                mass = round(mass,1)    # since we're extrapolating
        if( mass is None ) :
            massUnit = None
        else :
            massUnit = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( mass, 'amu' )
        if( name in info.transportables ) : particleQualifiers['transportable'] = 'true'
        if( not( levelIndex is None ) and ( levelIndex != 0 ) ) : particleQualifiers['levelIndex'] = levelIndex
        if( name in [ 'n', 'H1', 'H2', 'H3', 'He3', 'He4', 'gamma' ] ) : particleQualifiers['transportable'] = 'true'

        if( particleQualifiers.get( 'level' ) is not None ) :                   # Add new particle or level.
            # careful, always need to add ground state before excited level:
            baseName = name.split('_')[0]
            if( 'natural' in name ) :
                baseName = '_'.join( name.split('_')[:2] )
            elif( "FissionProduct" in name ) :
                baseName = '_'.join( name.split('_')[:2] )
            if( baseName not in info.particleList ) :
                info.particleList.addParticle( gnd.xParticle.xParticle( baseName, particleType, mass = massUnit, attributes={} ) )
            # now ground state is present, add excited state
            if( 'levelIndex' in particleQualifiers ) :
                index = particleQualifiers.pop('levelIndex')
            elif( particleQualifiers['level'].getValueAs( 'MeV' ) == '0.' ) :
                index = 0
            else :
                raise Exception( 'Index could not be determined for name = "%s", level = "%s", levelIndex = "%s"' % ( name, level, levelIndex ) )
            energy = particleQualifiers.pop( 'level' )
            particleOrLevel = gnd.xParticle.nuclearLevel( name, energy, index, attributes = particleQualifiers )
        else:
            particleOrLevel = gnd.xParticle.xParticle( name, particleType, mass = massUnit, attributes = particleQualifiers )
        info.particleList.addParticle( particleOrLevel )

    return( info.particleList.getParticle( name ) )

def newGNDParticle( info, particle, attributes = {}, multiplicity = 1, decayChannel = None, 
        ancestryName = 'product' ) :
    """The arguments to this function must match those of product.product excluding self and label."""

    label = info.indices.getLabel( particle.getToken( ) )
    return( gnd.product.product( particle, ancestryName = ancestryName, label = label, attributes = attributes,
        multiplicity = multiplicity, decayChannel = decayChannel ) )

class tokenAndIndices :

    suffixCharacters = 'abcdefghijklmnopqrstuvwxyz'

    def __init__( self ) :

        self.tokens = {}
        self.delayedFissionNeutron = False

    def getLabel( self, token ) :

        if( self.delayedFissionNeutron ) :
            if( token not in self.tokens ) : self.tokens[token] = -1
            self.tokens[token] += 1
            return( token + '__' + self.suffixCharacters[self.tokens[token]] )
        if( token not in self.tokens ) :
            self.tokens[token] = 0
            return( token )
        index = self.tokens[token]
        self.tokens[token] += 1
        suffix = ''
        while( True ) :
            index1, index = index % 26, index / 26
            suffix = self.suffixCharacters[index1] + suffix
            if( index == 0 ) : break
        return( token + '__' + suffix )

def returnConstantQ( Q_MeV ) :

    return( gnd.channelData.Q.component( gnd.channelData.Q.constant( physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( Q_MeV, 'MeV' ) ) ) )

def toGND( self, xenslIsotopes = None, verbose = 0, testing = False ) :
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
            if( S == None ) : S = data.S
            if( ( data.I == I ) and ( data.S == S ) ) : return( data )
        return( None )

    def getAllDistributionIs( yoDatas, IMax = 9, IExtras = [] ) :

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
                if( yoOld != None ) : mYos.append( [ n, yoOld ] )
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

            muProbabilityFrame = axes.labToken
            if( hasattr( data, 'muProbabilityFrame' ) ) : muProbabilityFrame = data.muProbabilityFrame
            axes_ = distributions.angular.pointwise.defaultAxes( energyUnit = 'MeV' )
            form = distributions.angular.pointwise( axes_, productFrame = muProbabilityFrame )
            axesMuP = axes.referenceAxes( form )
            for i, EMuP in enumerate( data.data ) :
                E, muP = EMuP
                form[i] = XYs.XYs( axesMuP, muP, accuracy = ENDL_Accuracy, value = E, index = i, parent = form )
            return( form )

        def energy( data ) :

            axes_ = distributions.energy.pointwise.defaultAxes( energyUnit = 'MeV' )
            form = distributions.energy.pointwise( axes_, productFrame = axes.labToken )
            axes_xy = axes.axes( )
            axes_xy[0] = axes_[1]
            axes_xy[1] = axes_[2]
            for E, EpP in data.data[0][1] : form.append( XYs.XYs( axes_xy, EpP, accuracy = ENDL_Accuracy, value = E ) )
            return( form )

        def LLNLLegendrePointwise( data ) :

            axes_ = distributions.Legendre.LLNLPointwise.defaultAxes( )
            form = distributions.Legendre.LLNLPointwise( axes_, productFrame = axes.labToken )
            axesRefXY = axes.referenceAxes( form )
            for l, EEpPs in data.data :
                EEpP = distributions.Legendre.LLNLPointwiseEEpP( int( l ) )
                xPrior = -1
                for index, E_EpPs in enumerate( EEpPs ) :
                    E, EpPs = E_EpPs
                    for xy in EpPs :
                        if( xy[0] == xPrior ) : xy[0] *= ( 1 + FUDGE_EPS )
                        xPrior = xy[0]
                    EEpP.append( XYs.XYs( axesRefXY, EpPs, ENDL_Accuracy, value = E, index = index ) )
                form.append( EEpP )
            return( form )

        def LLNLLegendrePointwise_L0_only( data ) :
            """Has only L=0 component, can be converted to uncorrelated angular (isotropic) and energy distributions."""

            isotropic = distributions.angular.isotropic( productFrame = axes.labToken )
            angular = distributions.angular.component( )
            angular.addForm( isotropic )
            angular.nativeData = isotropic.name

            # energy spectrum:
            axes_ = distributions.energy.pointwise.defaultAxes( energyUnit = 'MeV' )
            pw = distributions.energy.pointwise( axes_, productFrame = axes.labToken )
            for index, (energy, probability_list) in enumerate(data.data[0][1]):
                axes_ = axes.referenceAxes( parent=pw )
                pw.append( XYs.XYs( axes_, probability_list, accuracy=ENDL_Accuracy, index=index, value=energy, parent=pw) )
            energy = distributions.energy.component()
            energy.addForm( pw )
            energy.nativeData = pw.name
            return distributions.uncorrelated.component( angular, energy )

        def LLNLAngularEnergy( data ) :

            axes_ = distributions.angularEnergy.pointwise.defaultAxes( energyUnit = "MeV", energyInterpolationQualifier = axes.unitBaseToken,
                energy_outUnit = 'MeV', probabilityUnit = '1/MeV' )
            axes_[3].label = 'P(energy_out|energy_in,mu)'
            form = distributions.angularEnergy.pointwise( axes_, productFrame = axes.labToken )
            axesW_XY = axes.referenceAxes( form, dimension = 3 )
            axesXY = axes.referenceAxes( form )
            for i, E_MuEpPs in enumerate( data.data ) :
                w_xys = W_XYs.W_XYs( axesW_XY, index = i, value = E_MuEpPs[0] )
                for j, Mu_EpPs in enumerate( E_MuEpPs[1] ) :
                    w_xys.append( XYs.XYs( axesXY, Mu_EpPs[1], accuracy = ENDL_Accuracy, value = Mu_EpPs[0], index = j, parent = w_xys ) )
                form.append( w_xys )
            return( form )

        def addForm( component, I, datas, func ) :

            for i, data in enumerate( datas ) :
                if( data.I == I ) :
                    component.addForm( func( data ) )
                    del datas[i]
                    return
            raise Exception( 'I = %d not in datas' % I )

        if( particle.getName( ) == 'gamma' ) : particle.addAttribute( 'ENDFconversionFlag', 'MF6' )
        if( yo not in yosDatas ) : return
        datas = yosDatas[yo]
        Is = getAllDistributionIs( datas, IMax = 4, IExtras = [ 941, 942 ] )
        component = None
        if( Is == [] ) :
            if( len( datas ) ) :
                if( ( len( datas ) == 1 ) and ( datas[0].S == 7 ) and ( datas[0].I == 7 ) ) :   # Special case for endl I, S = 7, 7 with no distribution data.
                    component = distributions.base.referenceComponent( promptNeutronParticle.distributions )
                else :
                    print
                    for data in datas : print data
                    raise Exception( "Unsupported data for particle %s genre." % ( `Is`, particle.getName( ) ) )
        elif( ( 941 in Is ) or ( 942 in Is ) ) :
            raise Exception( 'I = 941 or 942 not supported' )
        elif( Is == [ 1 ] ) : 
            if( datas[0].S == 3 ) :
                component = distributions.angular.component( nativeData = distributions.base.pointwiseFormToken )
                addForm( component, 1, datas, angular )
            else :
                for data in datas :
                    if( data.I == 1 ) : data.muProbabilityFrame = axes.centerOfMassToken
                component = distributions.angular.twoBodyComponent( distributions.base.pointwiseFormToken )
                addForm( component, 1, datas, angular )
        elif( Is == [ 4 ] ) :
            for i in range(len(datas)):
                if datas[i].I==4: break
            if (len(datas[i].data) == 1 and datas[i].data[0][0] == 0 ):
                component = LLNLLegendrePointwise_L0_only( datas[i] )
                del datas[i]
            else:
                component = distributions.Legendre.component( distributions.base.LLNLLegendrePointwiseFormToken )
                addForm( component, 4, datas, LLNLLegendrePointwise )
        elif( Is == [ 1, 3 ] ) :
            angularComponent = distributions.angular.component( nativeData = distributions.base.pointwiseFormToken )
            addForm( angularComponent, 1, datas, angular )

            angularEnergyComponent = distributions.angularEnergy.LLNLComponent( distributions.base.pointwiseFormToken )
            addForm( angularEnergyComponent, 3, datas, LLNLAngularEnergy )

            component = distributions.angularEnergy.LLNL_withAngularComponent( angularComponent, angularEnergyComponent )
        elif( Is == [ 1, 4 ] ) :
            angularComponent = distributions.angular.component( nativeData = distributions.base.pointwiseFormToken )
            addForm( angularComponent, 1, datas, angular )

            energyComponent = distributions.energy.component( distributions.base.pointwiseFormToken )
            addForm( energyComponent, 4, datas, energy )

            component = distributions.uncorrelated.component( angularComponent, energyComponent )
        else :
            raise Exception( "I's = %s not supported for particle %s genre." % ( `Is`, particle.getName( ) ) )
        for i in xrange( len( datas ) - 1, -1, -1 ) :
            data = datas[i]
            if( data.I in [ 7, 9, 10, 13, 941, 942 ] ) :
                if( data.I == 7 ) :
                    if( data.S == 7 ) : 
                        particle.addAttribute( 'emissionMode', 'delayed' )
                        particle.addAttribute( 'decayRate', physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( data.getX1( ), "1/s" ) )
                    else :
                        particle.addAttribute( 'emissionMode', 'prompt' )
                xInterpolationStr, yInterpolationStr = getXYInterpolation( data )
                if( data.I in [ 7, 9 ] ) : 
                    axes_ = gnd.productData.multiplicity.pointwise.defaultAxes( energyUnit = 'MeV', energyInterpolation = xInterpolationStr, \
                            multiplicityInterpolation = yInterpolationStr )
                    particle.setMultiplicity( gnd.productData.multiplicity.pointwise( axes_, data.data, ENDL_Accuracy ) )
                else :          # for I = 10 and 13
                    if( data.I == 10 ) :
                        axes_ = gnd.productData.energyDeposition.pointwise.defaultAxes( energyUnit = 'MeV', energyInterpolation = xInterpolationStr, \
                            energyDepositionInterpolation = yInterpolationStr )
                        dataComponent = gnd.productData.energyDeposition.component( gnd.tokens.pointwiseFormToken )
                        dataComponent.addForm( gnd.productData.energyDeposition.pointwise( axes_, data.data, ENDL_Accuracy ) )
                    elif( data.I == 13 ) :
                        axes_ = gnd.productData.momentumDeposition.pointwise.defaultAxes( energyUnit = 'MeV', energyInterpolation = xInterpolationStr, \
                            momentumDepositionInterpolation = yInterpolationStr, momentumDepositionUnit = 'MeV/c' )
                        dataComponent = gnd.productData.momentumDeposition.component( gnd.tokens.pointwiseFormToken )
                        dataComponent.addForm( gnd.productData.momentumDeposition.pointwise( axes_, data.data, ENDL_Accuracy ) )
                    particle.addData( dataComponent )
                del datas[i]
        if( not component is None ) :
            particle.addDistributionComponent( component )
            particle.distributions.setNativeData( component.moniker )

    def makeNBodyChannelFrom_mYos( channelClass, mYos, ZAsYos, yosDatas, Q_MeV, specialCase = None, decayChannelInfo = [], levelIndex = None ) :
        """This routine does not support yo's (including the residual) in an excited state. This is consistent with endl99.
        This routine does not work if decayChannelInfo has data (see note Note_Be_9) due to added gndPath logic (gndPath has been
        removed and may be so should this last sentence)."""

        if( ( len( yosDatas[7] ) > 0 ) and ( hasattr( yosDatas[7][0], 'specialCase' ) ) ) :
            specialCase = yosDatas[7][0].specialCase
            if( specialCase == 'Am_242_m1' ) : level = yosDatas[7][0].level
        s = ' '
        channel = channelClass( returnConstantQ( Q_MeV ) )
        i = 0
        n = len( mYos ) - 1
        for m, ZA in mYos :
            decayChannel = None
            if( len( decayChannelInfo ) > 0 ) :
                if( decayChannelInfo[0] == i ) : decayChannel = decayChannelInfo[1]
            if( m > 1 ) :
                particle = newGNDParticle( info, getTypeName( info, ZA ), multiplicity = m, decayChannel = decayChannel )
            elif( ( specialCase == 'unknown level' ) and ( i == 1 ) ) :
                particle = newGNDParticle( info, getTypeName( info, ZA, level = 0.0, levelIndex = levelIndex ), decayChannel = decayChannel )
            elif( ( specialCase == 'Am_242_m1' ) and ( ZA == 95242 ) ) :
                particle = newGNDParticle( info, getTypeName( info, ZA, level = level, levelIndex = 2 ), decayChannel = decayChannel )
            else :
                particle = newGNDParticle( info, getTypeName( info, ZA ), decayChannel = decayChannel )
            s = ' + '
            if( ZA in ZAsYos ) :
                if( ( specialCase == 'S8' ) and ( ZA == 1 ) ) :
                    yo = 11
                else :
                    yo = ZAsYos[ZA]
                    if( ( i == n ) and ( yosDatas[yo] == [] ) ) : yo += 10
                addDistributionDataAndRemove( particle, yo, yosDatas )
            channel.addProduct( particle )
            i += 1
        return( channel )

    if( self.ZA == 95241 ) :
        I0s = self.findDatas( C = 46, S = 0, I = 0 )
        if( len( I0s ) == 2 ) :
            I0m1, I0g = I0s
            if( I0m1.getQ( ) > I0g.getQ( ) ) : I0m1, I0g = I0g, I0m1
            level = I0g.getQ( ) - I0m1.getQ( )
            dataM1s = self.findDatas( C = 46, S = 0, Q = I0m1.getQ( ) )
            for dataM1 in dataM1s :
                dataM1.specialCase = 'Am_242_m1'
                dataM1.level = level
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
    if( not hasattr( self, 'evaluation' ) ) :
        if( '/yi0' in self.source ) :
            evaluation = os.path.basename( self.source.split( '/yi0' )[0] )
        else :
            evaluation = os.path.basename( self.source )
    evaluation = "ENDL"                     # Override above for now. Note, this is used in reactionSuite.toENDF6
    evaluatedStyle = gnd.miscellaneous.style( 'evaluated', attributes = { 'version' : evaluation } )
    info = infos( xenslIsotopes, transportables = [ 'n', 'H1', 'H2', 'H3', 'He3', 'He4', 'gamma' ] )
    projectile = getTypeName( info, endl2.yoToZA( self.yi ) )

    CsAndSs, residualExcitationIndexLevels = {}, {}
    if( level is not None ) : residualExcitationIndexLevels[self.ZA] = [ [ None, level ] ]
    reactionsDatas = self.findReactionsDatas( )
    compoundZA = endl2.yoToZA( self.yi ) + self.ZA
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
        indexLevels = {}
        for index, indexAndLevel in enumerate( residualExcitationIndexLevels[residualZA] ) : 
            indexLevel = index + 1
            if( indexAndLevel[1] == 0 ) : indexLevel = 0
            indexLevels[indexAndLevel[1]] = indexLevel
        residualExcitationIndexLevels[residualZA] = indexLevels

    levelIndex = None
    if( level is not None ) : levelIndex = residualExcitationIndexLevels[self.ZA][level]
    target = getTypeName( info, self.ZA, level = level, levelIndex = levelIndex )

    reactionSuite = gnd.reactionSuite.reactionSuite( projectile, target, particleList = info.particleList, style = evaluatedStyle )
    info.setReactionSuite( reactionSuite )
    reactionSuite.setAttribute( 'temperature', physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( I0.getTemperature( ), 'MeV' ) )
    if( self.ZA == 95242 ) :
        reactionSuite.addNuclearMetaStableAlias( target.getName( ).split( '_' )[0], target.getName( ), 1 )

    I20s = []                                           # Some I20 data has bad Q values, so let's ignore all I20 data.
    for i, reactionDatas in enumerate( reactionsDatas ) :
        for j, reactionData in enumerate( reactionDatas ) :
            if( reactionData.I == 20 ) : I20s.append( [ i, j ] )
    I20s.reverse( )
    for i, j in I20s :
        del reactionsDatas[i][j]
        if( len( reactionsDatas[i] ) == 0 ) : del reactionsDatas[i]

    channelList, maxDate, delayedNeutrons = [], 0, {}
    for iChannel, reactionDatas in enumerate( reactionsDatas ) :
        C, S = reactionDatas[0].C, reactionDatas[0].S   # Phase 1
        c_Or_s_Level = 'c'
        if( ( S == 0 ) and ( 1 not in CsAndSs[C] ) ) : c_Or_s_Level = 's'
        info.newIndices( )
        maxReactionDate = 0
        I0, I12, I20 = None, None, None
        specialCase = None
        for yo in yosDatas : yosDatas[yo] = []
        for data in reactionDatas :                         # yosDatas should only contain I = 1, 3, 4, 7, 9, 10 and 13 type data.
            data.frame = ( endlmisc.getNumberOfColumns_( data.I, '' ) * ( '%s,' % axes.labToken ) )[:-1]
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
                    if( ( yo1I10 != None ) and ( yo12I10 != None ) ) :
                        yosDatas[11].append( yo12I10 )
                        yo1I10.set( yo1I10 - yo12I10 )
        for yo in yosDatas :
            if( ( len( yosDatas[yo] ) == 1 ) and ( yosDatas[yo][0].I == 10 ) ) : yosDatas[yo] = []
        if( I0 == None ) : raise Exception( 'Missing cross section data for reaction C = %d S = %d.' % ( C, S ) )   # End of Phase 1

        if( ( C == 12 ) and ( self.ZA == 4009 ) and ( self.yi == 1 ) ) :  # Special case for Be_9(n,2n) which should really be Be_9(n,2n a)
            specialCase = 'Be_9(n, 2n)'
            yos = ( 1, 1, 2004 )
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
            residuals = [ self.ZA ]
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
            residuals = [ self.ZA - productsZA ]        # End of phase 2

        if( 4008 in residuals ) :                                       # Special case for Be_8 breakup into two alphas
            if( ( 2004 not in yos ) and ( yosDatas[6] != [] ) and ( yosDatas[16] == [] ) ) :
                yosDatas[16] = yosDatas[6]
                yosDatas[6] = []

        if( 99000 < self.ZA < 99200 ) :
            residuals = [ self.ZA ]
        elif( ( C == 30 ) and ( endl2.yoToZA( self.yi ) in [ 1002, 1003 ] ) and ( self.ZA in [ 1002, 1003 ] ) ) :   # Special treatment for C = 30 as 
            if( yosDatas[6] != [] ) : raise Exception( 'C = 30 as yo = 6 data' )                    # yos should be [ g, n ] and not [ g, n, a ].
            yos = ( 7, 1 )
            residuals = ( 2004, )

        if( verbose > 0 ) : print '   ', C, S,
        if( ( yos[0] > 0 ) and ( C not in [ 71, 72, 73, 74 ] ) ) :  # Reactions that are two-body, two-body + break-up or discrete N-body except C = 71 - 74.
            outputChannel = None
            Is = getAllDistributionIs( yosDatas[ZAsYos[yos[0]]] )
            if( Is == [ 1 ] ) :     # Two-body initial product state.
                if( verbose > 0 ) : print '2Body',
                firstParticle = newGNDParticle( info, getTypeName( info, yos[0] ) )    # Assume first particle is in the ground state.
                addDistributionDataAndRemove( firstParticle, ZAsYos[yos[0]], yosDatas )
                Q = returnConstantQ( I0.getQ( ) )
                if( len( yos ) == 1 ) :                                 # Simple two-body, yos[0] + Residual [ + gamma].
                    resYo = residuals[0]
                    if( resYo in ZAsYos ) : resYo = ZAsYos[resYo] + 10
                    if( S == 1 ) :
                        Q = returnConstantQ( I0.getQ( ) - I0.getX1( ) )
                        resParticle = newGNDParticle( info, getTypeName( info, residuals[0] ) )
                        addDistributionDataAndRemove( resParticle, resYo, yosDatas )

                        decayChannel = None
                        if( 7 in yos ) : print '7 in yos'
                        if( ( 7 not in yos ) and ( I0.getX1( ) != 0. ) ) :
                            decayChannel = gnd.channels.NBodyDecayChannel( returnConstantQ( I0.getX1( ) ) )     # Assume residual returns to ground
                            decayChannel.addProduct( resParticle )                                              # state, hence Q = getX1( )
                            gammaParticle = newGNDParticle( info, getTypeName( info, 7 ) )
                            if( yosDatas[7] != [] ) :
                                addDistributionDataAndRemove( gammaParticle, 7, yosDatas )
                            decayChannel.addProduct( gammaParticle )

                        level = I0.getX1( )
                        levelIndex = residualExcitationIndexLevels[residuals[0]][level]
                        if( levelIndex is None ) : level = None
                        secondParticle = newGNDParticle( info, getTypeName( info, residuals[0], level = level, levelIndex = levelIndex ), 
                            decayChannel = decayChannel )
                    else :
                        decayChannel, level, levelIndex, resName = None, None, None, None
                        if( S == 0 ) :
                            level = I0.getELevel( )
                            if( C not in [ 8, 9, 10 ] ) :
                                levelIndex = c_Or_s_Level
                                if( ( yosDatas[7] != [] ) and ( residuals[0] != 7 ) ) :
                                    decayChannel = gnd.channels.NBodyDecayChannel( returnConstantQ( I0.getX1( ) ) )
                                    resParticle = newGNDParticle( info, getTypeName( info, residuals[0] ) )       # Residual in ground state.
                                    addDistributionDataAndRemove( resParticle, resYo, yosDatas )
                                    decayChannel.addProduct( resParticle )
                                    gammaParticle = newGNDParticle( info, getTypeName( info, 7 ) )
                                    addDistributionDataAndRemove( gammaParticle, 7, yosDatas )
                                    decayChannel.addProduct( gammaParticle )
                            elif( I0.getELevel( ) != 0 ) :            # Elastic scattering off an isomer
                                levelIndex = residualExcitationIndexLevels[residuals[0]][level]
                                resName = self.name
                            if( ( level == 0.0 ) and ( levelIndex is None ) ) : level = None
                        if( residuals[0] <= 2004 ) : level, levelIndex = None, None
                        if( self.name == 'Am242_m1' ) : resName = reactionSuite.aliases[self.name].getValue( )
                        secondParticle = newGNDParticle( info, getTypeName( info, residuals[0], level = level, levelIndex = levelIndex, name = resName ), 
                            decayChannel = decayChannel )
                        addDistributionDataAndRemove( secondParticle, resYo, yosDatas )
                else :                                          # Two-body with residual breaking up.
                    if( verbose > 0 ) : print 'w/breakup',
                    if( self.A == 0 ) : raise Exception( 'For C = %d, breakup of product not supported for two-body reaction for natural target.' )
                    mYos = getMultiplicityYos( self, yos[1:], residuals, yosDatas )
                    residualZA = self.ZA + endl2.yoToZA( self.yi ) - yos[0]
                    levelIndex, level = c_Or_s_Level, 0.0
                    if( I0.S == 1 ) :
                        level = I0.getX1( )

                        levelIndex = residualExcitationIndexLevels[residualZA][level]
                    elif( I0.S == 8 ) : 
                        levelIndex, level = 1, I0.getX1( )
                    Q_MeV = endl2.reactionQByZAs( [ endl2.yoToZA( self.yi ), self.ZA ], [ yos[0], residualZA ] ) - level
                    residualName = getTypeName( info, residualZA, level = level, levelIndex = levelIndex )
                    decayChannel = makeNBodyChannelFrom_mYos( gnd.channels.NBodyDecayChannel, mYos, ZAsYos, yosDatas, I0.getQ( ) - Q_MeV, 
                        specialCase = specialCase )
                    secondParticle = newGNDParticle( info, residualName, decayChannel = decayChannel )
                    Q = returnConstantQ( Q_MeV )
                outputChannel = gnd.channels.twoBodyOutputChannel( Q )
                outputChannel.addProduct( firstParticle )
                outputChannel.addProduct( secondParticle )
            else :                                              # Mutli-particle breakup (i.e., not two body).
                if( verbose > 0 ) : print 'NBody',
                mYos = getMultiplicityYos( self, yos, residuals, yosDatas )
                gammaPresent = False
                for m, yo in mYos : gammaPresent = gammaPresent or ( yo == 7 )
                if( ( not gammaPresent ) and ( yosDatas[7] != [] ) ) : mYos.append( [ 1, 7 ] )
                decayChannelInfo = []
                if( specialCase == 'Be_9(n, 2n)' ) :                                            # Special case for endl Be_9(n, 2n)
                    if( yosDatas[6] == [] ) : mYos = [[2, 1], [2, 2004] ]
                elif( ( mYos[-1] == [ 1, 4008 ] ) and ( yosDatas[16] != [] ) ) :                # Note_Be_9: Special case for breakup of Be_8 into two He_4's.
                    raise Exception( 'See note "Note_Be_9"' )
                    if( ( yosDatas[6] != [] ) and ( yosDatas[16] != [] ) ) :                    # This elif is probably not used anymore, replaced by
                        decayChannel = gnd.channels.NBodyDecayChannel( )                        # prior if.

                        He4Particle = newGNDParticle( info, getTypeName( info, 2004 ) )
                        addDistributionDataAndRemove( He4Particle, 6, yosDatas )
                        decayChannel.addProduct( He4Particle )

                        He4Particle = newGNDParticle( info, getTypeName( info, 2004 ) ) 
                        addDistributionDataAndRemove( He4Particle, 16, yosDatas )
                        decayChannel.addProduct( He4Particle )
                        decayChannelInfo = [ len( mYos ) - 1, decayChannel ]
                    else :
                        He4Particle = newGNDParticle( info, getTypeName( info, 2004 ), multiplicity = 2 )
                        addDistributionDataAndRemove( He4Particle, 16, yosDatas )
                        decayChannel = gnd.channels.NBodyDecayChannel( )
                        decayChannel.addProduct( He4Particle )
                        decayChannelInfo = [ len( mYos ) - 1, decayChannel ]
                        specialCase = None
                if( ( S == 0 ) and ( C in [ 11, 40, 41, 42, 44, 45 ] ) ) : specialCase = 'unknown level'
                outputChannel = makeNBodyChannelFrom_mYos( gnd.channels.NBodyOutputChannel, mYos, ZAsYos, yosDatas, I0.getQ( ), 
                    specialCase = specialCase, decayChannelInfo = decayChannelInfo, levelIndex = c_Or_s_Level )
        else :  # Reactions that are not two-body, two-body + break-up or discrete N-body except C = 71 to 74.
            if( I0.C == 15 ) :                              # Fission
                if( verbose > 0 ) : print 'Fission',
                outputChannel = gnd.channels.fissionChannel( returnConstantQ( I0.getQ( ) ), fissionGenre = gnd.channels.fissionGenreTotal )
                particle = newGNDParticle( info, getTypeName( info, 1 ) )
                addDistributionDataAndRemove( particle, 1, yosDatas )
                outputChannel.addProduct( particle )
                decayRates = delayedNeutrons.keys( )
                decayRates.sort( )
                decayRates.reverse( )
                for decayRate in decayRates:
                    delayedData = { 1 : delayedNeutrons[decayRate] }
                    delayedParticle = newGNDParticle( info, getTypeName( info, 1 ) )
                    addDistributionDataAndRemove( delayedParticle, 1, delayedData, promptNeutronParticle = particle )
                    outputChannel.addProduct( delayedParticle )
                if( yosDatas[7] != [] ) :
                    gammaParticle = newGNDParticle( info, getTypeName( info, 7 ) )
                    addDistributionDataAndRemove( gammaParticle, 7, yosDatas )
                    outputChannel.addProduct( gammaParticle )
            elif( I0.C in [ 71, 72, 74 ] ) :
                if( verbose > 0 ) : print "7[124]",
                qualifier = { 71 : 'coherent', 72 : 'incoherent', 74 : 'pair production' }[I0.C]
                for i in xrange( len( yosDatas[7] ) - 1, -1, -1 ) :
                    if( ( yosDatas[7][i].I == 4 ) or ( ( I0.C == 74 ) and ( yosDatas[7][i].I == 10 ) ) ) : del yosDatas[7][i]
                firstParticle = newGNDParticle( info, getTypeName( info, 7 ), attributes = { 'scattering' : qualifier } )
                addDistributionDataAndRemove( firstParticle, 7, yosDatas )
                if( C == 74 ) : 
                    raise Exception( 'A distribution class is needed here' )
                    firstParticle.setGenre( distributions.base.NBodyGenre )
                secondParticle = newGNDParticle( info, getTypeName( info, residuals[0] ) )
                if( I0.C in [ 71, 72 ] ) :
                    outputChannel = gnd.channels.twoBodyOutputChannel( returnConstantQ( I0.getQ( ) ) )
                else :
                    outputChannel = gnd.channels.NBodyOutputChannel( returnConstantQ( I0.getQ( ) ) )
                outputChannel.addProduct( firstParticle )
                outputChannel.addProduct( secondParticle )
            elif( I0.C == 73 ) :
                if( verbose > 0 ) : print "73",
                firstParticle = newGNDParticle( info, getTypeName( info, 9 ), attributes = { 'scattering' : 'photo-electric' } )
                secondParticle = newGNDParticle( info, getTypeName( info, residuals[0] ) )
                outputChannel = gnd.channels.NBodyOutputChannel( returnConstantQ( I0.getQ( ) ) )
                outputChannel.addProduct( firstParticle )
                outputChannel.addProduct( secondParticle )
            else :
                if( ( I0.C >= 50 ) and ( I0.C < 57 ) ) : 
                    if( verbose > 0 ) : print "production",
                    if( yos[0] != -4 ) : raise Exception( 'Does not appear to be production channel. Internal error. C = %d, yos = %s' % ( C, `yos` ) )
                    attributes = {}
                    if( (I0.S == 3) and (yos[1] == 7) ) : attributes['discrete'] = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( I0.getX1( ), 'MeV' )
                    multiplicity = gnd.productData.multiplicity.partialProduction( )
                    particle = newGNDParticle( info, getTypeName( info, yos[1] ), attributes, multiplicity = multiplicity )
                    addDistributionDataAndRemove( particle, ZAsYos[yos[1]], yosDatas )
                    outputChannel = gnd.channels.productionChannel( )
                    outputChannel.addProduct( particle )
                else :
                    if( verbose > 0 ) : print
                    print "NONO:", I0
                    continue
        if( verbose > 0 ) : print
        if( testing ) : print 'label %2d %2s: %-100s  ### %s' % ( C, S, reactionSuite.inputParticlesToReactionString( suffix = ' --> ' ) + \
            str( outputChannel ), outputChannel.genre )
        s = str( outputChannel )
        if( S == 2 ) :
            s += ' S2'
        elif( C == 9 ) :
            s += ' n+i'
        elif( C == 55 ) :
            s += ' discreteEnergy = %s' % I0.getX1( )
        if( s in channelList ) :
            endlmisc.printWarning( 'WARNING: %s already in list for %d: C=%s, S=%s, X1=%s\n' % ( s, self.ZA, C, S, I0.getX1( ) ) )
        else :
            channelList.append( s )

        if( I12 is not None ) : 
            energyDependentQ = gnd.channelData.Q.component( )
            axes_ = gnd.channelData.Q.pointwise.defaultAxes( energyUnit = 'MeV', energyInterpolation = xInterpolationStr, QUnit = 'MeV', QInterpolation = yInterpolationStr )
            energyDependentQ.addForm( gnd.channelData.Q.pointwise( axes_, I12.data, ENDL_Accuracy ) )
            outputChannel.setQ( energyDependentQ )
        crossSection = gnd.reactionData.crossSection.component( nativeData = gnd.tokens.pointwiseFormToken )
        xInterpolationStr, yInterpolationStr = getXYInterpolation( I0 )
        axes_ = gnd.reactionData.crossSection.pointwise.defaultAxes( energyUnit = 'MeV', energyInterpolation = xInterpolationStr, \
            crossSectionInterpolation = yInterpolationStr )
        crossSection.addForm( gnd.reactionData.crossSection.pointwise( axes_, I0.data, ENDL_Accuracy ) )

        CCounts = len( self.findDatas( C = I0.C, I = 0 ) )
        MT = ENDLCS_To_ENDFMT.getMTFromCS( I0.C, I0.S, CCounts = CCounts )
        if( I20 is not None ) : endlmisc.printWarning( 'WARNING: I20 URR data not currently supported.' )
        reaction = gnd.reaction.reaction( outputChannel, "%s" % iChannel, ENDF_MT = MT, crossSection = crossSection )

        temperature = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( I0.getTemperature( ), 'MeV' )
        if( temperature != reactionSuite.getAttribute( 'temperature' ) ) : raise Exception( 'temperatures do not match for reaction %s' % `reaction` )

        # date to string:
        maxReactionDate = str(maxReactionDate)
        maxReactionDate = '%s-%s-%s' % (maxReactionDate[:4],maxReactionDate[4:6],maxReactionDate[6:])
        
        reaction.setAttribute( 'date', maxReactionDate )
        if( S == 8 ) : reaction.setAttribute( 'process', "ENDL_S8" )
        if( C == 8 ) :
            reaction.setAttribute( 'process', "largeAngleCoulombScattering" )
        elif( C == 9 ) :
            reaction.setAttribute( 'process', "nuclearInterferenceForLargeAngleCoulombScattering" )
        elif( specialCase == 'Be_9(n, 2n)' ) :
            reaction.setAttribute( 'process', specialCase )

        if( C == 55 ) :
            reaction.__class__ = gnd.partialGammaProduction.partialGammaProduction
            reaction.moniker = gnd.reactions.base.partialGammaProductionToken
            reactionSuite.addPartialGammaProduction( reaction, iChannel )
        else:
            reactionSuite.addReaction( reaction, iChannel )
        for yo in yosDatas :
            if( len( yosDatas[yo] ) ) : endlmisc.printWarning( 'WARNING: unused data for ZA = %d, yo = %d\n    %s' % ( self.ZA, yo, `yosDatas[yo]` ) )

    maxDate = str(maxDate)
    maxDate = '%s-%s-%s' % (maxDate[:4],maxDate[4:6],maxDate[6:])
    evaluatedStyle.setAttribute( 'date', maxDate )

    documentation = self.getDocumentation()
    if( documentation is None ) : documentation = "No documentation!"
    docs = gnd.documentation.documentation( "ENDL", documentation )
    reactionSuite.addDocumentation( docs )

    return( reactionSuite )
