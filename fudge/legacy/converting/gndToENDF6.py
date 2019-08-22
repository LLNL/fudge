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

from pqu import physicalQuantityWithUncertainty
from fudge.core.math.xData import axes, XYs
from fudge.core.utilities import fudgeZA
from fudge.gnd import tokens
from fudge.gnd.productData import multiplicity
from fudge.gnd import miscellaneous
import endfFormats
import fudge
import fudge.gnd.productData.distributions as distributionModule

def angularPointwiseEnergy2ENDF6( self, targetInfo ) :

    interpolation = axesToEndfInterpolationFlag( self.axes )
    if( targetInfo['doMF4AsMF6'] ) :
        interpolation += 10
        ENDFDataList = [ endfFormats.endfContLine( 0, self.value, interpolation, 0, 2 * len( self ), len( self ) ) ]
        ENDFDataList += endfFormats.endfNdDataList( self )
    else :
        ENDFDataList = [ endfFormats.endfContLine( 0, self.value, 0, 0, 1, len( self ) ) ]
        ENDFDataList += endfFormats.endfInterpolationList( [ len( self ), interpolation ] )
        ENDFDataList += endfFormats.endfNdDataList( self )
    return( ENDFDataList )

def gammasToENDF6_MF6( MT, endfMFList, flags, tempInfo, gammas ) :

    if( len( gammas ) == 1 ) :          # Check for LLNL legacy data
        distributions = gammas[0].distributions
        component = distributions.components[distributions.nativeData]
        if( component.nativeData == distributionModule.base.LLNLLegendrePointwiseFormToken ) :
            component.toENDF6( MT, endfMFList, flags, tempInfo )
            return
    LANG, LEP, discreteLEP, interpolationEIn, frame, nOrders, discreteGammasVsEIn, continuum = 1, -2, -2, 2, None, -1, {}, None
    massRatio = tempInfo['mass'] / ( tempInfo['mass'] + 1. )
    for gamma in gammas :
        isPrimary = 'primary' in gamma.attributes
        if( isPrimary ) :
            gammaEnergy = gamma.getAttribute(  'primary' ).getValueAs( 'eV' )
        elif( 'discrete' in gamma.attributes ) :
            gammaEnergy = gamma.getAttribute( 'discrete' ).getValueAs( 'eV' )
        else :
            if( continuum is None ) :
                continuum = gamma
                continue
            else :
                raise Exception( 'Multiple continuum gammas detected for MT=%s' % ( MT ) )
        distributions = gamma.distributions
        component = distributions.components[distributions.nativeData]
        form = component.forms[component.nativeData]
        if( gamma.multiplicity.nativeData == tokens.weightedReferenceFormToken ) : # special case, must convert back to Legendre form for ENDF
            form = gamma.multiplicity.forms[tokens.weightedReferenceFormToken].toENDF6( flags, tempInfo )
        if( discreteLEP < 0 ) :
            discreteLEP = 2
            if( len( form.axes ) == 3 ) :
                independent, dependent, dummy = form.axes[1].interpolation.getInterpolationTokens( )
                if( dependent == axes.flatToken ) : discreteLEP = 1
        if( frame is None ) : frame = form.getProductFrame( )
        incident_e_vals = []
        for EInCl in form :
            EIn = EInCl.value
            while EIn in incident_e_vals: EIn += 1e-8
            incident_e_vals.append( EIn )
            realGammaEnergy = gammaEnergy
            if( isPrimary ) : realGammaEnergy = -( gammaEnergy + massRatio * EIn )
            if( EIn not in discreteGammasVsEIn ) : discreteGammasVsEIn[EIn] = []
            if( nOrders < 0 ) : nOrders = len( EInCl.coefficients )
            if( nOrders != len( EInCl.coefficients ) ) :
                raise Exception( "nOrders = %d != len( EInCl.coefficients ) = %d for gammaEnergy = %s, incident energy = %s" % \
                    ( nOrders, len( EInCl.coefficients ), gammaEnergy, EIn ) )
            discreteGammasVsEIn[EIn].append( [ realGammaEnergy ] + EInCl.coefficients )
    continuumGammasVsEIn = {}
    if( continuum is not None ) :
        distributions = continuum.distributions
        if( distributions.nativeData != distributionModule.base.noneFormToken ) :
            component = distributions.components[distributions.nativeData]
            if( component.nativeData != distributionModule.base.noneFormToken ) :
                if( isinstance( component, distributionModule.uncorrelated.component ) and continuum.getAttribute( 'ENDFconversionFlag' ) == 'MF6' ) :
                    tempInfo["gammaToENDF6"] = True                             # Special case, must be converted back to Legendre form.
                    tempInfo['product'] = continuum
                    form = component.toENDF6( MT, endfMFList, flags, tempInfo )
                    tempInfo.dict.pop( "gammaToENDF6" )
                else :
                    form = component.forms[component.nativeData]
                if( form is None ) :                    # None when component.toENDF6 above was able to write it out.
                    if( frame is not None ) : raise Exception( "Oops, cannot handle continuum and ( discrete and/or primary )." )
                else :
                    interpolationEIn = axisToEndfInterpolationFlag( form.axes[0] )
                    if( frame is None ) : frame = form.getProductFrame( )
                    if( LEP < 0 ) :
                        LEP = 2
                        if( len( form.axes ) == 3 ) :
                            independent, dependent, dummy = form.axes[1].interpolation.getInterpolationTokens( )
                            if( dependent == axes.flatToken ) : LEP = 1
                    incident_e_vals = []
                    for energy_in in form :
                        EIn = energy_in.value
                        continuumGammas = []
                        for EoutCl in energy_in :
                            if( nOrders < 0 ) : nOrders = len( EoutCl )
                            if( nOrders != len( EoutCl ) ) :
                                raise Exception( "nOrders = %d != len( EoutCl ) = %d for continuum gammas, incident energy = %s" % \
                                    ( nOrders, len( EoutCl ), EIn ) )
                            continuumGammas += [ EoutCl.value ] + EoutCl.coefficients
                        while( EIn in incident_e_vals ) : EIn += 1e-8
                        incident_e_vals.append( EIn )
                        continuumGammasVsEIn[EIn] = continuumGammas
    if( frame is None ) : return
    Eins = discreteGammasVsEIn.keys( )
    for Ein in continuumGammasVsEIn :
        if( Ein not in Eins ) : Eins.append( Ein )
    wordsPerEout = nOrders + 1
    if( LEP < 0 ) : LEP = discreteLEP
    ENDFDataList = [ endfFormats.endfContLine( 0, 0, LANG, abs( LEP ), 1, len( Eins ) ) ]
    ENDFDataList.append( endfFormats.endfInterpolationLine( [ len( Eins ), interpolationEIn ] ) )
    for Ein in sorted( Eins ) :
        ND = 0
        gammaData = []
        if( Ein in discreteGammasVsEIn ) :
            discreteGammas = sorted( discreteGammasVsEIn[Ein], reverse = True )
            for discreteGamma in discreteGammas : gammaData += discreteGamma
            ND = len( gammaData ) / wordsPerEout
        if( Ein in continuumGammasVsEIn ) : gammaData += continuumGammasVsEIn[Ein]
        NW = len( gammaData )
        NEP = len( gammaData ) / wordsPerEout
        ENDFDataList += [ endfFormats.endfContLine( 0, Ein, ND, nOrders - 1, NW, NEP ) ]
        ENDFDataList += endfFormats.endfDataList( gammaData )

    toENDF6_MF6( MT, endfMFList, flags, tempInfo, 1, frame, ENDFDataList )

def toENDF6_MF6( MT, endfMFList, flags, tempInfo, LAW, frame, MF6 ) :

    tempInfo['MF6LCTs'].append( { axes.labToken : 1, axes.centerOfMassToken : 2 }[frame] )
    LIP = 0
    if tempInfo['product'].particle.name in tempInfo['metastables'] and not (50<=MT<=91 or 600<=MT<=850):
        # set LIP to metastable index of product UNLESS excited state of product is encoded in MT number:
        alias, = [ alias for alias in tempInfo['reactionSuite'].aliases.values() if
                alias.getValue() == tempInfo['product'].particle.name ]
        LIP = int( alias.getAttribute('nuclearMetaStable') )
    if( MT not in endfMFList[6] ) : endfMFList[6][MT] = []
    particleID = tempInfo[ 'zapID' ]    # get ZAP for the outgoing particle
    if( particleID == 'gamma' ) :
        ZAP = 0
    else :
        Z, A, suffix, ZAP = fudgeZA.gndNameToZ_A_Suffix( particleID )
    productMass = tempInfo['particleMass']
    if( ( particleID == 'n' ) and ( MT in [ 18 ] ) ) :                            # Special case when fission has MF = 6 data.
        nPoints, multiplicityList = 2, [ endfFormats.endfInterpolationLine( [ 2, 2 ] ) ] + \
            endfFormats.endfNdDataList( [ [ tempInfo['EMin'], 1. ], [ tempInfo['EMax'], 1. ] ] )
    else :
        nPoints, multiplicityList = tempInfo['multiplicity'].endfMultiplicityList( tempInfo )
    if( particleID == 'gamma' ) and 'primaryGammaEnergy' in tempInfo.keys():
        AWP = tempInfo.dict.pop('primaryGammaEnergy')
    else: AWP = productMass / tempInfo['neutronMass']
    ENDFHeaderList = [ endfFormats.endfContLine( ZAP, AWP, LIP, LAW, 1, nPoints ) ] + multiplicityList
    endfMFList[6][MT] += ENDFHeaderList + MF6

def upDateENDFMF8Data( endfMFList, tempInfo ) :

    for MT in tempInfo['MF8'] : upDateENDF_MT_MF8Data( MT, endfMFList, tempInfo )

def upDateENDF_MT_MF8Data( MT, endfMFList, tempInfo ) :

    from fudge.gnd.reactionData import crossSection

    MF8Channels = tempInfo['MF8'][MT]
    ZA, mass = tempInfo[ 'ZA' ], tempInfo[ 'mass' ]
    LIS, LISO, NO, NS = tempInfo['LIS'], tempInfo['LISO'], 1, len( MF8Channels )
    firstReaction = MF8Channels[0]
    residual = firstReaction.outputChannel[0]

    crossSection_ = firstReaction.crossSection.getNativeData( )  # LMF determines which MF section to write to.
    if( isinstance( crossSection_, crossSection.reference ) ) : 
        multiplicity_ = residual.multiplicity.getNativeData( )
        if(   isinstance( multiplicity_, multiplicity.constant ) ) :
            LMF = 3
        elif( isinstance( multiplicity_, multiplicity.reference ) ) :
            LMF = 6
        else :
            LMF = 9
    else :
        LMF = 10

    level = residual.getLevelAsFloat( 'eV' )
    endfMFList[8][MT]  = [ endfFormats.endfHeadLine( ZA, mass, LIS, LISO, NS, NO ) ]
    if( LMF != 3 ) : endfMFList[LMF][MT] = [ endfFormats.endfHeadLine( ZA, mass, LIS,    0, NS,  0 ) ]

    for reaction in MF8Channels :
        outputChannel = reaction.outputChannel
        if( len( outputChannel ) != 1 ) : raise Exception( 'Currently, production channel can only have one product; not %d' % len( outputChannel ) )
        product = outputChannel[0]

        Z, A, suffix, ZAP = product.particle.getZ_A_SuffixAndZA( )
        QI = outputChannel.getConstantQAs( 'eV', final = True )
        LFS2, level2 = product.particle.getLevelIndex(), product.particle.getLevelAsFloat( 'eV' )
        QM = QI + level2

        endfMFList[8][MT].append( endfFormats.endfHeadLine( ZAP, level2, LMF, LFS2, 0, 0 ) )

        if( LMF == 9 ) :        # Currenty, multiplicity_.toENDF6List only has one interpolation and its linear.
            multiplicity_ = product.multiplicity.getNativeData( )
            nRegions, nPoints, dataListAsStrings = multiplicity_.toENDF6List( )
            endfMFList[LMF][MT].append( endfFormats.endfHeadLine( QM, QI, ZAP, LFS2, nRegions, nPoints ) )
            endfMFList[LMF][MT] += dataListAsStrings
        elif( LMF == 10 ) :
            interpolationFlatData, flatData = reaction.crossSection.getNativeData( ).toENDF6Data( MT, endfMFList, tempInfo, level2 )
            endfMFList[LMF][MT].append( endfFormats.endfHeadLine( QM, QI, ZAP, LFS2, len( interpolationFlatData ) / 2, len( flatData ) / 2 ) )
            endfMFList[LMF][MT] += endfFormats.endfInterpolationList( interpolationFlatData )
            endfMFList[LMF][MT] += endfFormats.endfDataList( flatData )

    endfMFList[8][MT].append( endfFormats.endfSENDLineNumber( ) )
    if( LMF != 3 ) : endfMFList[LMF][MT].append( endfFormats.endfSENDLineNumber( ) )

def gammasToENDF6_MF12_13( MT, MF, endfMFList, flags, tempInfo, gammas ) :

    doMF4AsMF6 = tempInfo['doMF4AsMF6']
    tempInfo['doMF4AsMF6'] = False
    ZA, mass, MFGammas, NI, continuum, NK = tempInfo[ 'ZA' ], tempInfo[ 'mass' ], [], 0, None, len( gammas )
    total, piecewise = None, None
    for gammaIndex, gamma in enumerate( gammas ) :
        distributions = gamma.distributions
        component = distributions.components[distributions.nativeData]

        originationLevel = 0
        if( 'originationLevel' in gamma.attributes ) : 
            originationLevel = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( gamma.getAttribute( 'originationLevel' ) ).getValueAs( 'eV' )

        isPrimary = 'primary' in gamma.attributes
        LP, LF, levelEnergy = 0, 2, 0.
        if( isPrimary ) :
            gammaEnergy = gamma.getAttribute(  'primary' ).getValueAs( 'eV' )
            LP = 2
        elif( 'discrete' in gamma.attributes ) :
            gammaEnergy = gamma.getAttribute( 'discrete' ).getValueAs( 'eV' )
            if( originationLevel != 0 ) : LP = 1
        else :
            if( continuum is None ) :
                continuum = gamma
            else :
                raise Exception( 'Multiple continuum gammas detected for MT=%s' % ( MT ) )
            gammaEnergy, LP, LF = 0, 0, 1
            energyForm = component.energyComponent.forms[component.energyComponent.nativeData]
            NC, MF15 = energyForm.toENDF6( flags, tempInfo )
            endfMFList[15][MT] = [ endfFormats.endfContLine( ZA, mass, 0, 0, NC, 0 ) ]
            endfMFList[15][MT] += MF15
            endfMFList[15][MT].append( endfFormats.endfSENDLineNumber( ) )

        isNotIsotropic = 0
        if( LF == 2 ) :
            angularForm = component.forms[component.nativeData]
        else :
            angularForm = component.angularComponent.forms[component.angularComponent.nativeData]
        if( angularForm.moniker == distributionModule.base.isotropicFormToken ) :
            MF14 = None
            LI, LTT = 1, 0
        elif( angularForm.moniker == distributionModule.base.LegendrePointwiseFormToken ) :
            isNotIsotropic = 1
            dummy, dummy, frame, MF14 = angularForm.toENDF6( flags, tempInfo )
            del MF14[-1]
            MF14[0] = endfFormats.floatToFunky( gammaEnergy ) + endfFormats.floatToFunky( originationLevel ) + MF14[0][22:]
            LI, LTT = 0, 1
        else :
            raise Exception( 'Angular form = %s is not supported: MT=%s MF=%s' % ( angularForm.moniker, MT, MF ) )

        if( MF == 12 ) :
            LO = 1
            multiplicity_ = gamma.multiplicity.getFormByToken( gamma.multiplicity.getNativeDataToken( ) )
            if( isinstance( multiplicity_, multiplicity.piecewise ) ) :
                NR, NP, MF12or13Data = multiplicity_.toENDF6List( )
            else :
                NR, NP, MF12or13Data = multiplicity_.toENDF6List( )
            if( type( NR ) == type( [] ) ) :
                if( gammaEnergy == 0. ) :
                    for index, NRsub in enumerate( endfFormats.endfInterpolationList( NR ) ) : MF12or13Data.insert( index, NRsub )
                NR = len( NR ) / 2
            if( isinstance( multiplicity_, multiplicity.piecewise ) ) :
                if( piecewise is not None ) :
                    if len(piecewise) == len(multiplicity_) and all( [piecewise[i].getDomain() == multiplicity_[i].getDomain()
                            for i in range(len(piecewise))] ):
                        for i in range(len(piecewise)):
                            piecewise[i].setData( (piecewise[i] + multiplicity_[i]).copyDataToXYs() )
                    else:
                        raise Exception( 'MF=12 piecewise multiplicities must all have same energy boundaries' )
                else: piecewise = multiplicity_
            elif( total is None )  :
                total = multiplicity_
            else :
                total = total + multiplicity_
            MF12 = [ endfFormats.endfContLine( gammaEnergy, originationLevel, LP, LF, NR, NP ) ] + MF12or13Data
            MFGammas.append( [ isNotIsotropic, gammaEnergy, originationLevel, MF12, MF14 ] )
        elif( MF == 13 ) :
            LO = 0
            multiplicity_ = gamma.multiplicity.getNativeData( )
            if( isinstance( multiplicity_, multiplicity.constant ) ) :
                try :
                    crossSectionMF13 = tempInfo['crossSectionMF13']
                except :
                    crossSectionMF13 = None
                if( 'crossSectionMF13' is not None ) :
                    crossSection = crossSectionMF13[gammaIndex]
                else :
                    crossSection = tempInfo['crossSection']
                axes_ = crossSection.axes
                xsec_mult = [ [ x, y * multiplicity_.getValue( x ) ] for x, y in crossSection ]
            else :
                axes_ = multiplicity_.axes
                xsec_mult = [ [ x, y * tempInfo['crossSection'].getValue( x ) ] for x, y in multiplicity_ ]
            xsec_mult = XYs.XYs( axes_, xsec_mult, 1e-3 )
            if( total is None ) :
                total = xsec_mult
            else :
                total, xsec_mult_ = total.mutualify( 1e-6, 1e-6, True, xsec_mult, 1e-6, 1e-6, True )
                total = total + xsec_mult_
            NR = 1; NP = len( xsec_mult )
            MF13 = [ endfFormats.endfContLine( gammaEnergy, originationLevel, LP, LF, NR, NP ) ]
            MF13 += [ endfFormats.endfInterpolationLine( [ len( xsec_mult ), endfFormats.twoAxesToENDFInterpolation( axes_, 0 ) ] ) ]
            MF13 += endfFormats.endfNdDataList( xsec_mult )
            MFGammas.append( [ isNotIsotropic, gammaEnergy, originationLevel, MF13, MF14 ] )
        else:
            pass

    endfMFList[MF][MT] = [ endfFormats.endfContLine( ZA, mass, LO, 0, NK, 0 ) ]
    if( NK > 1 ) :      # Need to calculate total multiplicity. This section is a kludge, only works with lin-lin data?????
        if( piecewise is None ) :
            NR, NP, interpolationList = 1, len( total ), [ len( total ), 2 ]
            totalMultiplicityData = total
        elif( total is None ):  # all multiplicities are piecewise
            interpolationList, NP, MF12or13Data = piecewise.toENDF6List( )
            NR = len( interpolationList ) / 2
            totalMultiplicityData = []
            for region in piecewise: totalMultiplicityData.extend( region.copyDataToXYs() )
        else :          # This is not very robust, it is actually a kludge.
            axes_ = piecewise[0].axes.copy( standAlone = True )
            axes_[0].unit = 'eV'
            piecewise_ = multiplicity.piecewise( axes_ )
            t_xys = [ [ tx, ty ] for tx, ty in total.copyDataToXYs( xUnit = 'eV' ) ]
            t_y = t_xys[0][1]
            for region in piecewise :
                r_xys = []
                for rx, ry in region.copyDataToXYs( xUnit = 'eV' ) :
                    while( ( len( t_xys ) > 0 ) and ( rx > t_xys[0][0] ) ) : del t_xys[0]
                    if( len( t_xys ) > 0 ) :
                        if( t_xys[0][0] <= rx ) : t_y = t_xys[0][1]
                    r_xys.append( [ rx, ry + t_y ] )
                axes_ = region.axes.copy( standAlone = True )
                axes_[0].unit = 'eV'
                piecewise_[region.index] = XYs.XYs( axes_, r_xys, region.getAccuracy( ) )

            interpolationList, NP, MF12or13Data = piecewise_.toENDF6List( )
            NR = len( interpolationList ) / 2
            totalMultiplicityData = []
            for region in piecewise_ :
                data = region.copyDataToXYs( )
                if( len( totalMultiplicityData ) > 0 ) :
                    if( totalMultiplicityData[-1][0] == data[0][0] ) : data = data[1:]
                totalMultiplicityData += data
        endfMFList[MF][MT].append( endfFormats.endfContLine( 0, 0, 0, 0, NR, NP ) )
        endfMFList[MF][MT] += endfFormats.endfInterpolationList( interpolationList )
        endfMFList[MF][MT] += endfFormats.endfNdDataList( totalMultiplicityData )

    MFGammas12_13 = [ [ gammaEnergy, originationLevel, MF12_13 ] for isNotIsotropic, gammaEnergy, originationLevel, MF12_13, MF14 in MFGammas ]
    MFGammas12_13.sort( reverse = True )
    for gammaEnergy, originationLevel, MF12_13 in MFGammas12_13 : endfMFList[MF][MT] += MF12_13  # Store MF 12 or 13 before sorting.
    MFGammas.sort( reverse = True )
    for isNotIsotropic, gammaEnergy, originationLevel, MF12_13, MF14 in MFGammas :
        if( isNotIsotropic ) : break
        NI += 1
    if( LI == 1 ) : NI = 0
    endfMFList[MF][MT].append( endfFormats.endfSENDLineNumber( ) )

    # if any gammas are anisotropic, MF14 data must be supplied for all gammas:
    isotropic = [gamma for gamma in MFGammas if not gamma[0]]
    anisotropic = [gamma for gamma in MFGammas if gamma[0]]
    if anisotropic:
        LI, LTT, NI = 0, 1, len(isotropic)
        endfMFList[14][MT] = [ endfFormats.endfContLine( ZA, mass, LI, LTT, NK, NI ) ]
        for gamma in isotropic:
            endfMFList[14][MT].append( endfFormats.endfContLine( gamma[1], gamma[2], 0,0,0,0 ) )
        for gamma in anisotropic:
            endfMFList[14][MT] += gamma[-1]
        endfMFList[14][MT].append( endfFormats.endfSENDLineNumber( ) )
    else:
        endfMFList[14][MT] = [ endfFormats.endfContLine( ZA, mass, LI, LTT, NK, NI ),
                endfFormats.endfSENDLineNumber() ]

    tempInfo['doMF4AsMF6'] = doMF4AsMF6

def gndToEndfInterpolationFlag( xInterpolationStr, yInterpolationStr ) :

    endf = { axes.linearToken : 0, axes.logToken : 1 }
    if( yInterpolationStr == axes.flatToken ) : return( 1 )
    if( yInterpolationStr == axes.chargedParticleToken ) : return( 6 )
    return( endf[xInterpolationStr] + 2 * ( 1 + endf[yInterpolationStr] ))

def axesToEndfInterpolationFlag( axes ) :

    if( len( axes ) != 2 ) : raise Exception( "axes must only be for a 2 axes object, instance has %s axes" % len( axes ) )
    independent, dependent, dummy = axes[0].interpolation.getInterpolationTokens( )
    return( gndToEndfInterpolationFlag( independent, dependent ) )

def axisToEndfInterpolationFlag( axis ) :

    from fudge.core.math.xData import axes
    if( axis.interpolation is None ) : raise Exception( 'axis is dependent and must be independent: label = %s' % axis.getLabel( ) )
    independent, dependent, qualifier = axis.interpolation.getInterpolationTokens( )
    interpolation = gndToEndfInterpolationFlag( independent, dependent )
    if( qualifier is not None ) : 
        if( qualifier == axes.correspondingPointsToken ) :
            interpolation += 10
        elif( qualifier == axes.unitBaseToken ) :
            interpolation += 20
        else :
            raise Exception( 'Unsupported interpolation qualfier = %s' % qualifier )
    return( interpolation )
