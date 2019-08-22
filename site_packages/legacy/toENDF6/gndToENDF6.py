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

from PoPs import misc as miscPoPsModule
from PoPs import IDs as IDsPoPsModule
from PoPs.families import nuclearLevel as nuclearLevelModule

from fudge.particles import nuclear

from fudge.gnd import sums as sumsModule
from fudge.gnd.productData import multiplicity as multiplicityModule
from fudge.gnd.reactionData import crossSection as crossSectionModule
from fudge.gnd.productData.distributions import base as distributionBaseModule
from fudge.gnd.productData.distributions import energy as energyModule
from fudge.gnd.productData.distributions import angular as angularModule
from fudge.gnd.productData.distributions import uncorrelated as uncorrelatedModule
from fudge.gnd.productData.distributions import energyAngular as energyAngularModule
from fudge.gnd.productData.distributions import unspecified as unspecifiedModule

import endfFormats as endfFormatsModule

def angularPointwiseEnergy2ENDF6( self, targetInfo, EinFactor = 1 ) :

    interpolation = gndToENDFInterpolationFlag( self.interpolation )
    if len(self) > 201:
        accuracy, angularData = self.thinToNumberOfPoints( 201 )
        assert len(angularData) <= 201
        angularData = angularData.normalize()
        print("          WARNING: ENDF limited to 201 mu values. %d points removed for incident energy %g %s" %
                (len(self) - len(angularData), self.value, self.axes[-1].unit) )
    else:
        angularData = self
    if( targetInfo['doMF4AsMF6'] ) :
        interpolation += 10
        ENDFDataList = [ endfFormatsModule.endfContLine( 0, self.value * EinFactor, interpolation, 0,
            2 * len( angularData ), len( angularData ) ) ]
        ENDFDataList += endfFormatsModule.endfNdDataList( angularData )
    else :
        ENDFDataList = [ endfFormatsModule.endfContLine( 0, self.value * EinFactor, 0, 0, 1, len( angularData ) ) ]
        ENDFDataList += endfFormatsModule.endfInterpolationList( [ len( angularData ), interpolation ] )
        ENDFDataList += endfFormatsModule.endfNdDataList( angularData )
    return( ENDFDataList )

def angularLegendreEnergy2ENDF6( self, targetInfo, EinFactor = 1 ) :

    data = self
    if( 'doProductionGamma' in targetInfo ) : data = data.coefficients[1:]
    ENDFDataList = [ endfFormatsModule.endfContLine( 0, self.value * EinFactor, 0, 0, len( data ), 0 ) ]
    ENDFDataList += endfFormatsModule.endfNdDataList( data )
    return( ENDFDataList )

def gammaType( component ) :

    energySubform, angularSubform = None, None
    isPrimary, isDiscrete = False, False
    if( isinstance( component, uncorrelatedModule.form ) ) :
        energySubform = component.energySubform.data
        if( isinstance( energySubform, energyModule.primaryGamma ) ) :
            isPrimary = True
        elif( isinstance( energySubform, energyModule.discreteGamma ) ) :
            isDiscrete = True
        angularSubform = component.angularSubform.data
    else :
        raise 'hell - fix me'
    return( isPrimary, isDiscrete, energySubform, angularSubform )

def gammasToENDF6_MF6_oneGamma( MT, endfMFList, flags, targetInfo, gamma, LANG, LEP, frame ) :

    component = gamma.distribution[targetInfo['style']]
    isPrimary, isDiscrete, energySubform, angularSubform = gammaType( component )
    if( isPrimary or isDiscrete ) : return( False )

    targetInfo['multiplicity'] = gamma.multiplicity[targetInfo['style']]

    if( isinstance( component, uncorrelatedModule.form ) ) :
        angularSubform = component.angularSubform.data
        energySubform = component.energySubform.data
        if( not( isinstance( angularSubform, angularModule.isotropic ) ) ) :
            raise Exception( 'Unsupport angular form = "%s"' % angularSubform.moniker )
        if( not( isinstance( energySubform, energyModule.regions2d ) ) ) : energySubform = [ energySubform ]
        interpolationEIn = gndToENDF2PlusDInterpolationFlag( energySubform[0].interpolation, energySubform[0].interpolationQualifier )
        EInFactor = PQUModule.PQU( 1, energySubform[0].axes[2].unit ).getValueAs( 'eV' )
        EOutFactor = PQUModule.PQU( 1, energySubform[0].axes[1].unit ).getValueAs( 'eV' )

        NE, NR = 0, 1
        ENDFDataList = []
        for i1, region in enumerate( energySubform ) :
            interpolationEIn = gndToENDF2PlusDInterpolationFlag( region.interpolation, region.interpolationQualifier )
            for energyInData in region :
                NE += 1
                EIn = energyInData.value
                data = []
                if( not( isinstance( energyInData, energyModule.regions1d ) ) ) : energyInData = [ energyInData ]
                for region2 in energyInData :
                    for Ep, probability in region2 :
                        data.append( Ep * EOutFactor )
                        data.append( probability / EOutFactor )
                ENDFDataList += [ endfFormatsModule.endfContLine( 0, EIn, 0, 0, len( data ), len( data ) / 2 ) ]
                ENDFDataList += endfFormatsModule.endfDataList( data )
    else :
        raise Exception( 'Unsupport continuum gamma distribution = "%s"' % component.moniker )

    interpolations = [ NE, interpolationEIn ]
    ENDFDataList.insert( 0, endfFormatsModule.endfContLine( 0, 0, LANG, LEP, NR, NE ) )
    ENDFDataList.insert( 1, endfFormatsModule.endfInterpolationLine( interpolations ) )

    toENDF6_MF6( MT, endfMFList, flags, targetInfo, 1, frame, ENDFDataList )
    return( True )

def gammasToENDF6_MF6( MT, endfMFList, flags, targetInfo, gammas ) :
# FIXME - Still need to convert energies to eV.

    def checkOrder( nOrders, data, gammaEnergy, EIn ) :

        if( nOrders < 0 ) : nOrders = len( data )
        if( nOrders != len( data ) ) :
            raise Exception( "nOrders = %d != len( data ) = %d for gammaEnergy = %s, incident energy = %s" % \
                    ( nOrders, len( data ), gammaEnergy, EIn ) )

    def getTotalMultiplicityRegion( region, multiplicityRegions ) :

        for i2, totalRegion in enumerate( multiplicityRegions ) :
            if( totalRegion.domainMin <= region.domainMin ) :
                if( totalRegion.domainMax > region.domainMin ) : return( i2, totalRegion )
        for i2, totalRegion in enumerate( multiplicityRegions ) : print i2, totalRegion.domainMin, totalRegion.domainMax
        raise Exception( 'Could not find total region with domain %s, %s' % ( self.domainMin, self.domainMax ) )

    if( len( gammas ) == 1 ) :          # Check for LLNL legacy data
        distribution = gammas[0].distribution
        component = distribution[targetInfo['style']]
        if( not( isinstance( component, uncorrelatedModule.form ) ) ) :
            raise 'hell - fix me: was check for LLNL special data'
            component.toENDF6( MT, endfMFList, flags, targetInfo )
            return

    LANG, interpolationEIn, nOrders, discreteGammasVsEIn, continuum = 1, 2, -1, {}, None
    massRatio = targetInfo['mass'] / ( targetInfo['mass'] + 1. )

    LEP, discreteLEP = -2, -2
    frame = None
    for gamma in gammas :               # Determine LEP and frame
        component = gamma.distribution[targetInfo['style']]
        isPrimary, isDiscrete, energySubform, angularSubform = gammaType( component )
        if( isPrimary or isDiscrete ) :
            if( not( isinstance( angularSubform, ( angularModule.XYs2d, angularModule.isotropic ) ) ) ) : raise 'hell - fix me'
            if( discreteLEP < 0 ) :
                discreteLEP = 2
                if( isinstance( angularSubform, angularModule.XYs2d ) ) :
                    if( len( angularSubform.axes ) == 3 ) :
                        if( angularSubform.interpolation == standardsModule.interpolation.flatToken ) : discreteLEP = 1
        else :
            if( LEP < 0 ) :
                if( isinstance( component, energyAngularModule.form ) ) :
                    energy1 = component.energyAngularForm[0]
                    LEP = gndToENDF2PlusDInterpolationFlag( energy1.interpolation, energy1.interpolationQualifier )
                else :
                    EpForm = energySubform[0]
                    if( isinstance( EpForm, energyModule.regions1d ) ) : EpForm = EpForm[0]
                    LEP = gndToENDF2PlusDInterpolationFlag( EpForm.interpolation, standardsModule.interpolation.noneQualifierToken )
        if( frame is None ) : frame = component.productFrame
    if( LEP < 0 ) : LEP = discreteLEP
    if( frame is None ) : return

    if( len( gammas ) == 1 ) :
        if( gammasToENDF6_MF6_oneGamma( MT, endfMFList, flags, targetInfo, gammas[0], LANG, abs( LEP ), frame ) ) : return

    total = [ tmp for tmp in targetInfo['reactionSuite'].sums.multiplicities if tmp.ENDF_MT == MT ]
    if( len( total ) > 1 ) :
            raise Exception( "Multiple total gamma multiplicities for MT=%d" % MT )
    elif( len( total ) == 1 ) :
        total = total[0]
        if( isinstance( total, sumsModule.multiplicitySum ) ) :
            totalMultiplicity = total.multiplicity[targetInfo['style']].copy( )
        elif( isinstance( total, multiplicityModule.component ) ) :
            totalMultiplicity = total[targetInfo['style']].copy( )
        else :
            raise TypeError( 'Unsupport multiplitiy type = "%s"' % total.moniker )
    else :                              # Total multiplicity is missing in sums section, compute from parts.
        multiplicities = [ gamma.multiplicity[targetInfo['style']] for gamma in gammas ]
        totalMultiplicity = multiplicities.pop( 0 )
        for multiplicity in multiplicities : totalMultiplicity += multiplicity
    if( isinstance( totalMultiplicity, multiplicityModule.regions1d ) ) :
        totalMultiplicityRegions = totalMultiplicity.copy( )
    else :
        totalMultiplicityRegions = multiplicityModule.regions1d( axes = totalMultiplicity.axes )
        totalMultiplicityRegions.append( totalMultiplicity )

# This data was originally a Legendre expansions with only L=0 terms, was converted to uncorrelated.
# Also, original data stored one multiplicity for all gammas, with weights for individual gammas.
# Convert back to Legendre.

    discreteWeightsRegions = [ {} for region in totalMultiplicityRegions ]
    for gamma in gammas :
        multiplicity = gamma.multiplicity[targetInfo['style']]
        component = gamma.distribution[targetInfo['style']]
        isPrimary, isDiscrete, energySubform, angularSubform = gammaType( component )
        if( not( isPrimary or isDiscrete ) ) :
            continuum = gamma
            continue

        if( isinstance( multiplicity, multiplicityModule.regions1d ) ) :
            multiplicityRegions = multiplicity.copy( )
        else :
            if( isinstance( multiplicity, multiplicityModule.XYs1d ) ) :
                multiplicityRegions = multiplicityModule.regions1d( axes = multiplicity.axes )
                multiplicityRegions.append( multiplicity )
            else :
                raise Exception( 'Unsupported multiplicity "%s"' % multiplicity.moniker )
        for i1, region in enumerate( multiplicityRegions ) :
            i2, totalRegion = getTotalMultiplicityRegion( region, totalMultiplicityRegions )
            for EIn, probability  in region :
                total = totalRegion.evaluate( EIn )
                if( total != 0 ) : probability /= total
                gammaEnergy = PQUModule.PQU( energySubform.value, energySubform.axes[1].unit ).getValueAs( 'eV' )
                realGammaEnergy = gammaEnergy
                if( isPrimary ) : realGammaEnergy = -( gammaEnergy + massRatio * EIn )
                if( EIn not in discreteWeightsRegions[i2] ) : discreteWeightsRegions[i2][EIn] = []
                discreteWeightsRegions[i2][EIn].append( [ realGammaEnergy, probability ] )
    if( continuum is None ) :               # Check if we need to fix probability for lowest energy where multiplicity is 0.
        region = discreteWeightsRegions[0]
        EIn = sorted( region.keys( ) )[0]
        total = sum( [ probability for energy, probability in region[EIn] ] )
        if( total == 0 ) :
            norm = 1 / len( region[EIn] )
            for xy in region[EIn] : xy[1] = norm

    continuumGammasRegions = [ {} for region in totalMultiplicityRegions ]
    if( continuum is not None ) :
        multiplicity = continuum.multiplicity[targetInfo['style']]
        if( isinstance( multiplicity, multiplicityModule.regions1d ) ) :
            multiplicityRegions = multiplicity.copy( )
        else :
            if( isinstance( multiplicity, multiplicityModule.XYs1d ) ) :
                multiplicityRegions = multiplicityModule.regions1d( axes = multiplicity.axes )
                multiplicityRegions.append( multiplicity )
            else :
                raise Exception( 'Unsupported multiplicity "%s"' % multiplicity.moniker )

        component = continuum.distribution[targetInfo['style']]
        if( isinstance( component, uncorrelatedModule.form ) ) :
            angularSubform = component.angularSubform.data
            energySubform = component.energySubform.data
            if( not( isinstance( angularSubform, angularModule.isotropic ) ) ) :
                raise Exception( 'Unsupport angular form = "%s"' % angularSubform.moniker )
            if( not( isinstance( energySubform, energyModule.regions2d ) ) ) : energySubform = [ energySubform ]
            interpolationEIn = gndToENDF2PlusDInterpolationFlag( energySubform[0].interpolation, energySubform[0].interpolationQualifier )
            EInFactor = PQUModule.PQU( 1, energySubform[0].axes[2].unit ).getValueAs( 'eV' )
            EOutFactor = PQUModule.PQU( 1, energySubform[0].axes[1].unit ).getValueAs( 'eV' )
            for i1, region in enumerate( energySubform ) :
                i2, multiplicityRegion = getTotalMultiplicityRegion( region, multiplicityRegions )
                i2, totalRegion = getTotalMultiplicityRegion( region, totalMultiplicityRegions )
                for energyInData in region :
                    EIn = energyInData.value
                    total = totalRegion.evaluate( EIn )
                    continuumMultiplicity = multiplicityRegion.evaluate( EIn )
                    data = []
                    if( isinstance( energyInData, energyModule.regions1d ) ) :
                        for region2 in energyInData :
                            for Ep, probability in region2 :
                                if( total != 0 ) : probability *= continuumMultiplicity / total
                                data.append( Ep * EOutFactor )
                                data.append( probability )
                    else :
                        for Ep, probability in energyInData :
                            if( total != 0 ) : probability *= continuumMultiplicity / total
                            data.append( Ep * EOutFactor )
                            data.append( probability )
                    continuumGammasRegions[i2][EIn * EInFactor] = data

        elif( isinstance( component, energyAngularModule.form ) ) :
            raise 'hell - FIXME'
            if( len( continuumGammasRegions ) > 0 ) :
                raise Exception( 'Multiple regions not supported for %s' % component.moniker )
            form = component.energyAngularForm
            interpolationEIn = gndToENDF2PlusDInterpolationFlag( form.interpolation, form.interpolationQualifier )
            EInFactor = PQUModule.PQU( 1, component.axes[-1].unit ).getValueAs( 'eV' )
            EOutFactor = PQUModule.PQU( 1, component.axes[-2].unit ).getValueAs( 'eV' )
            incident_e_vals = []
            for energy_in in form :
                EIn = energy_in.value
                continuumGammas = []
                for EoutCl in energy_in :
                    nOrders = checkOrder( nOrders, EoutCl, gammaEnergy, EIn )
                    continuumGammas += [ EoutCl.value * EOutFactor ] + EoutCl.coefficients
                incident_e_vals.append( EIn * EInFactor )
                continuumGammasRegions[0][EIn * EInFactor ] = continuumGammas
        else :
            raise Exception( 'Unsupport continuum gamma distribution = "%s"' % component.moniker )

    NE = 0
    interpolations = []
    for i1, discreteWeightsRegion in enumerate( discreteWeightsRegions ) :
        continuumGammasRegion = continuumGammasRegions[i1]
        N = len( sorted( set( discreteWeightsRegion.keys( ) + continuumGammasRegion.keys( ) ) ) )
        NE += N
        interpolations += [ NE, interpolationEIn ]
    NR = len( totalMultiplicityRegions )

    ENDFDataList = [ endfFormatsModule.endfContLine( 0, 0, LANG, abs( LEP ), NR, NE ) ]
    ENDFDataList.append( endfFormatsModule.endfInterpolationLine( interpolations ) )

    nOrders = 1                         # FIXME
    wordsPerEout = nOrders + 1
    for i1, discreteWeightsRegion in enumerate( discreteWeightsRegions ) :
        continuumGammasRegion = continuumGammasRegions[i1]
        EIns = sorted( set( discreteWeightsRegion.keys( ) + continuumGammasRegion.keys( ) ) )
        for EIn in EIns :
            ND = 0
            gammaData = []
            if( EIn in discreteWeightsRegion ) :
                discreteGammas = sorted( discreteWeightsRegion[EIn], reverse = True )
                for discreteGamma in discreteGammas : gammaData += discreteGamma
                ND = len( gammaData ) / wordsPerEout
            if( EIn in continuumGammasRegion ) : gammaData += continuumGammasRegion[EIn]
            NW = len( gammaData )
            NEP = len( gammaData ) / wordsPerEout
            ENDFDataList += [ endfFormatsModule.endfContLine( 0, EIn, ND, nOrders - 1, NW, NEP ) ]
            ENDFDataList += endfFormatsModule.endfDataList( gammaData )
    targetInfo['multiplicity'] = totalMultiplicity
    toENDF6_MF6( MT, endfMFList, flags, targetInfo, 1, frame, ENDFDataList )

def toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, frame, MF6 ) :

    reactionSuite = targetInfo['reactionSuite']
    MF6or26 = { 3 : 6, 23 : 26 }[targetInfo['crossSectionMF']]
    targetInfo['MF6LCTs'].append( { standardsModule.frames.labToken : 1, standardsModule.frames.centerOfMassToken : 2 }[frame] )
    LIP = 0
    if( ( targetInfo['product'].id in targetInfo['metastables'] ) and not( ( 50 <= MT <= 91 ) or ( 600 <= MT <= 850 ) ) ) :
            # set LIP to metastable index of product UNLESS excited state of product is encoded in MT number:
        alias = targetInfo['metastables'][ targetInfo['product'].id ]
        LIP = int( alias.metaStableIndex )
    if( MT not in endfMFList[MF6or26] ) : endfMFList[MF6or26][MT] = []
    particleID = targetInfo['zapID']    # get ZAP for the outgoing particle
    if( particleID == IDsPoPsModule.electron ) :
        ZAP = 11
    else :
        ZAP = miscPoPsModule.ZA( reactionSuite.PoPs[particleID] )
    productMass = targetInfo['particleMass']
    if( ( particleID == IDsPoPsModule.neutron ) and ( MT in [ 18 ] ) ) :                            # Special case when fission has MF = 6 data.
        nPoints = 2
        interpolationFlatData = [ nPoints, 2 ]
        multiplicityList = endfFormatsModule.endfDataList( [ [ targetInfo['EMin'], 1. ], [ targetInfo['EMax'], 1. ] ] )
    else :
        multiplicity = targetInfo['multiplicity']
        if( isinstance( multiplicity, multiplicityModule.component ) ) : multiplicity = multiplicity[targetInfo['style']]
        interpolationFlatData, nPoints, multiplicityList = multiplicity.toENDF6List( targetInfo )
    if( ( particleID == IDsPoPsModule.photon ) and ( LAW == 2 ) ) :
        projectileMass = reactionSuite.PoPs[reactionSuite.projectile].getMass( 'eV/c**2' )
        targetMass = reactionSuite.PoPs[reactionSuite.projectile].getMass( 'eV/c**2' )
        AWP = targetInfo['Q']
        AWP /= 1. + 0.5 * AWP / ( projectileMass + targetMass )             # Should be divided by residual mass, but this is close.
    else :
        AWP = productMass
    ENDFHeaderList = [ endfFormatsModule.endfContLine( ZAP, AWP, LIP, LAW, len( interpolationFlatData ) / 2, nPoints ) ]
    ENDFHeaderList += endfFormatsModule.endfInterpolationList( interpolationFlatData )
    ENDFHeaderList += multiplicityList
    endfMFList[MF6or26][MT] += ENDFHeaderList + MF6

def upDateENDFMF8Data( endfMFList, targetInfo ) :

    for MT in targetInfo['MF8'] : upDateENDF_MT_MF8Data( MT, endfMFList, targetInfo )

def upDateENDF_MT_MF8Data( MT, endfMFList, targetInfo ) :

    reactionSuite = targetInfo['reactionSuite']
    MF8Channels = targetInfo['MF8'][MT]
    ZA, mass = targetInfo[ 'ZA' ], targetInfo[ 'mass' ]
    LIS, LISO, NO, NS = targetInfo['LIS'], targetInfo['LISO'], 1, len( MF8Channels )
    firstReaction = MF8Channels[0]
    residual = firstReaction.outputChannel[0]

    crossSection_ = firstReaction.crossSection[targetInfo['style']]  # LMF determines which MF section to write to.
    if( isinstance( crossSection_, crossSectionModule.reference ) ) : 
        multiplicity = residual.multiplicity[targetInfo['style']]
        if(   isinstance( multiplicity, multiplicityModule.constant1d ) ) :
            LMF = 3
        elif( isinstance( multiplicity, multiplicityModule.reference ) ) :
            LMF = 6
        else :
            LMF = 9
    else :
        LMF = 10

    level = 0
    if( hasattr( residual, 'getLevelAsFloat' ) ) : level = residual.getLevelAsFloat( 'eV' )
    endfMFList[8][MT]  = [ endfFormatsModule.endfHeadLine( ZA, mass, LIS, LISO, NS, NO ) ]
    if( LMF not in [3,6] ) : endfMFList[LMF][MT] = [ endfFormatsModule.endfHeadLine( ZA, mass, LIS,    0, NS,  0 ) ]

    for reaction in MF8Channels :
        outputChannel = reaction.outputChannel
        if( len( outputChannel ) != 1 ) : raise Exception( 'Currently, production channel can only have one product; not %d' % len( outputChannel ) )
        product = outputChannel[0]

        particle = reactionSuite.PoPs[product.id]
        ZAP = miscPoPsModule.ZA( particle )

        QI = outputChannel.getConstantQAs( 'eV', final = True )
        LFS2, level2 = 0, 0
        if( isinstance( particle, nuclearLevelModule.particle ) ) :
            LFS2 = particle.intIndex
            level2 = particle.energy[0].float( 'eV' )
        QM = QI + level2

        endfMFList[8][MT].append( endfFormatsModule.endfHeadLine( ZAP, level2, LMF, LFS2, 0, 0 ) )

        if( LMF == 9 ) :        # Currenty, multiplicity.toENDF6List only has one interpolation and its linear.
            multiplicity = product.multiplicity[targetInfo['style']]
            interpolationFlatData, nPoints, multiplicityList = multiplicity.toENDF6List( targetInfo )
            endfMFList[LMF][MT].append( endfFormatsModule.endfHeadLine( QM, QI, ZAP, LFS2, len( interpolationFlatData ) / 2, nPoints ) )
            endfMFList[LMF][MT] += endfFormatsModule.endfInterpolationList( interpolationFlatData )
            endfMFList[LMF][MT] += multiplicityList
        elif( LMF == 10 ) :
            interpolationFlatData, flatData = reaction.crossSection[targetInfo['style']].toENDF6Data( MT, endfMFList, targetInfo, level2 )
            endfMFList[LMF][MT].append( endfFormatsModule.endfHeadLine( QM, QI, ZAP, LFS2, len( interpolationFlatData ) / 2, len( flatData ) / 2 ) )
            endfMFList[LMF][MT] += endfFormatsModule.endfInterpolationList( interpolationFlatData )
            endfMFList[LMF][MT] += endfFormatsModule.endfDataList( flatData )

    endfMFList[8][MT].append( endfFormatsModule.endfSENDLineNumber( ) )
    if( LMF not in [3,6] ) : endfMFList[LMF][MT].append( endfFormatsModule.endfSENDLineNumber( ) )

def adjustMF13Multiplicity( interpolationInfo, xsec_mult, multiplicity, crossSection ) :

        xsec_mult2 = []
        for energyIn, multiplicityValue in multiplicity :
            xsec_mult2.append( energyIn )
            if( crossSection is not None ) : multiplicityValue *= crossSection.evaluate( energyIn )
            xsec_mult2.append( multiplicityValue )
        if( len( xsec_mult ) > 0 ) :
            if( xsec_mult[-2:] == xsec_mult2[:2] ) : xsec_mult2 = xsec_mult2[2:]
        xsec_mult += xsec_mult2
        interpolationInfo += [ len( xsec_mult ) / 2, gndToENDFInterpolationFlag( multiplicity.interpolation ) ]

def gammasToENDF6_MF12_13( MT, MF, endfMFList, flags, targetInfo, gammas ) :
    """
    MF=12 stores the multiplicity for gamma products, can be converted directly to/from GND.
    MF=13 stores 'gamma production cross section', equal to cross section * multiplicity.
    ENDF-to-GND translator converts to multiplicity (dividing by reaction cross section).
    When writing back to ENDF-6 MF=13, need to convert back using adjustMF13Multiplicity.
    """

    doMF4AsMF6 = targetInfo['doMF4AsMF6']
    targetInfo['doMF4AsMF6'] = False
    ZA, mass, MFGammas, NI, continuum, NK, NKp = targetInfo[ 'ZA' ], targetInfo[ 'mass' ], [], 0, None, len( gammas ), 0

    total, piecewise, recomputeTotal = None, None, False

    crossSection = None
    if( MF == 13 ) :
        crossSection = targetInfo['crossSection']
        crossSection = crossSection.toPointwise_withLinearXYs( accuracy = 1e-3, upperEps = 1e-8 )

    if( NK > 1 ) :
        # If more than one gamma is present, both MF=12 and MF=13 start with the sum over all gammas.
        # Search for the correct multiplicitySum in reactionSuite/sums section

        total = [tmp for tmp in targetInfo['reactionSuite'].sums.multiplicities if tmp.ENDF_MT == MT]
        if len(total) > 1:
            raise Exception("Cannot find unique gamma multiplicity sum for MT%d" % MT)
        elif not total:     # total multiplicity is missing from sums section, need to recompute from parts
            recomputeTotal = True
        else:
            total = total[0].multiplicity.evaluated

    for gammaIndex, gamma in enumerate( gammas ) :
        distribution = gamma.distribution
        component = distribution[targetInfo['style']]
        if( isinstance( component, unspecifiedModule.form ) ) : continue
        NKp += 1

        originationLevel = 0
        if( 'originationLevel' in gamma.attributes ) : 
            originationLevel = PQUModule.PQU( gamma.getAttribute( 'originationLevel' ) ).getValueAs( 'eV' )

        isPrimary, isDiscrete = False, False
        if( isinstance( component, uncorrelatedModule.form ) ) :
            energySubform = component.energySubform.data
            if( isinstance( energySubform, energyModule.primaryGamma ) ) :
                isPrimary = True
            elif( isinstance( energySubform, energyModule.discreteGamma ) ) :
                isDiscrete = True
            angularSubform = component.angularSubform.data
        else :
            raise 'hell - fix me'
        LP, LF, levelEnergy = 0, 2, 0.
        if( isPrimary ) :
            gammaEnergy = PQUModule.PQU( energySubform.value, energySubform.axes[1].unit ).getValueAs( 'eV' )
            LP = 2
        elif( isDiscrete ) :
            gammaEnergy = PQUModule.PQU( energySubform.value, energySubform.axes[1].unit ).getValueAs( 'eV' )
            if( originationLevel != 0 ) : LP = 1
        else :
            if( continuum is not None ) : raise Exception( 'Multiple continuum gammas detected for MT=%s' % ( MT ) )
            continuum = gamma
            energySubform = component.energySubform.data
            gammaEnergy, LP, LF = 0, 0, 1
            NC, MF15 = energySubform.toENDF6( flags, targetInfo )
            endfMFList[15][MT] = [ endfFormatsModule.endfContLine( ZA, mass, 0, 0, NC, 0 ) ]
            endfMFList[15][MT] += MF15
            endfMFList[15][MT].append( endfFormatsModule.endfSENDLineNumber( ) )

        isNotIsotropic = 0
        if( isinstance( angularSubform, angularModule.isotropic ) ) :
            MF14 = None
            LI, LTT = 1, 0
        elif( isinstance( angularSubform, angularModule.XYs2d ) ) :
            isNotIsotropic = 1
            targetInfo['doProductionGamma'] = True
            dummy, dummy, MF14 = angularSubform.toENDF6( flags, targetInfo )
            del targetInfo['doProductionGamma']
            del MF14[-1]
            MF14[0] = endfFormatsModule.floatToFunky( gammaEnergy ) + endfFormatsModule.floatToFunky( originationLevel ) + MF14[0][22:]
            LI, LTT = 0, 1
        else :
            raise 'hell - fix me'

        multiplicity = gamma.multiplicity.evaluated

        if recomputeTotal:
            if( isinstance( multiplicity, multiplicityModule.regions1d ) ) :
                if( piecewise is not None ) :
                    if( len( piecewise ) != len( multiplicity ) ) :
                        raise Exception( 'MF=12/13 piecewise multiplicities must all have same number of regions' )
                    for i1, region in enumerate( piecewise ) :
                        region.setData( ( region + multiplicity[i1] ).copyDataToXYs( ) )
                else :
                    piecewise = multiplicity
                total = piecewise[0]            # BRB, FIXME, this is a kludge until total is put into sums.
            elif( total is None )  :
                total = multiplicity
            else :
                total = total + multiplicity

        if( MF == 12 ) :
            LO = 1
            NR, NP, MF12or13Data = multiplicity.toENDF6List( targetInfo )
            if( type( NR ) == type( [] ) ) :
                for index, NRsub in enumerate( endfFormatsModule.endfInterpolationList( NR ) ) : MF12or13Data.insert( index, NRsub )
                NR = len( NR ) / 2
            MF12 = [ endfFormatsModule.endfContLine( gammaEnergy, originationLevel, LP, LF, NR, NP ) ] + MF12or13Data
            MFGammas.append( [ isNotIsotropic, gammaEnergy, originationLevel, MF12, MF14 ] )

        elif( MF == 13 ) :
            LO = 0

            interpolationInfo = []
            xsec_mult = []
            if( isinstance( multiplicity, multiplicityModule.XYs1d ) ) :
                adjustMF13Multiplicity( interpolationInfo, xsec_mult, multiplicity, crossSection )
            elif( isinstance( multiplicity, multiplicityModule.regions1d ) ) :
                for region in multiplicity : adjustMF13Multiplicity( interpolationInfo, xsec_mult, region, crossSection )
            else :
                raise Exception( 'Unsupport multiplicity "%s" for MF=13' % multiplicity.moniker )

            NR = len( interpolationInfo ) / 2
            NP = len( xsec_mult ) / 2
            MF13 = [ endfFormatsModule.endfContLine( gammaEnergy, originationLevel, LP, LF, NR, NP ) ]
            MF13 += [ endfFormatsModule.endfInterpolationLine( interpolationInfo ) ]
            MF13 += endfFormatsModule.endfNdDataList( xsec_mult )
            MFGammas.append( [ isNotIsotropic, gammaEnergy, originationLevel, MF13, MF14 ] )
        else:
            pass
    if( NKp == 0 ) : return

    endfMFList[MF][MT] = [ endfFormatsModule.endfContLine( ZA, mass, LO, 0, NK, 0 ) ]
    if( NK > 1 ) :          # If more than one gamma is present, both MF=12 and MF=13 start with the sum over all gammas.

        totalInterpolationInfo = []
        xsec_mult = []
        if( isinstance( total, multiplicityModule.XYs1d ) ) :
            adjustMF13Multiplicity( totalInterpolationInfo, xsec_mult, total, crossSection )
        elif( isinstance( total, multiplicityModule.regions1d ) ) :
            for region in total : adjustMF13Multiplicity( totalInterpolationInfo, xsec_mult, region, crossSection )
        else:
            raise Exception( 'Unsupported total multiplicity type "%s" for MF=13' % total.moniker )

        NR = len( totalInterpolationInfo ) / 2
        NP = len( xsec_mult ) / 2
        endfMFList[MF][MT].append( endfFormatsModule.endfContLine( 0, 0, 0, 0, NR, NP ) )
        endfMFList[MF][MT] += [ endfFormatsModule.endfInterpolationLine( totalInterpolationInfo ) ]
        endfMFList[MF][MT] += endfFormatsModule.endfNdDataList( xsec_mult )

    MFGammas12_13 = [ [ gammaEnergy, originationLevel, MF12_13 ] for isNotIsotropic, gammaEnergy, originationLevel, MF12_13, MF14 in MFGammas ]
    MFGammas12_13.sort( reverse = True )
    for gammaEnergy, originationLevel, MF12_13 in MFGammas12_13 : endfMFList[MF][MT] += MF12_13  # Store MF 12 or 13 before sorting.
    MFGammas.sort( reverse = True )
    for isNotIsotropic, gammaEnergy, originationLevel, MF12_13, MF14 in MFGammas :
        if( isNotIsotropic ) : break
        NI += 1
    if( LI == 1 ) : NI = 0
    endfMFList[MF][MT].append( endfFormatsModule.endfSENDLineNumber( ) )

    # if any gammas are anisotropic, MF14 data must be supplied for all gammas:
    isotropic = [gamma for gamma in MFGammas if not gamma[0]]
    anisotropic = [gamma for gamma in MFGammas if gamma[0]]
    if anisotropic:
        LI, LTT, NI = 0, 1, len(isotropic)
        endfMFList[14][MT] = [ endfFormatsModule.endfContLine( ZA, mass, LI, LTT, NK, NI ) ]
        for gamma in isotropic:
            endfMFList[14][MT].append( endfFormatsModule.endfContLine( gamma[1], gamma[2], 0,0,0,0 ) )
        for gamma in anisotropic:
            endfMFList[14][MT] += gamma[-1]
        endfMFList[14][MT].append( endfFormatsModule.endfSENDLineNumber( ) )
    else:
        endfMFList[14][MT] = [ endfFormatsModule.endfContLine( ZA, mass, LI, LTT, NK, NI ),
                endfFormatsModule.endfSENDLineNumber() ]

    targetInfo['doMF4AsMF6'] = doMF4AsMF6

def gndToENDFInterpolationFlag( interpolation ) :

    if( interpolation == standardsModule.interpolation.flatToken ) : return( 1 )
    if( interpolation == standardsModule.interpolation.linlinToken ) : return( 2 )
    if( interpolation == standardsModule.interpolation.linlogToken ) : return( 3 )
    if( interpolation == standardsModule.interpolation.loglinToken ) : return( 4 )
    if( interpolation == standardsModule.interpolation.loglogToken ) : return( 5 )
    if( interpolation == standardsModule.interpolation.chargedParticleToken ) : return( 6 )
    raise ValueError( 'unsupported interpolation = "%s"' % interpolation )

def gndToENDF2PlusDInterpolationFlag( interpolation, interpolationQualifier ) :

    interpolationFlag = gndToENDFInterpolationFlag( interpolation )
    if( interpolationQualifier == standardsModule.interpolation.noneQualifierToken ) :
        offset = 0
    elif( interpolationQualifier == standardsModule.interpolation.correspondingPointsToken ) :
        offset = 10
    elif( interpolationQualifier == standardsModule.interpolation.unitBaseToken ) :
        offset = 20
    else :
        raise ValueError( 'unsupport interpolationQualifier = "%s"' % str( interpolationQualifier )[:40] )
    return( interpolationFlag + offset )
