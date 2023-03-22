# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import enums as xDataEnumsModule
from xData import XYs1d as XYs1dModule
from xData import regions as regionsModule

from PoPs import specialNuclearParticleID as specialNuclearParticleIDPoPsModule
from PoPs import IDs as IDsPoPsModule
from PoPs.families import nuclide as nuclideModule
from PoPs.chemicalElements import misc as chemicalElementMiscModule

from brownies.legacy.converting import endf_endl as endf_endlModule

from fudge import enums as enumsModule
from fudge import outputChannel as outputChannelModule
from fudge.reactions import base as baseModule
from fudge.productData import multiplicity as multiplicityModule
from fudge.reactionData.doubleDifferentialCrossSection import thermalNeutronScatteringLaw
from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import CoulombPlusNuclearElastic
from fudge.productData.distributions import unspecified as unspecifiedModule
from fudge.productData.distributions import angular as angularModule
from fudge.productData.distributions import branching3d as branching3dModule
from fudge.productData.distributions import energyAngular as energyAngularModule

from .. import gndsToENDF6 as gndsToENDF6Module
from .. import endfFormats as endfFormatsModule

def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ) :

    def checkDecayProducts( parent ) :

        if( parent.outputChannel is None ) : return( False )
        doMF4AsMF6 = False
        for product in parent.outputChannel :
            conversionFlag = targetInfo['ENDFconversionFlags'].get(product,"")
            if( len( product.distribution ) and not isinstance( product.distribution.evaluated, unspecifiedModule.Form ) ) :
                if ( product.pid not in ( IDsPoPsModule.neutron, IDsPoPsModule.photon ) ) :   # MF4/5 can only handle neutrons
                    doMF4AsMF6 = True
                elif ( conversionFlag == 'MF6' ) :  # sometimes neutrons still go in MF6
                    doMF4AsMF6 = True
        return( doMF4AsMF6 )

    def divideIgnoring0DividedBy0( self, other ) :

        if( isinstance( self, XYs1dModule.XYs1d ) and isinstance( other, XYs1dModule.XYs1d ) ) :
            d = self.union( other.domainSlice( domainMin = self.domainMin, domainMax = self.domainMax ) )
            result = []
            for p in d :
                vp = other.evaluate( p[0] )
                if( vp == 0 ) :
                    if( p[1] != 0 ) : raise Exception( 'Divide non-zero number by zero at %e' % p[0] )
                else :
                    p[1] = p[1] / vp
                result.append( [ p[0], p[1] ] )
            return( XYs1dModule.pointwiseXY( data = result ) )

        elif ( isinstance( self, regionsModule.Regions1d ) and isinstance( other, regionsModule.Regions1d )
                and len(self)==len(other) ):
            regions = regionsModule.Regions1d( )
            for idx in range(len(self)):
                regions.append( XYs1dModule.XYs1d( data = divideIgnoring0DividedBy0( self[idx], other[idx] ) ) )
            return regions

        else:
            raise NotImplementedError( 'Divide XYs1d by regions or vice-versa' )

    def thinWeights( weights ) :

        i, n, thinnedWeights = 1, len( weights ), [ weights[0] ]
        if( n == 2 ) :
            thinnedWeights.append( weights[-1] )
        else :
            while( i < n ) :
                y = thinnedWeights[-1][1]
                while( True ) :
                    i += 1
                    if( i == n ) : break
                    if( abs( y - weights[i][1] ) > 1e-8 * y ) : break
                    if( abs( y - weights[i-1][1] ) > 1e-8 * y ) : break
                thinnedWeights.append( weights[i-1] )
        return( thinnedWeights )

    targetInfo['reaction'] = self

    if hasattr(self, 'doubleDifferentialCrossSection'):
        if self.doubleDifferentialCrossSection:
            differentialForm = gndsToENDF6Module.getForm( targetInfo['style'], self.doubleDifferentialCrossSection )
            if isinstance( differentialForm, CoulombPlusNuclearElastic.Form ) and differentialForm.data is None:
                differentialForm.RutherfordScattering.toENDF6(endfMFList, flags, targetInfo, verbosityIndent)
                return  # pure Rutherford scattering
            elif isinstance( differentialForm, (thermalNeutronScatteringLaw.coherentElastic.Form,
                                                thermalNeutronScatteringLaw.incoherentElastic.Form,
                                                thermalNeutronScatteringLaw.incoherentInelastic.Form) ):
                endfMFList.setdefault(7, {})
                differentialForm.toENDF6(endfMFList, flags, targetInfo, verbosityIndent)
                return  # ignore cross sections, multiplicities, etc.

    reactionSuite = targetInfo['reactionSuite']
    outputChannel = self.outputChannel
    if outputChannel.genre == enumsModule.Genre.production:
        print( '        toENDF6 does not support writing of "%s" channel' % outputChannel.genre )
        return
    MT = self.ENDF_MT
    if( flags['verbosity'] >= 10 ) : print('%s%s' % (verbosityIndent, self.label))
    targetInfo['Q'] = self.getQ( 'eV', final = False )
    conversionFlags = targetInfo['ENDFconversionFlags'].get( self.outputChannel.Q, "" )
    if( conversionFlags != "" ) : targetInfo['Q'] = float( conversionFlags.split( '=' )[1] )
    if( MT in [ 515, 517 ] ) : targetInfo['Q'] = 0.0
    targetInfo['QM'] = None

    LR, tryLR = 0, False
    if outputChannel.genre == enumsModule.Genre.twoBody:
        tryLR = True
    elif outputChannel.genre == enumsModule.Genre.NBody:
        if MT == 91: tryLR = True
    if tryLR and len(outputChannel) > 1:
        secondProduct = outputChannel[1]
        primaryResidualName, decayProducts, decayChannel, = secondProduct.pid.split( '_' )[0], [], secondProduct.outputChannel
        numberOfDistributions = 0
        if decayChannel is not None:
            for decayProduct in decayChannel :
                if not isinstance(decayProduct.distribution[0], unspecifiedModule.Form): numberOfDistributions += 1
                decayProductName = decayProduct.pid
                if decayProductName not in [primaryResidualName, IDsPoPsModule.photon]:
                    decayProducts.append(decayProductName)
        if len(decayProducts) == 1:   # Kludge for Carbon breakup into 3 alphas.
            if (primaryResidualName in ('C0', 'C12')) and (decayProducts == ['He4']): LR = 23
            else: LR = 1   # FIXME may want to use other more specific LR flags
        elif len(decayProducts) > 1:                                        # This must be a breakup reaction.
            if numberOfDistributions > 0:
                if numberOfDistributions != len(decayProducts):
                    raise Exception('Do not know what to do here: MT="%s' % MT)
                LR = 1
            else:
                MTProducts = endf_endlModule.ENDF_MTtoC_ProductList(0, '')
                for decayProduct in decayProducts[:-1]:
                    decayProductID = specialNuclearParticleIDPoPsModule.specialNuclearParticleID(
                        decayProduct, specialNuclearParticleIDPoPsModule.Mode.nuclide)
                    MTProducts.productCounts[decayProductID] += 1
                for MT_LR in [22, 23, 24, 25, 28, 29, 30, 32, 33, 34, 35, 36]:   # 39 and 40 not allowed in ENDF6
                    if endf_endlModule.endfLRtoC_ProductLists[MT_LR].productCounts == MTProducts.productCounts:
                        LR = MT_LR
                        break
                if (LR == 32) and (primaryResidualName == 'B10') and (decayProducts[-1] == 'He4'):
                    LR = 35   # Kludge for bad data.
        if LR != 0:
            QM = gndsToENDF6Module.getForm(targetInfo['style'], outputChannel.Q).evaluate(0) + \
                 gndsToENDF6Module.getForm(targetInfo['style'], decayChannel.Q).evaluate(0)
            targetInfo['QM'] = QM
    targetInfo['LRs'][MT] = LR

    level = 0.
    for product in outputChannel :
        baseID, chemicalElementSymbol, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( product.pid, qualifierAllowed = True )
        if( A is None ) : continue
        productID = product.pid.split( '{' )[0]
        _product = reactionSuite.PoPs[productID]
        if( isinstance( _product, nuclideModule.Particle ) ) :
            if( _product.index > 0 ) : level = _product.energy[0].float( 'eV' )

    targetInfo['EMin'], targetInfo['EMax'] = self.crossSection.domainMin, self.crossSection.domainMax     # Need to convert to eV.
    self.crossSection.toENDF6( MT, endfMFList, targetInfo, level, LR )

    if hasattr(self, 'doubleDifferentialCrossSection'):
        self.doubleDifferentialCrossSection.toENDF6(MT, endfMFList, targetInfo)

    doMF4AsMF6 = False
    if reactionSuite.projectile not in (IDsPoPsModule.neutron, IDsPoPsModule.photon): doMF4AsMF6 = True
    for product in outputChannel:
        if len(product.distribution) and not isinstance(product.distribution.evaluated, unspecifiedModule.Form):
            conversionFlag = targetInfo['ENDFconversionFlags'].get(product, "")
            if product.pid not in (IDsPoPsModule.neutron, IDsPoPsModule.photon) and conversionFlag != 'implicitProduct':
                if self.outputChannel.genre != enumsModule.Genre.twoBody:
                    doMF4AsMF6 = True
        doMF4AsMF6 = doMF4AsMF6 or checkDecayProducts( product )
    if (not doMF4AsMF6 and
            self.outputChannel.genre == enumsModule.Genre.twoBody and
            self.outputChannel[0].pid == IDsPoPsModule.neutron and
            isinstance(self.outputChannel[1].distribution.evaluated, angularModule.TwoBody) and
            isinstance(self.outputChannel[1].distribution.evaluated.angularSubform, angularModule.Recoil)):
        recoil = self.outputChannel[1]
        if recoil not in targetInfo['ENDFconversionFlags']:
            # make it implicit unless any decay product has explicit distribution
            if recoil.outputChannel is None or all([isinstance(product.distribution.evaluated,
                        (unspecifiedModule.Form, branching3dModule.Form)) for product in recoil.outputChannel]):
                targetInfo['ENDFconversionFlags'][product] = 'implicitProduct'
                if recoil.outputChannel is not None:
                    for recoilProduct in recoil.outputChannel:
                        targetInfo['ENDFconversionFlags'][recoilProduct] = 'implicitProduct'

    targetInfo['isFission'] = False
    targetInfo['totalFission'] = False
    if self.isFission():
        targetInfo['isFission'] = True
        targetInfo['totalFission'] = len(outputChannel.fissionFragmentData.delayedNeutrons) == 0 and reactionSuite.projectile != IDsPoPsModule.photon
        doMF4AsMF6 = False
        targetInfo['doMF4AsMF6'] = doMF4AsMF6
        if( len( outputChannel.fissionFragmentData.fissionEnergyReleases ) > 0 ) :
            outputChannel.fissionFragmentData.fissionEnergyReleases.evaluated.toENDF6( endfMFList, flags, targetInfo )

        delayedNubar = targetInfo.get( 'totalDelayedNubar', 0 )
        computeDelayed = ( delayedNubar == 0 )
        for productIndex, delayedNeutron in enumerate( outputChannel.fissionFragmentData.delayedNeutrons ) :
            if computeDelayed: delayedNubar += gndsToENDF6Module.getForm( targetInfo['style'], delayedNeutron.product.multiplicity )
            try:
                rate = delayedNeutron.rate[targetInfo['style']]
            except:
                rate = delayedNeutron.rate[0]
            targetInfo['delayedRates'].append(float( rate.pqu( unit = '1/s' ) ))
        for delayedNeutron in outputChannel.fissionFragmentData.delayedNeutrons :
            multiplicity = gndsToENDF6Module.getForm(targetInfo['style'], delayedNeutron.product.multiplicity)
            if( isinstance( multiplicity, multiplicityModule.Unspecified ) ) : continue
            weight = divideIgnoring0DividedBy0( multiplicity, delayedNubar )
            delayedNeutron.product.ENDF6_delayedNubarWeights = thinWeights( weight )
            delayedNeutron.product.toENDF6( 455, endfMFList, flags, targetInfo, verbosityIndent = verbosityIndent + '    ' )

    if MT in [ 526, 527 ]: doMF4AsMF6 = True
    targetInfo['doMF4AsMF6'] = doMF4AsMF6
    targetInfo['MF6LCTs'], targetInfo['gammas'] = [], []
    numberOfGammas = len([product for product in outputChannel.products if product.pid == IDsPoPsModule.photon])
    for productIndex, product in enumerate(outputChannel.products):
        if product.pid == IDsPoPsModule.photon and outputChannel.genre != enumsModule.Genre.twoBody and MT != 527:
            if numberOfGammas != 1 or not isinstance(product.distribution[0], energyAngularModule.Form):
                targetInfo['gammas'].append(product)
                continue
        targetInfo['productIndex'] = str(productIndex)
        targetInfo['productToken'] = product.pid
        targetInfo['productLabel'] = product.label
        product.toENDF6(MT, endfMFList, flags, targetInfo, verbosityIndent = verbosityIndent + '    ')

    gammas = []
    for gamma in targetInfo['gammas']:
        multiplicity = gndsToENDF6Module.getForm(targetInfo['style'], gamma.multiplicity)
        if isinstance(multiplicity, multiplicityModule.Branching1d): continue
        gammas.append(gamma)
    if len(gammas):
        gamma = gammas[0]
        targetInfo['zapID'] = gamma.pid
        targetInfo['particleMass'] = gamma.getMass('eV/c**2')
        multiplicity = gamma.multiplicity
        conversionFlag = targetInfo['ENDFconversionFlags'].get(gamma,"")
        if conversionFlag == 'MF6':
            targetInfo['multiplicity'] = gamma.multiplicity
            gndsToENDF6Module.gammasToENDF6_MF6( MT, endfMFList, flags, targetInfo, gammas )
        elif 'MF13' in conversionFlag:
            targetInfo['crossSection'] = gndsToENDF6Module.getForm( targetInfo['style'], self.crossSection )
            gndsToENDF6Module.gammasToENDF6_MF12_13( MT, 13, endfMFList, flags, targetInfo, gammas )
        elif isinstance(multiplicity, multiplicityModule.Constant1d):
            pass
        elif isinstance(multiplicity, multiplicityModule.Unspecified):
            pass
        else:
            gndsToENDF6Module.gammasToENDF6_MF12_13( MT, 12, endfMFList, flags, targetInfo, gammas )

    if len(targetInfo['MF6LCTs']) > 0:
        lcts = targetInfo['MF6LCTs']
        LCT = lcts[0]
        if LCT is None: LCT = 2
        for i in lcts:
            if i is None: continue
            if i != LCT: LCT = 3
        if 500 <= MT <= 572: LCT = 0
        if LR > 0:  # break-up reaction, check if we should use LCT=4
            if('implicitProduct' in targetInfo['ENDFconversionFlags'].get(self.outputChannel[1], '') and
                    self.outputChannel[1].distribution.evaluated.productFrame == xDataEnumsModule.Frame.centerOfMass):
                if lcts[0] == 2 and all( [lct==1 for lct in lcts[1:]] ):
                    LCT=4
            elif lcts[:2] == [2,2] and all( [lct==1 for lct in lcts[2:]] ):
                LCT=4

        MF6or26 = { 3 : 6, 23 : 26 }[targetInfo['crossSectionMF']]
        endfMFList[MF6or26][MT].insert( 0, endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], 0, LCT, len( targetInfo['MF6LCTs'] ), 0 ) )
        endfMFList[MF6or26][MT].append( endfFormatsModule.endfSENDLineNumber( ) )

    targetInfo['reaction'] = None

baseModule.Base_reaction.toENDF6 = toENDF6
