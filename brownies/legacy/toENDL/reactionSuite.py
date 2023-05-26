# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

from pqu import PQU as PQUModule

from PoPs import specialNuclearParticleID as specialNuclearParticleIDPoPsModule
from PoPs import IDs as IDsPoPsModule
from PoPs.chemicalElements import misc as miscGroupsPoPsModule
from PoPs.chemicalElements import chemicalElement as chemicalElementPoPsModule
from PoPs.decays import misc as miscDecaysPoPsModule
from PoPs.families import nuclide as nuclideModule

from fudge import enums as enumsModule
from fudge import styles as stylesModule
from fudge import outputChannel as outputChannelModule
from fudge import reactionSuite as reactionSuiteModule
from fudge.reactionData import crossSection as crossSectionModule
from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import CoulombPlusNuclearElastic as CoulombPlusNuclearElasticModule
from fudge.outputChannelData import Q as QModule
from fudge.productData import multiplicity as multiplicityModule
from fudge.productData.distributions import energy as energyModule
from fudge.productData.distributions import uncorrelated as uncorrelatedModule
from fudge.productData.distributions import reference as referenceModule

from brownies.legacy.endl import fudgeParameters
from brownies.legacy.endl import endlZA as endlZAClass
from brownies.legacy.endl import endl2 as endl2Module
from brownies.legacy.endl import endlmisc as endlmiscModule
from ..converting import endf_endl

fudgeParameters.VerboseMode = 0

energyUnit = 'MeV'

def getENDLFile( endlZA, yo, C, S, I ) :

    files = endlZA.findFiles( yo, C, I, S )
    if( len( files ) > 0 ) :
        if( len( files ) != 1 ) : raise Exception( 'Multiple yo, C, S, I = %d, %d, %d, %d files.' % ( yo, C, S, I ) )
        return( files[0] )
    return( endlZA.addFile( yo, C, I, S, halflife = None, printWarnings = False ) )

def processCoulombElastic(reactionSuite, reaction, endlZA, yi, temperature, muCutoff=0.94, accuracy=1e-3):

    def yoEqualYiOutput(yo, CValue, angualarData):
        """For internal use."""

        if yo > 1:
            I1File = getENDLFile(endlZA, yo+10, CValue, 0, 1)
            recoil = []
            for energyIn, POfMu in angualarData:
                POfMuRecoil = list(reversed([[-x if x != 0.0 else x, y] for x, y in POfMu]))
                recoil.append([energyIn, POfMuRecoil])
            I1File.addData(recoil, Q=0.0, X1=0.0, temperature=temperature, X4=0.0)

    particle = reactionSuite.PoPs[reaction.outputChannel.products[1].pid]
    if isinstance(particle, nuclideModule.Alias): particle = reactionSuite.PoPs[particle.pid]
    residualZA = miscGroupsPoPsModule.ZA(particle)
    try:
        yo = endl2Module.ZAToYo(residualZA)
    except:
        yo = 0

    muCutoffStyle = stylesModule.CoulombPlusNuclearElasticMuCutoff('C9', 'dummy', muCutoff)
    for doubleDifferentialCrossSection in reaction.doubleDifferentialCrossSection:
        if isinstance(doubleDifferentialCrossSection, CoulombPlusNuclearElasticModule.Form):
            crossSection, angular = doubleDifferentialCrossSection.processCoulombPlusNuclearMuCutoff(muCutoffStyle, excludeRutherfordScattering=True)
            if crossSection is not None:
                crossSectionData = [[x, y] for x, y in crossSection]
                I0File = getENDLFile(endlZA, 0, 9, 0, 0)
                I0File.addData(crossSectionData, Q=0.0, X1=0.0, temperature=temperature, X4=0.0)

                angualarData = []
                for xys1d in angular:
                    POfMu = [[x, y] for x, y in xys1d]
                    if POfMu[0][0] > -1.0: POfMu.insert(0, [-1.0, 0])
                    if POfMu[-1][0] < 1.0: POfMu.append([1.0, 0])
                    if POfMu[-3][1] == 0.0: POfMu.pop(-2)
                    angualarData.append([xys1d.outerDomainValue, POfMu])
                I1File = getENDLFile(endlZA, yi, 9, 0, 1)
                I1File.addData(angualarData, Q=0.0, X1=0.0, temperature=temperature, X4=0.0)

                yoEqualYiOutput(yo, 9, angualarData)

            crossSection = doubleDifferentialCrossSection.RutherfordScattering.crossSectionVersusEnergy(muCutoff, accuracy=accuracy)
            I0File = getENDLFile(endlZA, 0, 8, 0, 0)
            crossSectionData = [[x, y] for x, y in crossSection]
            I0File.addData(crossSectionData, Q=0.0, X1=0.0, temperature=temperature, X4=0.0)

            dSigma_dMu = doubleDifferentialCrossSection.RutherfordScattering.dSigma_dMuVersusEnergy(muCutoff, accuracy=accuracy)
            angualarData = []
            for xys1d in dSigma_dMu:
                xys1d = xys1d.normalize()
                POfMu = [[x, y] for x, y in xys1d]
                angualarData.append([xys1d.outerDomainValue, POfMu])
            I1File = getENDLFile(endlZA, yi, 8, 0, 1)
            I1File.addData(angualarData, Q=0.0, X1=0.0, temperature=temperature, X4=0.0)
            yoEqualYiOutput(yo, 8, angualarData)

def processC55S3Photons( reactionSuite, multiplicity, crossSection, C55S3Photons, photonBranchingDatas ) :

    initial_pid = multiplicity.product.parentProduct.pid
    final_pid = 'FIX ME'
    isotope = reactionSuite.PoPs[multiplicity.product.parentProduct.pid].isotope
    photonBranchingData = photonBranchingDatas[isotope.symbol]
    processC55S3Photons2( initial_pid, final_pid, 1.0, crossSection, C55S3Photons, photonBranchingData )

def processC55S3Photons2( initialLevel, finalLevel, probability, crossSection, C55S3Photons, photonBranchingData, level = 1 ) :

    if( initialLevel == finalLevel ) : return
    for branchingRatio, energy, endLevel, photonEmissionProbability in photonBranchingData[initialLevel]['photons'] :
        _probability = probability * branchingRatio
        energy = float( '%.8e' % energy.getValueAs( 'MeV' ) )            # Fix rounding issues.
        if( energy not in C55S3Photons ) :
            C55S3Photons[energy] = { 'crossSection' : crossSectionModule.XYs1d( data = [], axes = crossSectionModule.defaultAxes( energyUnit ) ), 'angular' : None }
        C55S3Photons[energy]['crossSection'] += _probability * crossSection
        processC55S3Photons2( endLevel, finalLevel, _probability, crossSection, C55S3Photons, photonBranchingData, level + 1 )

def processOutputChannel( reactionSuite, outputChannel, yos, endlZA, temperature, C, S, Q, X1, crossSection, C55S3Photons, photonBranchingDatas ) :

    superProducts = list( outputChannel )
    residual = None
    if( C != 15 ) :
        for superProduct in superProducts :
            if( superProduct.pid != IDsPoPsModule.photon ) : residual = superProduct 

    superProducts += list( outputChannel.fissionFragmentData.delayedNeutrons )

    for superProduct in superProducts :
        multiplicityCount = 1
        if( superProduct is residual ) :
            multiplicity = superProduct.multiplicity[0]
            multiplicityCount = min( 2, max( 1, int( multiplicity.evaluate( multiplicity.domainMin ) ) ) )
        for multiplicityIndex in range( 1, multiplicityCount + 1 ) :
            if superProduct.moniker == 'product' :
                product = superProduct
            elif superProduct.moniker == 'delayedNeutron' : 
                product = superProduct.product
                
            _S = S
            _X1 = X1
            if( ( C == 15 ) and ( product.pid == IDsPoPsModule.neutron ) ) :
                if( superProduct.moniker == 'delayedNeutron' ) :
                    _S = 7
                    _X1 = PQUModule.PQU( superProduct.rate[0].value, superProduct.rate[0].unit ).getValueAs( '1/s' )
            particle = reactionSuite.PoPs[product.pid]
            if( isinstance( particle, nuclideModule.Alias ) ) : particle = reactionSuite.PoPs[particle.pid]
            if( isinstance( particle, chemicalElementPoPsModule.ChemicalElement ) ) : continue
            nuclearLevelIndex = 0
            if isinstance( particle, nuclideModule.Particle): nuclearLevelIndex = particle.index
            productZA = miscGroupsPoPsModule.ZA( particle )
            if productZA < 2005 and nuclearLevelIndex == 0:
                productID = product.pid
                if( productID == IDsPoPsModule.photon ) : productID = 'gamma'
                yo = endl2Module.ZAToYo(productZA)
                if( _S != 7 ) :
                    if( ( C != 15 ) and ( superProduct is residual ) and ( multiplicityIndex == multiplicityCount ) ) : yo += 10
                    if( yo in yos ) :
                        print( 'WARNNIG: yo %s already exists: skipping.' % yo )
                        continue
                    yos.append( yo )
                if outputChannel.genre != enumsModule.Genre.twoBody:
                    multiplicity = product.multiplicity.toENDL( )
                    if( isinstance( multiplicity, multiplicityModule.Branching1d ) ) :
                        processC55S3Photons( reactionSuite, multiplicity, crossSection, C55S3Photons, photonBranchingDatas )
                    elif( multiplicity is not None ) :
                        I = 9
                        if( productID == IDsPoPsModule.neutron ) : I = 7
                        file = getENDLFile( endlZA, yo, C, _S, I )
                        data = file.addData( multiplicity, Q = Q, X1 = _X1, temperature = temperature )
                if( _S == 7 ) :
                    if( isinstance( product.distribution[0], referenceModule.Form ) ) :
                        distributionData = None
                    else :
                        distributionData = product.distribution.toENDL( )
                else :
                    distributionData = product.distribution.toENDL( )
                if( distributionData is not None ) :
                    for I in distributionData :
                        if( distributionData[I] is not None ) :
                            file = getENDLFile( endlZA, yo, C, _S, I )
                            data = file.addData( distributionData[I], Q = Q, X1 = _X1, temperature = temperature )
            else :
                if( product.outputChannel is not None ) :
                    processOutputChannel( reactionSuite, product.outputChannel, yos, endlZA, temperature, C, S, Q, X1, crossSection, C55S3Photons, photonBranchingDatas )

def toENDL( self, directory, verbose = 0 ) :

    self.convertUnits( { 'eV' : energyUnit } )

    temperature = self.styles[0].temperature.getValueAs( 'MeV/k' )

    projectileId = self.projectile
    if( projectileId == 'photon' ) : projectileId = 'gamma'
    yi = endlmiscModule.incidentParticleTags( projectileId )[0]

    suffix = ''
    if( 'FissionProductENDL9912' in self.target ) :
        ZA = int( self.target[-5:] )
    else :
        target = self.PoPs[self.target]
        if( isinstance( target, chemicalElementPoPsModule.ChemicalElement ) ) :
            ZA = 1000 * target.Z
        else :
            ZA = miscGroupsPoPsModule.ZA( target )
        if( isinstance( target, nuclideModule.Alias ) ) : suffix = 'm'

    endlZA = endlZAClass( ZA, yi, workDir = directory, suffix = suffix )

    halflife = endlZA.bdflsFile.halflife( endlZA.ZA )
    if( halflife is None ) : endlZA.bdflsFile.addZAHalflife( endlZA.ZA, 1e50, warning = False )

    os.system( 'rm -rf %s/*' % endlZA.workDir )

    if( verbose > 0 ) : print( '%5s    %-56s  C  S %-8s %10s %10s %4s' % ( 'MT', 'reaction label', 'genre', 'Q', 'X1', 'MT' ) )

    specialMTs = [ -21, -29, -1102, -1091 ]
    C55S3Photons = {}
    photonBranchingDatas = {}
    for chemicalElement in self.PoPs.chemicalElements :
        for isotope in chemicalElement.isotopes : photonBranchingDatas[isotope.symbol] = miscDecaysPoPsModule.photonBranchingData( self.PoPs, isotope.symbol )

    totalCrossSection = None
    if( yi == 1 ) : totalCrossSection = crossSectionModule.XYs1d( data = [], axes = crossSectionModule.defaultAxes( energyUnit ) )
    for reaction in self.reactions :
        ENDF_MT = reaction.ENDF_MT
        if ENDF_MT == 2 and yi in [2, 3, 4, 5, 6] and self.target != IDsPoPsModule.neutron:
            processCoulombElastic(self, reaction, endlZA, yi, temperature)
            continue
        if( verbose > 0 ) : print( '% 5s ' % ENDF_MT, end = '' )
        outputChannel = reaction.outputChannel
        if( ENDF_MT < 0 ) :
            if( ENDF_MT == -49 ) :
                CS_ENDF_MT = { 'MT' : ENDF_MT, 'C' : -ENDF_MT, 'S' : 0 }
            elif ENDF_MT == -20:
                CS_ENDF_MT = {'MT': ENDF_MT, 'C': 21, 'S': 0}
                if outputChannel.genre == enumsModule.Genre.twoBody:
                    CS_ENDF_MT['S'] = 1
            elif ENDF_MT == -30:
                CS_ENDF_MT = {'MT': ENDF_MT, 'C': 30, 'S': 0}
                if outputChannel.genre == enumsModule.Genre.twoBody:
                    CS_ENDF_MT['S'] = 1
            elif( ENDF_MT not in specialMTs ) :
                print( 'ERROR, conversion of MT = %s to C, S currently not supported -- skipping this reaction.' % ENDF_MT )
                continue
            elif( ENDF_MT < -999 ) :
                CS_ENDF_MT = { 'MT' : ENDF_MT, 'C' : 46, 'S' : 0 }
            else :
                CS_ENDF_MT = { 'MT' : ENDF_MT, 'C' : -ENDF_MT, 'S' : 0 }
        else:
            C, S = endf_endl.getCSFromMT( ENDF_MT )
            if outputChannel.genre == enumsModule.Genre.twoBody:
                residual = self.PoPs[outputChannel[1].pid]
                if isinstance(residual, nuclideModule.Particle) and ENDF_MT != 2: S = 1
            if ENDF_MT in (600, 601):
                if (yi == 3 and ZA == 4007) or (yi == 3 and ZA == 4007):
                    C = 48
            CS_ENDF_MT = { 'C' : C, 'S' : S, 'MT' : ENDF_MT }

        MT = CS_ENDF_MT['MT']
        C = CS_ENDF_MT['C']
        S = CS_ENDF_MT['S']
        if( MT == 2 ) : S = 0                       # Special case for meta-stables for which CS_ENDF_MT['S'] is 1.
        if( outputChannel.process is not None ) :
            if( outputChannel.process == 'ENDL:S2' ) : S = 2
        pid = specialNuclearParticleIDPoPsModule.specialNuclearParticleID(outputChannel.products[0].pid, specialNuclearParticleIDPoPsModule.Mode.nuclide)
        if( MT == 20 ) :
            if pid == "H1": C = 21
        elif( MT == 32 ) :
            if pid == "H2": C, S = 35, 1

        if( C == 15 ) :
            Q = 0.0
            if( len( outputChannel.Q ) > 0 ) :
                Q = outputChannel.Q.getConstant( )
                Q_form = outputChannel.Q[0]
                if( isinstance( Q_form, QModule.XYs1d ) ) :
                    file = getENDLFile( endlZA, 0, C, S, 12 )
                    data = file.addData( [ [ x, y ] for x, y in Q_form ], Q = Q, X1 = X1, temperature = temperature )
        else :
            Q = reaction.thresholdQAs( energyUnit )
        X1 = 0.0
        if outputChannel.genre == enumsModule.Genre.twoBody:
            residual = self.PoPs[outputChannel[1].pid]
            if( isinstance( residual, nuclideModule.Particle ) ) : X1 = residual.nucleus.energy[0].value
        if( verbose > 0 ) : print( '   %-56s % 2d %1d %-8s % 10g % 10g %4s' % ( reaction, C, S, outputChannel.genre, Q, X1, MT ) )

        if( C < 1 ) :
            print( 'WARNNIG: unsupported C = %d, S = %d, MT = %d' % ( C, S, MT ) )
            continue

        crossSectionData = None
        if( C in [ 71, 72, 73, 74 ] ) :
            crossSection = reaction.crossSection[0]
            if( isinstance( crossSection, crossSectionModule.Regions1d ) ) :
                crossSectionData = []
                for regionIndex, region in enumerate( crossSection ) :
                    if( regionIndex > 0 ) :
                        if( crossSectionData[-1] == region[0] ) : crossSectionData.pop( )
                    for x, y in region : crossSectionData.append( [ x, y ] )
        else :
            crossSection = reaction.crossSection.toPointwise_withLinearXYs( lowerEps = 1e-6, upperEps = 1e-6 )
        if( crossSectionData is None ) : crossSectionData = [ [ x, y ] for x, y in crossSection ]
        file = getENDLFile( endlZA, 0, C, S, 0 )
        X4 = 0.0
        if( C == 46 ) :
            residual = self.PoPs[outputChannel[1].pid]
            if( isinstance( residual, nuclideModule.Particle ) ) : X4 = residual.nucleus.energy[0].pqu( ).getValueAs( "MeV" )
        data = file.addData( crossSectionData, Q = Q, X1 = X1, temperature = temperature, X4 = X4 )

        if( totalCrossSection is not None ) :
            if( len( totalCrossSection ) > 0 ) :
                energyMin = crossSection.domainMin
                if( ( energyMin > totalCrossSection.domainMin ) and ( crossSection[0][1] != 0.0 ) ) :
                    crossSection.setValue( energyMin, 0.0 )
                elif( energyMin < totalCrossSection.domainMin ) :
                    crossSection = crossSection.trim( )
                    if( crossSection.domainMin < totalCrossSection.domainMin ) :
                        crossSection = crossSection.domainSlice( domainMin = totalCrossSection.domainMin )
                if( crossSection.domainMax > totalCrossSection.domainMax ) :
                    crossSection = crossSection.domainSlice( domainMax = totalCrossSection.domainMax )
            totalCrossSection += crossSection

        yos = []
        processOutputChannel( self, outputChannel, yos, endlZA, temperature, C, S, Q, X1, crossSection, C55S3Photons, photonBranchingDatas )

    if( totalCrossSection is not None ) :
        file = getENDLFile( endlZA, 0, 1, 0, 0 )
        totalCrossSection = [ [ x, y ] for x, y in totalCrossSection ]
        data = file.addData( totalCrossSection, Q = 0.0, X1 = 0.0, temperature = temperature )

    C55S0Photon = []
    if( len( self.orphanProducts ) > 0 ) :
        if( len( self.orphanProducts ) == 1 ) :
            orphanProduct = self.orphanProducts[0]
            crossSectionC55 = orphanProduct.crossSection[0].link[0].toPointwiseLinear( )
            for product in orphanProduct.outputChannel :
                multiplicity = product.multiplicity[0].toPointwiseLinear( ).trim( )
                crossSection = [ [ x, y * crossSectionC55.evaluate( x ) ] for x, y in multiplicity ]
                crossSection = crossSectionModule.XYs1d( data = crossSection, axes = crossSectionModule.defaultAxes( energyUnit ) )
                distribution = product.distribution[0]
                if( isinstance( distribution, uncorrelatedModule.Form ) ) :
                    energy = distribution.energySubform.data
                    angular = distribution.angularSubform.data
                    if( isinstance( energy, energyModule.DiscreteGamma ) ) :
                        C55S3Photons[energy.value] = { 'crossSection' : crossSection, 'angular' : angular, 'domainMin' : energy.domainMin, 'domainMax' : energy.domainMax }
                    elif( isinstance( energy, energyModule.XYs2d ) ) :
                        C55S0Photon.append( [ crossSection, distribution.toENDL( ) ] )
                    else :
                        print( '    WARNING: unsupported C55 photon energy distribution "%s".' % energy.moniker )
                else :
                    print( '    WARNING: unsupported C55 photon distribution "%s".' % distribution.moniker )
        else :
            print( '    WARNING: more than 1 orphanProduct is not supported: got %s' % len( self.orphanProducts ) )

    I1Isotropic = [ [ -1.0, 0.5 ], [ 1.0, 0.5 ] ]
    for energy in sorted( C55S3Photons ) :
        crossSection = C55S3Photons[energy]['crossSection']
        crossSectionData = [ [ x, y ] for x, y in crossSection ]
        file = getENDLFile( endlZA, 0, 55, 3, 0 )
        data = file.addData( crossSectionData, Q = 0.0, X1 = energy, temperature = temperature )

        file = getENDLFile( endlZA, 7, 55, 3, 1 )
        angular = C55S3Photons[energy]['angular']
        domainMin = C55S3Photons[energy]['domainMin']
        domainMax = C55S3Photons[energy]['domainMax']
        if( angular is None ) :
            I1Data = [ [ domainMin, I1Isotropic ], [ domainMax, I1Isotropic ] ]
        else :
            I1Data = angular.toENDL( )
        data = file.addData( I1Data, Q = 0.0, X1 = energy, temperature = temperature )

    if( len( C55S0Photon ) > 0 ) :
        if( len( C55S0Photon ) == 1 ) :
            crossSection, distributionData = C55S0Photon[0]
            file = getENDLFile( endlZA, 0, 55, 0, 0 )
            crossSectionData = [ [ x, y ] for x, y in crossSection ]
            data = file.addData( crossSectionData, Q = 0.0, X1 = 0.0, temperature = temperature )

            if( distributionData is not None ) :
                for I in distributionData :
                    if( distributionData[I] is not None ) :
                        file = getENDLFile( endlZA, 7, 55, 0, 4 )
                        data = file.addData( distributionData[I], Q = 0.0, X1 = 0.0, temperature = temperature )
        else :
            print( '    WARNING: more than 1 C=55 S=0 photon found: number = %s.' % len( C55S0Photon ) )

    for data in endlZA.findDatas( ) : data.setFormat( 15 )

    endlZA.source = endlZA.workDir

    sep = ''
    documentation = ''
    for style in self.styles:
        if isinstance(style, stylesModule.Evaluated):
            documentation += sep + style.documentation.body.body
            if len(documentation) > 0: sep = '\n'
    if len(documentation) == 0:
        sep = ''
        for style in self.styles:
            if isinstance(style, stylesModule.Evaluated):
                documentation += sep + style.documentation.endfCompatible.body
                sep = '\n'
    endlZA.setDocumentation(documentation)

    endlZA.save()

    return endlZA

reactionSuiteModule.ReactionSuite.toENDL = toENDL
