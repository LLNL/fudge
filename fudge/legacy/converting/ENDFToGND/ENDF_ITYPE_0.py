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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

import sys
import endfFileToGNDMisc as endfFileToGNDMiscModule
from ENDF_ITYPE_0_Misc import BadResonances, getTotalOrPromptFission, getDelayedFission, getFissionEnergies, readMF2, \
    readMF8, readMF8_454_459, parseReaction, parseCovariances
from fudge.particles import nuclear
from fudge.legacy.converting import endf_endl, toGNDMisc
from pqu import PQU as PQUModule
from fudge.gnd import tokens as tokensModule
from fudge.gnd import xParticle as xParticleModule
from fudge.gnd import sums as sumsModule
from fudge.gnd import channels as channelsModule
from fudge.gnd.reactions import reaction as reactionModule
from fudge.gnd.reactions import production as productionModule
from fudge.gnd.reactions import fissionComponent as fissionComponentModule
from fudge.gnd.reactionData import crossSection as crossSectionModule
from fudge.gnd.productData import multiplicity as multiplicityModule

import xData.XYs as XYsModule

def deriveMT3MF3FromMT1_2( info, reactionSuite ) :

    totalCrossSection, elasticCrossSection = None, None
    for reaction in reactionSuite.sums :
        if( reaction.ENDF_MT == 1 ) :
            totalCrossSection = reaction.crossSection.evaluated
            break

    for reaction in reactionSuite.reactions :
        if( reaction.ENDF_MT == 2 ) :
            elasticCrossSection = reaction.crossSection.evaluated
            break
    if( totalCrossSection is None ) : raise Exception( 'No total cross section for calculating non-elastic' )
    if( elasticCrossSection is None ) : raise Exception( 'No elastic cross section for calculating non-elastic' )

    form = totalCrossSection - elasticCrossSection
    form.label = info.style
    return( form )

def ITYPE_0( MTDatas, info, reactionSuite, singleMTOnly, MTs2Skip, parseCrossSectionOnly, doCovariances ) :

    warningList = []

    info.totalOrPromptFissionNeutrons = {}
    info.totalMF6_12_13Gammas = {}
    if( 452 in MTDatas ) :
        info.totalOrPromptFissionNeutrons['total'] = getTotalOrPromptFission( info, MTDatas[452][1], 'total', warningList )
        #MTDatas.pop( 452 ) # don't remove these yet, still need the covariance info
    if( 455 in MTDatas ) :
        info.indices.delayedFissionNeutron = True
        info.delayedFissionDecayChannel = getDelayedFission( info, MTDatas[455], warningList )
        info.indices.delayedFissionNeutron = False
        #MTDatas.pop( 455 )
    if( 456 in MTDatas ) :
        info.totalOrPromptFissionNeutrons[tokensModule.promptToken] = getTotalOrPromptFission( info, MTDatas[456][1], tokensModule.promptToken, warningList )
        #MTDatas.pop( 456 )
    if( 458 in MTDatas ) :
        info.fissionEnergies = getFissionEnergies( info, MTDatas[458], warningList )
        #MTDatas.pop( 458 )
    if ( 454 in MTDatas ) : 
        info.independentFissionYields = readMF8_454_459( info, 454, MTDatas[454], warningList )
    if ( 459 in MTDatas ) : 
        info.cumulativeFissionYields = readMF8_454_459( info, 459, MTDatas[459], warningList )
    sys.stdout.flush( )
    for warning in warningList : info.logs.write( "       WARNING: %s\n" % warning, stderrWriting = True )

    MTList = endfFileToGNDMiscModule.niceSortOfMTs( MTDatas.keys( ), verbose = False, logFile = info.logs )

    haveTotalFission = (18 in MTList)
    fissionMTs = [mt for mt in MTList if mt in (19,20,21,38)]

    summedReactions = {}
    summedReactionsInfo = { 4 : range( 50, 92 ), 103 : range( 600, 650 ), 104 : range( 650, 700 ), 105 : range( 700, 750 ), 106 : range( 750, 800 ), 107 : range( 800, 850 ) }
    for summedMT, partialReactions in summedReactionsInfo.items( ) :
        if( summedMT not in MTList ) : continue
        for MT in MTList :
            if( MT in partialReactions ) :
                summedReactions[summedMT] = None
                break

    for summedMT in ( 1, 3 ) :
        if( summedMT in MTList ) : summedReactions[summedMT] = None

    MT5Reaction = None
    reactions = []
    fissionComponents = []
    productions = []
    nonElastic = []
    delayInsertingSummedReaction = []

    for MT in MTList :
        if( MT in MTs2Skip ) : continue
        if( ( singleMTOnly is not None ) and ( MT != singleMTOnly ) ) : continue

        warningList = []
        MTData = MTDatas[MT]

        # Sometimes excited states are identified in MF8. Read this before reading distributions to make sure info is present.
        LMF, radioactiveDatas = readMF8( info, MT, MTData, warningList )

        doParseReaction = 3 in MTData
        if( not( doParseReaction ) ) :
            if( MT == 3 ) : doParseReaction = ( 12 in MTData ) or ( 13 in MTData )
        if( doParseReaction ) : # normal reaction, with cross section and distributions
            try :
                crossSection, outputChannel, MFKeys = parseReaction( info, info.target, info.projectileZA,
                        info.targetZA, MT, MTData, warningList, parseCrossSectionOnly = parseCrossSectionOnly )
            except KeyboardInterrupt:
                raise
            except:
                import traceback
                info.logs.write( traceback.format_exc( ), stderrWriting = True )
                info.doRaise.append( traceback.format_exc( ) )
                info.logs.write( '\n' )
                sys.stdout.flush( )
                continue

            info.logs.write( '\n' )
            sys.stdout.flush( )
            if( len( MFKeys ) ) :
                warningList.append( 'For reaction MT = %d, the following MFs were not converted: %s\n' % ( MT, MFKeys ) )
            if( outputChannel is None ) : break

            if( MT in summedReactions ) :
                summedReactions[MT] = [ crossSection, outputChannel ]
            else :
                if( MT != 2 ) : nonElastic.append( MT )
                reaction = reactionModule.reaction( outputChannel, "%s" % 0, ENDF_MT = MT, date = info.Date )
                reaction.crossSection.add( crossSection )
                if( MT == 5 ) :
                    MT5Reaction = reaction
                elif MT in fissionMTs and haveTotalFission: # this is 1st, 2nd, etc fission but total is also present
                    fissionComponents.append( reaction )
                else :
                    if( MT in summedReactionsInfo ) :
                        delayInsertingSummedReaction.append( reaction )
                    else :
                        reactions.append( [ MT, reaction ] )
        else :
            MFList = []
            for MF in [ 4, 5, 6, 12, 13, 14, 15 ] :
                if( MF in MTData ) : MFList.append( '%d' % MF )
            if( MFList != '' ) : warningList.append( 'MT = %d has MF = %s data and no MF 3 data' % ( MT, ', '.join( MFList ) ) )

        for radioactiveData in radioactiveDatas : # Get radioactive production data (if any) from MF 8-10. Cross section form depends on value of LMF.
            if( LMF in [ 3, 6, 9 ] ) :  # Cross section is reference to MF3.
                productionCrossSection = crossSectionModule.reference( link = reaction.crossSection, label = info.style )
            elif( LMF == 10 ) :         # MF10 data is cross section. Product's multipliticy is 1.
                productionCrossSection = radioactiveData[2]
            else :
                raise Exception( "Unknown LMF=%d encountered in MF=8 for MT=%d" % ( LMF, MT ) )

            Q = outputChannel.Q[info.style]
            if( LMF in [ 9, 10 ] ) :
                Q = toGNDMisc.returnConstantQ( info.style, radioactiveData[4] )

            if( LMF == 6 ) :      # Product multiplicity is in MF6, so production channel multiplicity needs to refer to it:
                MF6prod = outputChannel.getProductsWithName( radioactiveData[1].name )

                if( len( MF6prod ) != 1 ) : # CMM: I'm not sure if we need this test anymore
                    warningList.append( 'Unique MT%d radioactive product %s not found in product list!' %
                                ( MT, radioactiveData[1].name ) )
                    info.doRaise.append( warningList[-1] )
                    continue

                radioactiveData[1].multiplicity.add( multiplicityModule.reference( label = info.style, link = MF6prod[0].multiplicity ) )
                radioactiveData[1].multiplicity.remove( 'unknown' )   # 'unknown' was put in as a place-holder by readMF8

            productionOutputChannel = channelsModule.productionChannel( )
            productionOutputChannel.Q.add( Q )
            productionOutputChannel.products.add( radioactiveData[1] )

            production = productionModule.production( productionOutputChannel, label = -1, ENDF_MT = MT, date = info.Date )
            production.crossSection.add( productionCrossSection )
            productions.append( production )

        for warning in warningList : info.logs.write( "       WARNING: %s\n" % warning, stderrWriting = True )

    iChannel = 0
    for MT, reaction in reactions :
        reaction.setLabel( str( iChannel ) )
        reactionSuite.reactions.add( reaction )
        iChannel += 1

    for reaction in delayInsertingSummedReaction :
        reaction.setLabel( '%s' % iChannel )
        reactionSuite.reactions.add( reaction )
        iChannel += 1

    if( MT5Reaction is not None ) :
        MT5Reaction.setLabel( "%s" % iChannel )
        reactionSuite.reactions.add( MT5Reaction )
        iChannel += 1

# BRB, The channelIDs should be in a common area?????
    channelIDs = { 1 : 'total', 3 : 'nonelastic', 4 : '(z,n)', 103 : '(z,p)', 104 : '(z,d)', 105 : '(z,t)', 106 : '(z,He3)', 107 :'(z,alpha)' }
    if( 3 in summedReactions ) : summedReactionsInfo[3] = nonElastic
    if( ( 1 in summedReactions ) and ( 2 in MTList ) ) : summedReactionsInfo[1] = [ 2 ] + nonElastic
    summedReactionMTs = endfFileToGNDMiscModule.niceSortOfMTs( summedReactions.keys( ), verbose = False, logFile = info.logs )
    for MT in ( 4, 3, 1 ) :
        if( MT in summedReactionMTs ) :
            summedReactionMTs.remove( MT )
            summedReactionMTs.insert( 0, MT )
    for i1, MT in enumerate( summedReactionMTs ) :
        if( summedReactions[MT] is None ) : continue
        crossSection, outputChannel = summedReactions[MT]
        if( ( MT == 3 ) and ( crossSection is None ) ) : 
            crossSection = deriveMT3MF3FromMT1_2( info, reactionSuite )
        summands = [ sumsModule.add( link = r.crossSection ) for r in reactionSuite.reactions if r.ENDF_MT in summedReactionsInfo[MT] ]
        summedCrossSection = sumsModule.crossSectionSum( name=channelIDs[MT], label = str( i1 ), ENDF_MT = MT,
                summands = sumsModule.listOfSummands( summandList = summands ), date = info.Date )
        summedCrossSection.Q.add( outputChannel.Q[info.style] )
        summedCrossSection.crossSection.add( crossSection )
        reactionSuite.sums.add( summedCrossSection )

        gammas = []
        for product in outputChannel :
            if( isinstance( product.particle, xParticleModule.photon ) ) :
                gammas.append( product )
            else :
                if( product.decayChannel is not None ) :
                    for product2 in product.decayChannel :
                        if( isinstance( product2.particle, xParticleModule.photon ) ) : gammas.append( product2 )
        if( len( gammas ) > 0 ) :
            productChannel = channelsModule.NBodyOutputChannel( )
            for QForm in outputChannel.Q : productChannel.Q.add( QForm )
            for gamma in gammas : productChannel.products.add( gamma )
            productionReaction = productionModule.production( productChannel, str( iChannel ), MT, date = info.Date )
            crossSectionLink = crossSectionModule.reference( link = summedCrossSection.crossSection, label = info.style )
            productionReaction.crossSection.add( crossSectionLink )
            reactionSuite.reactions.add( productionReaction )
            iChannel += 1

    for i1, reaction in enumerate( fissionComponents ) :  # 1st-chance, 2nd-chance, etc. Convert them to fissionComponent instances:
        fissionComponent = fissionComponentModule.fissionComponent( reaction.outputChannel, str( i1 ), 
                reaction.ENDF_MT, date = reaction.date )
        for crossSection in reaction.crossSection : fissionComponent.crossSection.add( crossSection )
        reactionSuite.fissionComponents.add( fissionComponent )

    for i1, production in enumerate( productions ) :
        production.setLabel( "%s" % i1 )
        reactionSuite.productions.add( production )

    warningList = []
    try :               # Parse resonance section.
        mf2 = None
        if( 151 in MTDatas and not parseCrossSectionOnly ) :
            mf2 = MTDatas.get( 151 ).get( 2 )    # Resonance data available.
        if( mf2 ) :
            info.logs.write( '    Reading resonances (MF=2 MT=151)\n' )
            resonances, resonanceMTs = readMF2( info, mf2, warningList )
            resonances.reconstructCrossSection = ( info.LRP == 1 )   # LRP was read in from first line of ENDF file
            if resonances.unresolved and not resonances.resolved:
                if resonances.unresolved.tabulatedWidths.forSelfShieldingOnly:
                    resonances.reconstructCrossSection = False
            reactionSuite.addResonances( resonances )

            # add spins appearing in resonance region to the particle list
            for particle, spinParity in info.particleSpins.items():
                if particle=='target': particle = reactionSuite.target
                else: particle = reactionSuite.getParticle( particle )
                if isinstance(particle, xParticleModule.isotope) and particle.name not in ('gamma','n','H1'):
                    # spin should be associated with ground state level:
                    particle.addLevel( xParticleModule.nuclearLevel( name = particle.name + '_e0',
                        energy = PQUModule.PQU( '0 eV' ), label = 0 ) )
                    particle = particle.levels[0]
                particle.attributes['spin'] = spinParity[0]
                if spinParity[1]: particle.attributes['parity'] = spinParity[1]

            if resonances.reconstructCrossSection:
                # modify cross sections for relevant channels to indicate resonance contribution is needed:
                resonanceLink = crossSectionModule.resonanceLink( link = resonances )

                for MT in resonanceMTs :
                    MTChannels  = [ r1 for r1 in reactionSuite.reactions         if( r1.getENDL_CS_ENDF_MT()['MT'] == MT )
                                    and isinstance(r1, reactionModule.reaction) ]
                    MTChannels += [ r1 for r1 in reactionSuite.sums   if( r1.ENDF_MT == MT ) ]
                    MTChannels += [ r1 for r1 in reactionSuite.fissionComponents if( r1.getENDL_CS_ENDF_MT()['MT'] == MT ) ]
                    if( len( MTChannels ) == 0 ) :
                        if( MT in ( 3, 18, 19 ) ) :
                            continue
                        else :
                            warningList.append( 'Unable to find channel corresponding to resonance data for MT%d' % MT )
                    elif( len( MTChannels ) == 1 ) :
                        crossSectionComponent = MTChannels[0].crossSection
                        backgroundForm = crossSectionComponent[info.style]
                        crossSectionComponent.remove( backgroundForm.label )
                        crossSectionComponent.add( crossSectionModule.resonancesWithBackground( info.style, backgroundForm, resonanceLink ) )
                    else :
                        raise 'hell - FIXME'                # This may not be correct.
                        crossSectionComponent = MTChannels[0].crossSection
                        backgroundComponent = crossSectionComponent[info.style].crossSection
                        backgroundForm = backgroundComponent[info.style]
                        backgroundComponent.remove( backgroundForm.label )
                        referredXSecForm = crossSectionModule.resonancesWithBackground( info.style, backgroundForm, resonanceLink )
                        backgroundComponent.add( referredXSecForm )
    except BadResonances as e:
        warningList.append( '       ERROR: unable to parse resonances! Error message: %s' % e )
        info.doRaise.append( warningList[-1] )

    if( doCovariances ) :
        covarianceMFs = sorted( set( [mf for mt in MTDatas.values() for mf in mt.keys() if mf>30] ) )
        if covarianceMFs:
            info.logs.write( '    Reading covariances (MFs %s)\n' % ','.join(map(str,covarianceMFs) ) )
        try:
            """ parse covariances. This also requires setting up links from data to covariances, so we
            must ensure the links are synchronized """

            MTdict = {}
            for reaction in reactionSuite :
                MT = reaction.ENDF_MT
                if MT in MTdict:
                    MTdict[MT].append( reaction )
                else:
                    MTdict[MT] = [reaction]
            covarianceSuite, linkData = parseCovariances( info, MTDatas, MTdict, singleMTOnly = singleMTOnly,
                    resonances = getattr( reactionSuite, 'resonances', None ) )
            if( len( covarianceSuite ) > 0 ) :
                covarianceSuite.target = str(info.target)
                covarianceSuite.projectile = str(info.projectile)
                covarianceSuite.styles.add( info.evaluatedStyle )
                #covarianceSuite.removeExtraZeros() # disable for easier comparison to ENDF
            else :
                covarianceSuite = None
        except Exception as e:
            warningList.append( "Couldn't parse covariances! Error message: %s" % e )
            info.doRaise.append( warningList[-1] )
            covarianceSuite = None
            raise
    else :
        covarianceSuite = None
    
    for ZA in info.ZA_AWRMasses :
        mostCount = 0
        for AWR in info.ZA_AWRMasses[ZA] :
            if( info.ZA_AWRMasses[ZA][AWR] > mostCount ) :
                mostCount = info.ZA_AWRMasses[ZA][AWR]
                mostAWR = AWR
        info.ZA_AWRMasses[ZA] = mostAWR * info.masses.getMassFromZA( 1 )

    if( info.level > 0 ) : # AWR is for isomer mass. Adjust info.ZAMasses to GS mass:
        info.ZA_AWRMasses[info.targetZA] -= PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( info.level ),'eV/c**2').getValueAs('amu')

    for particle in reactionSuite.particles :                        # Fix up any particle whose mass is not defined (i.e., is None).
        if( particle.name == 'gamma' ) : continue
        ZA = nuclear.getZ_A_suffix_andZAFromName( particle.name )[3]
        mass = particle.getMass( 'amu' )
        if( ZA in info.ZA_AWRMasses ) :
            mass = info.ZA_AWRMasses[ZA]
        elif( ( mass is None ) or ( ZA in info.ZAMasses ) ) :
            if( info.ZAMasses[ZA] is None ) :
                mass = info.masses.getMassFromZA( ZA )
            else :
                mass = abs( info.ZAMasses[ZA] )
        particle.setMass( PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( mass ), 'amu' ) )
            # also add ground state level if there are discrete excited levels.
        levels = particle.levels.keys( )
        if( 's' in levels ) : levels.remove( 's' )
        if( 'c' in levels ) : levels.remove( 'c' )
        if( ( len( levels ) > 0 ) and ( 0 not in levels ) ) :
# BRB, Caleb, we should not hardwire the '_e0' here, need to implements naming in the xParticle.py module.
            particle.addLevel( xParticleModule.nuclearLevel( name = particle.name + '_e0',
                    energy = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( 0 ), 'eV' ), label = 0 ) )

    MF12BaseMTsAndRange = [ [ 50, 92 ], [ 600, 650 ], [ 650, 700 ], [ 700, 750 ], [ 750, 800 ], [ 800, 850 ] ]

    if( singleMTOnly is None ) :
        branchingData = None
        #if( len( info.MF12_LO2 ) > 0 ) : reactionSuite.gammaBranching = {}
        for MTLO2, MF12_LO2 in sorted(info.MF12_LO2.items()) :  # The logic below assumes MTs are in ascending order per base MT.
            branchingBaseMT = None
            for MTBase, MTEnd in MF12BaseMTsAndRange :             # Determine base MT for this MTLO2
                if( MTBase < MTLO2 < MTEnd ) :
                    branchingBaseMT = MTBase
                    break
            if( branchingBaseMT is not None ) :
                residualZA = endf_endl.ENDF_MTZAEquation( info.projectileZA, info.targetZA, branchingBaseMT )[0][-1]
                residual = toGNDMisc.getTypeNameENDF( info, residualZA, None )
                residualName = residual.name
                level = MTLO2 - branchingBaseMT
                levelName, levelEnergy = '_e%d' % level, MF12_LO2[0]['ES']
                fullName = residualName + levelName
                levelEnergy_eV = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( levelEnergy ), 'eV' )
                    # compare this value to level energy from the particle list (from MF3 Q-value).
                particleLevelEnergy_eV = reactionSuite.getParticle(fullName).energy
                if( levelEnergy_eV.value != particleLevelEnergy_eV.value ) :
                    if( particleLevelEnergy_eV.value < 1e-12 ) :
                        warningList.append( "MF12 parent level energy (%s) set to zero?" % particleLevelEnergy_eV )
                        info.doRaise.append( warningList[-1] )
                    elif( abs( levelEnergy_eV - particleLevelEnergy_eV ) / particleLevelEnergy_eV < 1e-4 ) :
                        MFLabel = '3'
                                                                                            # Value with most precision wins.
                        str1 = PQUModule.floatToShortestString( levelEnergy_eV * 1e-20 )          # 1e-20 to insure e-form is returned.
                        str2 = PQUModule.floatToShortestString( particleLevelEnergy_eV * 1e-20 )  # Want 1.23e-16 and not 12300 to differ
                        if( len( str1 ) > len( str2 ) ) :                                   # length from 1.2345e-16 and not 12345.
                            reactionSuite.getParticle( fullName ).energy = levelEnergy_eV
                            MFLabel = '12'
                        warningList.append( "MF12 level energy %s differs from MF3 value %s. Setting to MF%s value." % \
                                ( levelEnergy_eV, particleLevelEnergy_eV, MFLabel ) )
                    else :
                        warningList.append( "MF12 parent level energy (%s) doesn't match known level" % particleLevelEnergy_eV )
                        info.doRaise.append( warningList[-1] )
                for MF12 in MF12_LO2 :
                    try :
                        finalLevelEnergy = MF12['ESk']
                        if( finalLevelEnergy > 0. ) :
                                # find particle in the particleList with energy == finalLevelEnergy
                            finalParticles = [ lev for lev in reactionSuite.getParticle( residualName ).levels.values()
                                    if lev.energy.getValueAs('eV') == finalLevelEnergy ]
                            if( len( finalParticles ) == 1 ) : finalParticle = finalParticles[0]
                            else :   # no exact match, look for levels within .01% of the exact value.
                                idx = 0
                                while( True ) :
                                    idx += 1
                                    finalParticleName = residualName+'_e%i'%idx
                                    if( not reactionSuite.hasParticle( finalParticleName ) ) :
                                        warningList.append( "MF12 final level energy (%s eV) doesn't match known level when decaying out of level %s " % \
                                                ( finalLevelEnergy, MTLO2 ) )
                                        info.doRaise.append( warningList[-1] )
                                    thisLevelEnergy = reactionSuite.getParticle(finalParticleName).energy.getValueAs('eV')
                                    if( abs( thisLevelEnergy - finalLevelEnergy ) < 1e-4 * finalLevelEnergy ) :
                                        finalParticle = reactionSuite.getParticle(finalParticleName)
                                        break   # found it
                        else : finalParticle = reactionSuite.getParticle(residualName+'_e0')
                        gammaTransition = 1.
                        if( len( MF12['branching'] ) > 2 ) : gammaTransition = MF12['branching'][1]
                        gamma = xParticleModule.nuclearLevelGamma( finalParticle, MF12['angularSubform'], MF12['branching'][0], 1-gammaTransition )
                        reactionSuite.getParticle( fullName ).addGamma( gamma )
                    except Exception, err :
                        warningList.append( 'raise somewhere in "for MF12 in MF12_LO2" loop: MT%d, %s' % ( MT, str( err ) ) )
                        info.doRaise.append( warningList[-1] )
            else :
                raise Exception( "Could not determine base MT for MF=12's MT=%s" % MTLO2 )

    sys.stdout.flush( )
    for warning in warningList : info.logs.write( "       WARNING: %s\n" % warning, stderrWriting = True )

    if( reactionSuite.resonances is not None and reactionSuite.resonances.reconstructCrossSection ):
        info.logs.write( '    Reconstructing resonances\n' )
        reactionSuite.reconstructResonances( styleName=info.reconstructedStyle,
                accuracy=info.reconstructedAccuracy, verbose=True, thin=True )

    def adjustMF13Multiplicity2( multiplicity, crossSection ) :

        energyMultiplicity = []
        if( multiplicity.domainMax( ) > crossSection.domainMax( ) ) :
                multiplicity = multiplicity.domainSlice( domainMax = crossSection.domainMax( ) )
        for energyIn, multiplicityValue in multiplicity :
            crossSectionAtEnergy = crossSection.evaluate( energyIn )
            if( crossSectionAtEnergy != 0 ) : multiplicityValue /= crossSectionAtEnergy
            energyMultiplicity.append( [ energyIn, multiplicityValue ] )
        multiplicity.setData( energyMultiplicity )

    def adjustMF13Multiplicity( multiplicity, crossSection ) :

        if( isinstance( multiplicity, multiplicityModule.pointwise ) ) :
            adjustMF13Multiplicity2( multiplicity, crossSection )
        elif( isinstance( multiplicity, multiplicityModule.piecewise ) ) :
            for region in multiplicity : adjustMF13Multiplicity2( region, crossSection )
        else :
            raise Exception( 'Unsupported multiplicity type "%s"' % multiplicity.moniker )

    for reaction in reactionSuite.reactions :
        MT = reaction.ENDF_MT
        crossSection = None
        allproducts = list(reaction.outputChannel)
        for prod in reaction.outputChannel:
            if prod.decayChannel is not None:
                allproducts.extend( list(prod.decayChannel) )
        for product in allproducts :
            multiplicity = product.multiplicity[info.style]
            if( hasattr( multiplicity, '_temp_divideByCrossSection' ) ) :
                if( crossSection is None ) : crossSection = reaction.crossSection.toPointwise_withLinearXYs( )
                adjustMF13Multiplicity( multiplicity, crossSection )
                del multiplicity._temp_divideByCrossSection
        if( MT in info.totalMF6_12_13Gammas ) :
            MF, multiplicity = info.totalMF6_12_13Gammas[MT]
            if( MF == 13 ) : adjustMF13Multiplicity( multiplicity, crossSection )
# BRB, Caleb, why the ',' after gammaProduction.
            gammaProduction, = [ tmp for tmp in reactionSuite.reactions if tmp.ENDF_MT == MT ]
            summands = [ sumsModule.add( link = r.multiplicity ) for r in gammaProduction.outputChannel.getProductsWithName('gamma') ]
            if len(summands)==0:
                for _product in gammaProduction.outputChannel:
                    if _product.decayChannel is not None:
                        summands += [ sumsModule.add( link = r.multiplicity ) for r in _product.decayChannel.getProductsWithName('gamma') ]
            if MT in channelIDs: name = channelIDs[MT]
            else: name = str(gammaProduction.outputChannel)

            multiplicitySum = sumsModule.multiplicitySum( name = name + " total gamma multiplicity",
                    label = str( len( reactionSuite.sums ) + 1 ), ENDF_MT = MT, summands = sumsModule.listOfSummands( summands ) )
            multiplicitySum.multiplicity.add( multiplicity )
            reactionSuite.sums.add( multiplicitySum )

    return( covarianceSuite )
