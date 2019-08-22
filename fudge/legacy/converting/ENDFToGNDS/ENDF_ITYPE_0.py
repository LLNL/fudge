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

import sys
import fractions

from pqu import PQU as PQUModule

from PoPs import IDs as IDsPoPsModule
from PoPs.quantities import quantity as quantityModule
from PoPs.quantities import mass as massModule
from PoPs.quantities import spin as spinModule
from PoPs.quantities import parity as parityModule
from PoPs.decays import decayData as decayDataModule
from PoPs.decays import probability as probabilityModule
from PoPs.decays import product as productModule
from PoPs.groups import misc as chemicalElementMiscModule
from PoPs.families import nuclide as nuclideModule
from PoPs.families import gaugeBoson as gaugeBosonModule

from xData import standards as standardsModule

import endfFileToGNDSMisc as endfFileToGNDSMiscModule
from ENDF_ITYPE_0_Misc import BadResonances, getTotalOrPromptFission, getDelayedFission, readMF2, \
    readMF8, parseReaction, parseCovariances
from fudge.gnds import tokens as tokensModule
from fudge.gnds import externalFile as externalFileModule
from fudge.gnds import sums as sumsModule
from fudge.gnds import channels as channelsModule
from fudge.gnds.reactions import reaction as reactionModule
from fudge.gnds.reactions import production as productionModule
from fudge.gnds.reactions import incompleteReaction as incompleteReactionModule
from fudge.gnds.reactions import fissionComponent as fissionComponentModule
from fudge.gnds.reactionData import crossSection as crossSectionModule
from fudge.gnds.productData import multiplicity as multiplicityModule
from fudge.gnds.productData.distributions import unspecified as unspecifiedModule
from fudge.gnds.productData.distributions import branching3d as branching3dModule

from .. import toGNDSMisc as toGNDSMiscModule
from .. import endf_endl as endf_endlModule

def deriveMT3MF3FromMT1_2( info, reactionSuite ) :

    totalCrossSection, elasticCrossSection = None, None
    for reaction in reactionSuite.sums.crossSections :
        if( reaction.ENDF_MT == 1 ) :
            totalCrossSection = reaction.crossSection.evaluated
            break

    for reaction in reactionSuite.reactions :
        if( reaction.ENDF_MT == 2 ) :
            elasticCrossSection = reaction.crossSection.evaluated
            break
    if( totalCrossSection is None ) : raise Exception( 'No total cross section for calculating non-elastic' )
    if( elasticCrossSection is None ) : raise Exception( 'No elastic cross section for calculating non-elastic' )

    try :
        form = totalCrossSection - elasticCrossSection
    except :
        totalCrossSection = totalCrossSection.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 1e-6 )
        elasticCrossSection = elasticCrossSection.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 1e-6 )
        form = totalCrossSection - elasticCrossSection

    form.label = info.style
    return( form )

def ITYPE_0( MTDatas, info, reactionSuite, singleMTOnly, MTs2Skip, parseCrossSectionOnly, doCovariances, verbose, reconstructResonances=True ) :

    warningList = []

    info.totalOrPromptFissionNeutrons = {}
    info.totalMF6_12_13Gammas = {}
    if( 452 in MTDatas ) :
        info.totalOrPromptFissionNeutrons['total'] = getTotalOrPromptFission( info, MTDatas[452][1], 'total', warningList )
        #MTDatas.pop( 452 ) # don't remove these yet, still need the covariance info
    if( 455 in MTDatas ) :
        info.delayedFissionDecayChannel = getDelayedFission( info, MTDatas[455], warningList )
        #MTDatas.pop( 455 )
    if( 456 in MTDatas ) :
        info.totalOrPromptFissionNeutrons[tokensModule.promptToken] = getTotalOrPromptFission( info, MTDatas[456][1], tokensModule.promptToken, warningList )
        #MTDatas.pop( 456 )
    if( 458 in MTDatas ) :
        info.fissionEnergyReleaseData = MTDatas[458]
        #MTDatas.pop( 458 )
    sys.stdout.flush( )
    for warning in warningList : info.logs.write( "       WARNING: %s\n" % warning, stderrWriting = True )

    MTList = endfFileToGNDSMiscModule.niceSortOfMTs( MTDatas.keys( ), verbose = False, logFile = info.logs )

    haveTotalFission = (18 in MTList)
    fissionMTs = [mt for mt in MTList if mt in (19,20,21,38)]

    continuumMTs = { 91 : range( 50, 91 ), 649 : range( 600, 649 ), 699 : range( 650, 699 ), 749: range( 700, 749 ), 799 : range( 750, 799 ), 849 : range( 800, 849 ) }
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

    try :
        inelasticMT = { 1 : 91, 1001 : 649, 1002 : 699, 1003 : 749, 2003 : 799, 2004 : 849 }[info.projectileZA]
    except :
        inelasticMT = -1

    MT5Reaction = None
    reactions = []
    fissionComponents = []
    productions = []
    incompleteReactions = []
    nonElastic = []
    delayInsertingSummedReaction = []
    linksToCheck = []   # links that may need to be updated after reading resonances

    for MT in MTList :
        if( MT in MTs2Skip ) : continue
        if( ( singleMTOnly is not None ) and ( MT != singleMTOnly ) ) : continue

        channelProcess = None
        if( MT == inelasticMT ) :
            channelProcess = channelsModule.processes.continuum
        else :
            if( MT in continuumMTs ) :
                for discreteMT in continuumMTs[MT] :
                    if( discreteMT in MTList ) : channelProcess = channelsModule.processes.continuum

        warningList = []
        MTData = MTDatas[MT]

        # Sometimes excited states are identified in MF8. Read this before reading distributions to make sure info is present.
        LMF, radioactiveDatas = readMF8( info, MT, MTData, warningList )

        doParseReaction = 3 in MTData
        if( not( doParseReaction ) ) :
            if( MT == 3 ) : doParseReaction = ( 12 in MTData ) or ( 13 in MTData )
        if( doParseReaction ) : # normal reaction, with cross section and distributions
            try :
                crossSection, outputChannel, MFKeys = parseReaction( info, info.target, info.projectileZA, info.targetZA, 
                        MT, MTData, warningList, parseCrossSectionOnly = parseCrossSectionOnly, channelProcess = channelProcess )
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
                reaction = reactionModule.reaction( outputChannel, ENDF_MT = MT )
                if( hasattr( info, 'dSigma_form' ) ) :
                    reaction.doubleDifferentialCrossSection.add( info.dSigma_form )
                    del info.dSigma_form
                    crossSection = crossSectionModule.CoulombPlusNuclearElastic(
                            link = reaction.doubleDifferentialCrossSection.evaluated,
                            label = info.style, relative = True )
                reaction.crossSection.add( crossSection )
                if( MT == 5 ) :
                    MT5Reaction = reaction
                elif MT in fissionMTs and haveTotalFission: # this is 1st, 2nd, etc fission but total is also present

                    from fudge.gnds.channelData.fissionEnergyReleased import fissionEnergyReleased
                    if isinstance( reaction.outputChannel.Q.evaluated, fissionEnergyReleased ):
                        Qcomponent = reaction.outputChannel.Q
                        qval = toGNDSMiscModule.returnConstantQ( info.style,
                            Qcomponent.evaluated.nonNeutrinoEnergy.data.coefficients[0],
                            crossSection )
                        Qcomponent.remove( info.style )
                        Qcomponent.add( qval ) # just put the approximate constant Q-value on 1st-chance, 2nd-chance etc.

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
            if( MFList != [] ) : warningList.append( 'MT = %d has MF = %s data and no MF 3 data' % ( MT, ', '.join( MFList ) ) )

        for radioactiveData in radioactiveDatas : # Get radioactive production data (if any) from MF 8-10. Cross section form depends on value of LMF.
            if( LMF in [ 3, 6, 9 ] ) :  # Cross section is reference to MF3.
                productionCrossSection = crossSectionModule.reference( link = reaction.crossSection.evaluated, label = info.style )
                linksToCheck.append( productionCrossSection )
            elif( LMF == 10 ) :         # MF10 data is cross section. Product's multipliticy is 1.
                productionCrossSection = radioactiveData[4]
            else :
                raise Exception( "Unknown LMF=%d encountered in MF=8 for MT=%d" % ( LMF, MT ) )

            ZAP = radioactiveData[0]
            ELFS = radioactiveData[1]
            LFS = radioactiveData[2]

            Q = outputChannel.Q[info.style]
            if( LMF in [ 9, 10 ] ) :
                Q = toGNDSMiscModule.returnConstantQ( info.style, radioactiveData[6], productionCrossSection )

            if MT==18:
                outputChannel = channelsModule.fissionChannel( channelsModule.fissionGenreTotal )
                outputChannel.Q.add(Q)
                icr = incompleteReactionModule.incompleteReaction(outputChannel,MT,label="fission")
                icr.crossSection.add( productionCrossSection )
                incompleteReactions.append( icr )
                continue

            if( LMF == 6 ) :      # Product multiplicity is in MF6, so production channel multiplicity needs to refer to it:
                residual = toGNDSMiscModule.getTypeNameGamma( info, ZAP, level = ELFS, levelIndex = LFS )
                MF6prod = outputChannel.getProductsWithName( residual.id )

                if( len( MF6prod ) != 1 ) : # problem appears in JEFF-3.2 Y90 and Y91
                    warningList.append( 'Unique MT%d radioactive product %s not found in product list!' %
                                ( MT, residual.id ) )
                    info.doRaise.append( warningList[-1] )
                    continue

                multiplicity = multiplicityModule.reference( label = info.style, link = MF6prod[0].multiplicity )
            else :
                multiplicity = radioactiveData[3]

            try :
                residual = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getTypeNameGamma( info, ZAP, level = ELFS, levelIndex = LFS ),
                        crossSection, multiplicity = multiplicity )
            except :
                info.logs.write( '\nMT = %s\n' % MT )
                raise

            productionOutputChannel = channelsModule.productionChannel( )
            productionOutputChannel.Q.add( Q )
            productionOutputChannel.products.add( productionOutputChannel.products.uniqueLabel( residual ) )
            productionOutputChannel.process = "%s%s" % (reactionSuite.target,
                        endf_endlModule.endfMTtoC_ProductLists[MT].reactionLabel.replace('z,', reactionSuite.projectile+',') )

            production = productionModule.production( productionOutputChannel, ENDF_MT = MT )
            production.crossSection.add( productionCrossSection )
            productions.append( production )

        for warning in warningList : info.logs.write( "       WARNING: %s\n" % warning, stderrWriting = True )

    for MT, reaction in reactions :
        reactionSuite.reactions.add( reaction )

    for reaction in delayInsertingSummedReaction :
        reactionSuite.reactions.add( reaction )

    if( MT5Reaction is not None ) :
        reactionSuite.reactions.add( MT5Reaction )

# BRB, The channelIDs should be in a common area?????
    channelIDs = { 1 : 'total', 3 : 'nonelastic', 4 : '(z,n)', 103 : '(z,p)', 104 : '(z,d)', 105 : '(z,t)', 106 : '(z,He3)', 107 :'(z,alpha)' }
    if( 3 in summedReactions ) : summedReactionsInfo[3] = nonElastic
    if( ( 1 in summedReactions ) and ( 2 in MTList ) ) : summedReactionsInfo[1] = [ 2 ] + nonElastic
    summedReactionMTs = endfFileToGNDSMiscModule.niceSortOfMTs( summedReactions.keys( ), verbose = False, logFile = info.logs )
    for MT in ( 4, 3, 1 ) :
        if( MT in summedReactionMTs ) :
            summedReactionMTs.remove( MT )
            summedReactionMTs.insert( 0, MT )
    for i1, MT in enumerate( summedReactionMTs ) :
        if( summedReactions[MT] is None ) : continue
        crossSection, outputChannel = summedReactions[MT]
        omitWhenWritingENDF = False
        if( ( MT == 3 ) and ( crossSection is None ) ) : 
            crossSection = deriveMT3MF3FromMT1_2( info, reactionSuite )
            omitWhenWritingENDF = True
        summands = [ sumsModule.add( link = r.crossSection ) for r in reactionSuite.reactions if r.ENDF_MT in summedReactionsInfo[MT] ]
        summedCrossSection = sumsModule.crossSectionSum( label = channelIDs[MT], ENDF_MT = MT,
                summands = sumsModule.listOfSummands( summandList = summands ) )
        summedCrossSection.Q.add( outputChannel.Q[info.style] )
        summedCrossSection.crossSection.add( crossSection )
        if( omitWhenWritingENDF ):
            info.ENDFconversionFlags.add( summedCrossSection, "omit" )
        reactionSuite.sums.crossSections.add( summedCrossSection )

        gammas = []
        for product in outputChannel :
            particle = reactionSuite.PoPs[product.id]
            if( isinstance( particle, gaugeBosonModule.particle ) ) :
                gammas.append( product )
            else :
                if( product.outputChannel is not None ) :
                    for product2 in product.outputChannel :
                        particle = reactionSuite.PoPs[product2.id]
                        if( isinstance( particle, gaugeBosonModule.particle ) ) : gammas.append( product2 )
        if( len( gammas ) > 0 ) :
            productChannel = channelsModule.NBodyOutputChannel( )
            for QForm in outputChannel.Q : productChannel.Q.add( QForm )
            for gamma in gammas : productChannel.products.add( productChannel.products.uniqueLabel( gamma ) )
            productionReaction = reactionModule.reaction( productChannel, ENDF_MT = MT, label = str( i1 ) )
            crossSectionLink = crossSectionModule.reference( link = summedCrossSection.crossSection.evaluated, label = info.style )
            linksToCheck.append( crossSectionLink )
            productionReaction.crossSection.add( crossSectionLink )
            reactionSuite.orphanProducts.add( productionReaction )

    for i1, reaction in enumerate( fissionComponents ) :  # 1st-chance, 2nd-chance, etc. Convert them to fissionComponent instances:
        fissionComponent = fissionComponentModule.fissionComponent( reaction.outputChannel, reaction.ENDF_MT )
        for crossSection in reaction.crossSection : fissionComponent.crossSection.add( crossSection )
        reactionSuite.fissionComponents.add( fissionComponent )

    for production in productions:
        reactionSuite.productions.add( production )

    for incomplete in incompleteReactions:
        reactionSuite.incompleteReactions.add( incomplete )

    if hasattr( info, 'totalDelayedMultiplicity' ):
        prompt, delayed = [], []
        for neutron in reactionSuite.getReaction('fission').outputChannel.getProductsWithName( IDsPoPsModule.neutron ):
            link_ = sumsModule.add( link = neutron.multiplicity )
            if neutron.getAttribute('emissionMode') == tokensModule.delayedToken:
                delayed.append( link_ )
            else:
                prompt.append( link_ )

        delayedNubar = sumsModule.multiplicitySum( label = "delayed fission neutron multiplicity",
                ENDF_MT = 455, summands = sumsModule.listOfSummands(delayed) )
        delayedNubar.multiplicity.add( info.totalDelayedMultiplicity )
        reactionSuite.sums.multiplicities.add( delayedNubar )

        total = prompt + [sumsModule.add( link = delayedNubar.multiplicity )]
        totalNubar = sumsModule.multiplicitySum( label = "total fission neutron multiplicity",
                ENDF_MT = 452, summands = sumsModule.listOfSummands(total) )
        totalNubar.multiplicity.add( info.totalOrPromptFissionNeutrons['total'] )
        reactionSuite.sums.multiplicities.add( totalNubar )

    warningList = []
    try :               # Parse resonance section.
        mf2 = None
        if( 151 in MTDatas and not parseCrossSectionOnly ) :
            mf2 = MTDatas.get( 151 ).get( 2 )    # Resonance data available.
        if( mf2 ) :
            info.logs.write( '    Reading resonances (MF=2 MT=151)\n' )
            resonances, resonanceMTs = readMF2( info, mf2, warningList )
            if info.LRP == 2: # LRP was read in from first line of ENDF file
                if resonances.resolved: resonances.resolved.evaluated.useForSelfShieldingOnly = True
                if resonances.unresolved: resonances.unresolved.evaluated.useForSelfShieldingOnly = True
            reactionSuite.addResonances( resonances )

            if resonances.reconstructCrossSection:
                # modify cross sections for relevant channels to indicate resonance contribution is needed:
                resonanceLink = crossSectionModule.resonanceLink( link = resonances )

                for MT in resonanceMTs :
                    MTChannels  = [ r1 for r1 in reactionSuite.reactions         if( r1.ENDF_MT == MT ) ]
                    MTChannels += [ r1 for r1 in reactionSuite.sums.crossSections   if( r1.ENDF_MT == MT ) ]
                    MTChannels += [ r1 for r1 in reactionSuite.fissionComponents if( r1.ENDF_MT == MT ) ]
                    if( len( MTChannels ) == 0 ) :
                        if( MT in ( 3, 18, 19 ) ) :
                            continue
                        else :
                            warningList.append( 'Unable to find channel corresponding to resonance data for MT%d' % MT )
                    elif( len( MTChannels ) == 1 ) :
                        crossSectionComponent = MTChannels[0].crossSection
                        backgroundForm = crossSectionComponent[info.style]
                        backgroundForm.label = None
                        crossSectionComponent.remove( backgroundForm.label )

                        # break background up into resolved, unresolved and fast regions:
                        haveResolved = (resonances.resolved is not None and not
                                resonances.resolved.evaluated.useForSelfShieldingOnly)
                        haveUnresolved = (resonances.unresolved is not None and not
                                resonances.unresolved.evaluated.useForSelfShieldingOnly)

                        if isinstance(backgroundForm, crossSectionModule.XYs1d):
                            backgroundBoundaries = []
                        else:
                            backgroundBoundaries = [region.domainMax for region in backgroundForm[:-1]]

                        missingBoundaries = {}
                        if haveResolved and resonances.resolved.domainMax not in backgroundBoundaries:
                            missingBoundaries[resonances.resolved.domainMax] = ['mergeWithNextRegion']
                        if haveUnresolved and resonances.unresolved.domainMax not in backgroundBoundaries:
                            missingBoundaries[resonances.unresolved.domainMax] = ['mergeWithNextRegion']

                        if missingBoundaries:
                            # ENDF manual doesn't insist on breaking the background into regions, but GNDS does
                            bkForm2 = crossSectionModule.regions1d(axes=backgroundForm.axes)

                            domainMin = backgroundForm.domainMin
                            if isinstance(backgroundForm, crossSectionModule.XYs1d):
                                for boundary in sorted(missingBoundaries.keys()):
                                    if boundary not in backgroundForm.domainGrid:
                                        missingBoundaries[boundary].append('deletePoint')
                                    bkForm2.append( backgroundForm.domainSlice(domainMin, boundary) )
                                    domainMin = boundary
                                bkForm2.append( backgroundForm.domainSlice(domainMin, backgroundForm.domainMax) )

                            else:   # regions1d
                                startIdx = 0
                                for boundary in sorted(missingBoundaries.keys()):
                                    for idx in range(startIdx, len(backgroundForm)):
                                        if backgroundForm[idx].domainMax > boundary: break
                                        bkForm2.append( backgroundForm[idx] )
                                    # add the missing boundary:
                                    if boundary not in backgroundForm[idx].domainGrid:
                                        missingBoundaries[boundary].append('deletePoint')
                                    bkForm2.append(
                                        backgroundForm[idx].domainSlice(domainMin, boundary))
                                    domainMin = boundary
                                    startIdx = idx
                                bkForm2.append( backgroundForm[idx].domainSlice(domainMin, backgroundForm[idx].domainMax) )

                                for region in backgroundForm[idx + 1:]:
                                    bkForm2.append(region)

                            for link in linksToCheck:   # update links if necessary:
                                if link.link is backgroundForm:
                                    link.link = bkForm2
                            backgroundForm = bkForm2

                        RRBack = URRBack = fastBack = None

                        idx = 0
                        if haveResolved:
                            while backgroundForm[idx].domainMax < resonances.resolved.domainMax:
                                idx += 1
                            if backgroundForm[idx].domainMax > resonances.resolved.domainMax:
                                warningList.append("Domain mismatch between resolved resonances and background")
                                info.doRaise.append(warningList[-1])
                            if idx == 0:
                                RRBack = backgroundForm[idx]
                                RRBack.index = None
                            else:
                                RRBack = crossSectionModule.regions1d(axes=backgroundForm.axes)
                                for jdx in range(idx+1):
                                    RRBack.append(backgroundForm[jdx])
                            RRBack = crossSectionModule.resolvedRegion( RRBack )
                            if resonances.resolved.domainMax in missingBoundaries:
                                info.ENDFconversionFlags.add( RRBack, ",".join(missingBoundaries[resonances.resolved.domainMax]) )
                            idx += 1

                        startIdx = idx
                        if haveUnresolved:
                            while backgroundForm[idx].domainMax < resonances.unresolved.domainMax:
                                idx += 1
                            if backgroundForm[idx].domainMax != resonances.unresolved.domainMax:
                                warningList.append("Domain mismatch between unresolved resonances and background")
                                info.doRaise.append(warningList[-1])
                            if idx==startIdx:
                                URRBack = backgroundForm[idx]
                                URRBack.index = None
                            else:
                                URRBack = crossSectionModule.regions1d(axes=backgroundForm.axes)
                                for jdx in range(startIdx,idx+1):
                                    URRBack.append(backgroundForm[jdx])
                            URRBack = crossSectionModule.unresolvedRegion( URRBack )
                            if resonances.unresolved.domainMax in missingBoundaries:
                                info.ENDFconversionFlags.add( URRBack, ",".join(missingBoundaries[resonances.unresolved.domainMax]) )
                            idx += 1

                        startIdx = idx
                        if startIdx == len(backgroundForm)-1:
                            fastBack = backgroundForm[startIdx]
                            fastBack.index = None
                        else:
                            fastBack = crossSectionModule.regions1d(axes=backgroundForm.axes)
                            for jdx in range(startIdx, len(backgroundForm)):
                                fastBack.append(backgroundForm[jdx])
                        fastBack = crossSectionModule.fastRegion( fastBack )

                        background_ = crossSectionModule.background( RRBack, URRBack, fastBack )
                        crossSectionComponent.add( crossSectionModule.resonancesWithBackground(
                            info.style, resonanceLink, background_, backgroundForm.uncertainty ) )
                        for link in linksToCheck:
                            if link.link is backgroundForm:
                                link.link = crossSectionComponent[ info.style ]
                    else :
                        raise NotImplementedError("Multiple reactions match resonance MT%d" % MT)

    except BadResonances as e:
        warningList.append( '       ERROR: unable to parse resonances! Error message: %s' % e )
        info.doRaise.append( warningList[-1] )

    MF12BaseMTsAndRange = [ [ 50, 92 ], [ 600, 650 ], [ 650, 700 ], [ 700, 750 ], [ 750, 800 ], [ 800, 850 ] ]

    if( singleMTOnly is None ) :
        for MTLO2, MF12_LO2 in sorted(info.MF12_LO2.items()) :  # The logic below assumes MTs are in ascending order per base MT.
            branchingBaseMT = None
            for MTBase, MTEnd in MF12BaseMTsAndRange :             # Determine base MT for this MTLO2
                if( MTBase < MTLO2 < MTEnd ) :
                    branchingBaseMT = MTBase
                    break
            if( branchingBaseMT is not None ) :
                residualZA = endf_endlModule.ENDF_MTZAEquation( info.projectileZA, info.targetZA, branchingBaseMT )[0][-1]
                residual = toGNDSMiscModule.getTypeNameENDF( info, residualZA, None )
                residualName = residual.id
                level = MTLO2 - branchingBaseMT
                levelEnergy = MF12_LO2[0]['ES']
                fullName = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex( residualName, level )

                decayProduct = None
                for reaction in reactionSuite.reactions :
                    if( MTLO2 == reaction.ENDF_MT ) : 
                        if( isinstance( reaction.outputChannel, channelsModule.twoBodyOutputChannel ) ) :
                            _decayProduct = reaction.outputChannel[1]
                            if( _decayProduct.id == fullName ) :
                                if( _decayProduct.outputChannel is None ) :
                                    decayProduct = _decayProduct
                                    break
                if( decayProduct is None ) :
                    pass
                else :
                    crossSection = reaction.crossSection[0]
                    decayChannel = channelsModule.NBodyOutputChannel( )
                    decayChannel.Q.add( toGNDSMiscModule.returnConstantQ( info.style, -reaction.outputChannel.Q[0].constant, crossSection ) )

                    decayResidualName = reactionSuite.PoPs[fullName].isotope.symbol
                    decayChannel.products.add( toGNDSMiscModule.newGNDSParticle( info, decayResidualName, crossSection ) )

                    multiplicity = multiplicityModule.branching1d( info.style )
                    decayPhoton = toGNDSMiscModule.newGNDSParticle( info, IDsPoPsModule.photon, crossSection, multiplicity = multiplicity )
                    pids = branching3dModule.pids( decayProduct.id, decayResidualName )
                    decayPhoton.distribution.add( branching3dModule.form( info.style, pids, standardsModule.frames.labToken ) )
                    decayChannel.products.add( decayPhoton )

                    decayProduct.addOutputChannel( decayChannel )
                    reaction.updateLabel( )

                particleLevelEnergy_eV = reactionSuite.PoPs[fullName].energy[0].value                   # Compare this value to level energy from the particle list (from MF3 Q-value).
                if( levelEnergy != particleLevelEnergy_eV ) :
                    if( particleLevelEnergy_eV < 1e-12 ) :
                        warningList.append( "MF12 parent level energy (%s) set to zero?" % particleLevelEnergy_eV )
                        info.doRaise.append( warningList[-1] )
                    elif( abs( levelEnergy - particleLevelEnergy_eV ) < 1e-4 * particleLevelEnergy_eV ) :
                        MFLabel = '3'
                                                                                            # Value with most precision wins.
                        str1 = PQUModule.floatToShortestString( levelEnergy * 1e-20 )          # 1e-20 to insure e-form is returned.
                        str2 = PQUModule.floatToShortestString( particleLevelEnergy_eV * 1e-20 )  # Want 1.23e-16 and not 12300 to differ
                        if( len( str1 ) > len( str2 ) ) :                                   # length from 1.2345e-16 and not 12345.
                            reactionSuite.PoPs[fullName].energy[0].value = levelEnergy
                            MFLabel = '12'
                        if( str1 != str2 ) :
                            warningList.append( "MT%d MF12 level energy %s differs from MF3 value %s. Setting to MF%s value." %
                                    ( MTLO2, levelEnergy, particleLevelEnergy_eV, MFLabel ) )
                    else :
                        warningList.append( "MT%d MF12 parent level energy (%s) doesn't match known level" % ( MTLO2, particleLevelEnergy_eV ) )
                        info.doRaise.append( warningList[-1] )
                for i1, MF12 in enumerate( MF12_LO2 ) :
                    try :
                        finalLevelEnergy = MF12['ESk']
                        if( finalLevelEnergy > 0. ) :   # Find particle in the particleList with energy = finalLevelEnergy
                            finalParticles = [ lev for lev in reactionSuite.PoPs[residualName].ancestor
                                    if lev.energy.float('eV') == finalLevelEnergy ]
                            if( len( finalParticles ) == 1 ) :
                                finalParticle = finalParticles[0]
                            else :                      # No exact match, look for levels within .01% of the exact value.
                                idx = 0
                                while( True ) :
                                    idx += 1
                                    finalParticleName = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex( residualName, idx )
                                    if( not reactionSuite.hasParticle( finalParticleName ) ) :
                                        warningList.append( "MF12 final level energy (%s eV) doesn't match known level when decaying out of level %s " % \
                                                ( finalLevelEnergy, MTLO2 ) )
                                        info.doRaise.append( warningList[-1] )
                                    try :
                                        thisLevelEnergy = reactionSuite.PoPs[finalParticleName].energy.pqu( ).getValueAs( 'eV' )
                                    except KeyError :
                                        raise Exception( 'Could not find nuclide of %s with desired energy level of %s.' % ( residualName, finalLevelEnergy ) )
                                    except :
                                        raise
                                    if( abs( thisLevelEnergy - finalLevelEnergy ) < 1e-4 * finalLevelEnergy ) :
                                        finalParticle = reactionSuite.PoPs[finalParticleName]
                                        break   # found it
                        else :
                            finalParticle = reactionSuite.PoPs[residualName]
                        gammaTransition = 1.
                        if( len( MF12['branching'] ) > 2 ) : gammaTransition = MF12['branching'][1]

                        if( gammaTransition != 1 ) : raise Exception( 'Fix me' )
                        probability = probabilityModule.double( info.PoPsLabel, MF12['branching'][0] )

                        decayMode = decayDataModule.decayMode( str( i1 ), 'electroMagnetic' )
                        decayMode.probability.add( probability )
                        _decay = decayDataModule.decay( str( i1 ), decayDataModule.decayModesParticle)
                        _decay.products.add( productModule.product( IDsPoPsModule.photon, IDsPoPsModule.photon ) )
                        _decay.products.add( productModule.product( finalParticle.id, finalParticle.id ) )
                        decayMode.decayPath.add( _decay )
                        if MF12['LG'] == 2: # internal conversion competes with gamma emission
                            from PoPs.decays import spectrum as spectrumModule
                            Pgamma = spectrumModule.shell( MF12['branching'][1], label=spectrumModule.shell.total )
                            decayMode.photonEmissionProbabilities.add( Pgamma )
                        reactionSuite.PoPs[fullName].decayData.decayModes.add( decayMode )
                    except Exception as err :
                        warningList.append( 'raise somewhere in "for MF12 in MF12_LO2" loop: MT%d, %s' % ( MT, str( err ) ) )
                        info.doRaise.append( warningList[-1] )
            else :
                raise Exception( "Could not determine base MT for MF=12's MT=%s" % MTLO2 )

    if( doCovariances ) :
        covarianceMFs = sorted( set( [mf for mt in MTDatas.values() for mf in mt.keys() if mf>30] ) )
        if covarianceMFs:
            info.logs.write( '    Reading covariances (MFs %s)\n' % ','.join(map(str,covarianceMFs) ) )
        try:
            """ parse covariances. This also requires setting up links from data to covariances, so we
            must ensure the links are synchronized """

            MTdict = {}
            for reaction in ( list( reactionSuite.reactions ) + list( reactionSuite.sums.crossSections ) + list( reactionSuite.productions )
                    + list( reactionSuite.fissionComponents ) + list( reactionSuite.incompleteReactions ) ) :
                MT = reaction.ENDF_MT
                if MT in MTdict:
                    MTdict[MT].append( reaction )
                else:
                    MTdict[MT] = [reaction]
            covarianceSuite, linkData = parseCovariances( info, MTDatas, MTdict, singleMTOnly = singleMTOnly,
                    resonances = getattr( reactionSuite, 'resonances', None ) )
            if( len( covarianceSuite.covarianceSections ) > 0 or len( covarianceSuite.parameterCovariances ) > 0 ) :
                covarianceSuite.target = str(info.target)
                covarianceSuite.projectile = str(info.projectile)
                covarianceSuite.styles.add( info.evaluatedStyle )
                reactionSuite.externalFiles.add( externalFileModule.externalFile("covariances","test.endf6-covar.xml") )
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

    info.massTracker.useMostCommonAMUmasses()

    if( info.level > 0 ) : # AWR is for isomer mass. Adjust info.ZAMasses to GS mass:
        groundStateMass = info.massTracker.getMassAMU( info.targetZA ) - PQUModule.PQU(
            PQUModule.pqu_float.surmiseSignificantDigits( info.level ),'eV/c**2').getValueAs('amu')
        info.massTracker.addMassAMU( info.targetZA, groundStateMass )  # overwrite excited state mass

    for ZA in info.massTracker.amuMasses :
        if( ZA in [ 1 ] ) : continue
        mass = info.massTracker.amuMasses[ZA]
        elementSymbol = chemicalElementMiscModule.symbolFromZ[ZA/1000]
        name = chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( elementSymbol, ZA % 1000 )
        name = chemicalElementMiscModule.nuclideIDFromIsotopeSymbolAndIndex( name, 0 )
        mass = massModule.double( info.PoPsLabel, mass, quantityModule.stringToPhysicalUnit( 'amu' ) )
        if( name not in reactionSuite.PoPs ) : toGNDSMiscModule.getPoPsParticle( info, ZA, levelIndex = 0 )
        particle = reactionSuite.PoPs[name]
        particle.mass.add( mass )

    sys.stdout.flush( )
    for warning in warningList : info.logs.write( "       WARNING: %s\n" % warning, stderrWriting = True )

    targetID = reactionSuite.target
    if( targetID in reactionSuite.PoPs.aliases ) : targetID = reactionSuite.PoPs[targetID].pid
    ignoreID = None
    for particleID, spinParity in info.particleSpins.items( ) :
        if( particleID == 'target' ) : ignoreID = targetID
    for particleID, spinParity in info.particleSpins.items( ) :
        if( ignoreID == particleID ) : continue
        spin = spinModule.fraction( info.PoPsLabel, fractions.Fraction( spinParity[0] ), spinModule.baseUnit )
        if( particleID == reactionSuite.target ) : particleID = targetID
        if( particleID == 'target' ) :
            particle = reactionSuite.PoPs[targetID]
        else :
            particle = reactionSuite.PoPs[particleID]

        if( isinstance( particle, nuclideModule.particle ) ) : particle = particle.nucleus

        if( len( particle.spin ) == 0 ) : particle.spin.add( spin )
        if( spinParity[1] ) :
            parity = spinParity[1].value
            particle.parity.add( parityModule.integer( info.PoPsLabel, parity, parityModule.baseUnit ) )

    for reaction in reactionSuite.reactions :
        if( isinstance( reaction.outputChannel, channelsModule.twoBodyOutputChannel ) ) :
            residual = reaction.outputChannel[1]
            if( residual.outputChannel is None ) :
                particle = info.PoPs[residual.id]
                if( isinstance( particle, nuclideModule.particle ) ) :
                    if( ( particle.nucleus.index != 0 ) and ( len( particle.decayData.decayModes ) == 0 ) ) :   # gamma data should be in orphanProducts.
                        groundState = particle.isotope.nuclides[0]
                        ZA = chemicalElementMiscModule.ZA( groundState )
                        multiplicity = residual.multiplicity[0].copy( )
                        residual.addOutputChannel( channelsModule.NBodyOutputChannel( ) )
                        residual.outputChannel.Q.add( toGNDSMiscModule.returnConstantQ( info.style, particle.nucleus.energy[0].value, multiplicity ) )
                        product = toGNDSMiscModule.newGNDSParticle( info, toGNDSMiscModule.getPoPsParticle( info, ZA, groundState.id ),
                                    crossSection, multiplicity = multiplicity )
                        product.distribution.add( unspecifiedModule.form( info.style, standardsModule.frames.labToken ) )
                        info.ENDFconversionFlags.add( product, 'implicitProduct' )
                        residual.outputChannel.products.add( product )
                        reaction.updateLabel( )

    if( reconstructResonances and reactionSuite.resonances is not None and reactionSuite.resonances.reconstructCrossSection ):
        info.logs.write( '    Reconstructing resonances\n' )
        try:
            reactionSuite.reconstructResonances( style=info.reconstructedStyle,
                accuracy=info.reconstructedAccuracy, verbose = verbose, thin=True )
        except Exception as ex:
            warningList.append( "Resonance reconstruction failed: %s" % ex )
            info.doRaise.append( warningList[-1] )

    def adjustMF13Multiplicity2( multiplicity, crossSection ) :

        energyMultiplicity = []
        if( multiplicity.domainMax > crossSection.domainMax ) :
                multiplicity = multiplicity.domainSlice( domainMax = crossSection.domainMax )
        for energyIn, multiplicityValue in multiplicity :
            crossSectionAtEnergy = crossSection.evaluate( energyIn )
            if( crossSectionAtEnergy != 0 ) : multiplicityValue /= crossSectionAtEnergy
            energyMultiplicity.append( [ energyIn, multiplicityValue ] )
        multiplicity.setData( energyMultiplicity )

    def adjustMF13Multiplicity( multiplicity, crossSection ) :

        if( isinstance( multiplicity, multiplicityModule.XYs1d ) ) :
            adjustMF13Multiplicity2( multiplicity, crossSection )
        elif( isinstance( multiplicity, multiplicityModule.regions1d ) ) :
            for region in multiplicity : adjustMF13Multiplicity2( region, crossSection )
        else :
            raise Exception( 'Unsupported multiplicity type "%s"' % multiplicity.moniker )

    def adjustMF13Gammas( reaction ) :

        MT = reaction.ENDF_MT
        crossSection = None
        allproducts = list( reaction.outputChannel )
        for prod in reaction.outputChannel :
            if prod.outputChannel is not None :
                allproducts.extend( list( prod.outputChannel ) )

        for product in allproducts :
            multiplicity = product.multiplicity[info.style]
            if( hasattr( multiplicity, '_temp_divideByCrossSection' ) ) :
                if( crossSection is None ) :
                    crossSection = reaction.crossSection.toPointwise_withLinearXYs( accuracy = 1e-3, upperEps = 1e-8 )
                adjustMF13Multiplicity( multiplicity, crossSection )
                del multiplicity._temp_divideByCrossSection

        if( MT in info.totalMF6_12_13Gammas ) :
            MF, multiplicity = info.totalMF6_12_13Gammas[MT]
            if( MF == 13 ) : adjustMF13Multiplicity( multiplicity, crossSection )
            gammaProduction = [ tmp for tmp in reactionSuite.reactions if tmp.ENDF_MT == MT ] # raises ValueError if more than one match found
            gammaProduction += [ tmp for tmp in reactionSuite.orphanProducts if tmp.ENDF_MT == MT ]
            if( len( gammaProduction ) != 1 ) : raise ValueError( "No unique match found." )
            gammaProduction = gammaProduction[0]
            summands = [ sumsModule.add( link = r.multiplicity ) for r in gammaProduction.outputChannel.getProductsWithName( IDsPoPsModule.photon ) ]
            if len(summands)==0:
                for _product in gammaProduction.outputChannel:
                    if _product.outputChannel is not None:
                        summands += [ sumsModule.add( link = r.multiplicity ) for r in _product.outputChannel.getProductsWithName( IDsPoPsModule.photon ) ]
            if( MT in channelIDs ) :
                name = channelIDs[MT]
            else :
                name = str( gammaProduction.outputChannel )

            multiplicitySum = sumsModule.multiplicitySum( label = name + " total gamma multiplicity",
                    ENDF_MT = MT, summands = sumsModule.listOfSummands( summands ) )
            multiplicitySum.multiplicity.add( multiplicity )
            reactionSuite.sums.multiplicities.add( multiplicitySum )

    for reaction in reactionSuite.reactions : adjustMF13Gammas( reaction )
    for reaction in reactionSuite.orphanProducts : adjustMF13Gammas( reaction )

    # Short-lived light isotopes sometimes listed implicitly in ENDF-6, without associated mass:
    for pid,ZA in (('He5',2005), ('Li5',3005), ('Be8',4008)):
        if pid in reactionSuite.PoPs:
            particle = reactionSuite.PoPs[pid]
            try:
                particle.getMass('amu')
            except:
                from fudge.legacy.endl.structure import masses as AMEmasses
                particle.mass.add( massModule.double(info.PoPsLabel, AMEmasses.getMassFromZA(ZA), 'amu'))

    return( covarianceSuite )
