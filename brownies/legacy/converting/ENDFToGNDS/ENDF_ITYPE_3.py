# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
For translating ENDF ITYPE=3 data (photo-atomic and electro-atomic interactions)
"""
#   500 Total charged-particle stopping power       # Currently not supported (executes a raise).
#   501 total photon interaction                    # For gamma, sum of 502, 504, 516, 522
#   502 photon coherent scattering                  # Used with MT 505 and 506
#   504 photon incoherent scattering
#   505 imaginary scattering factor
#   506 real scattering factor
#   515 pair production, electron field
#   516 pair production                             # Sum of MT 515 and 517
#   517 pair production, nuclear field
#   522 photoelectric absorption                    # For gamma, sum of 534 to 572
#   523 photo-excition cross section
#   526 electro-atomic scattering                   # Large angle Coulomb scattering
#   527 electro-atomic bremsstrahlung
#   528 electro-atomic excitation cross section
#   534 to 572 subshell photoelectric- or electro-atomic cross section

# 2-body 
#   T(g,g)T    502, 505, 506 and 504
#   T(g,e-)T{}        522           For gamma, sum of 534 to 572
#   T(e-,e-)T{?}      526           called 'electro-atomic elastic scattering' on page 198 and 'electro-atomic scattering' on page 315
#
# Uncorrelated bodies
#   T(g,e- e+)T       515, 516, 517
#   T(g,)T{}          523
#   T(e-,e- g)T       527
#   T(e-,e-?)T{}      528

from PoPs import IDs as IDsPoPsModule
from PoPs.groups import misc as chemicalElementMiscPoPsModule

import xData.standards as standardsModule
import brownies.legacy.converting.toGNDSMisc as toGNDSMiscModule

from fudge import outputChannel as outputChannelModule
import fudge.sums as sumsModule

import fudge.reactions.reaction as reactionModule
import fudge.reactions.incompleteReaction as incompleteReactionModule

from fudge.reactionData.doubleDifferentialCrossSection import base as baseModule
from fudge.reactionData.doubleDifferentialCrossSection.photonScattering import incoherent as incoherentModule
from fudge.reactionData.doubleDifferentialCrossSection.photonScattering import coherent as coherentModule

from fudge.productData import multiplicity as multiplicityModule
from fudge.productData.distributions import photonScattering as photonScatteringModule
from fudge.productData.distributions import unspecified as unspecifiedModule

from . import endfFileToGNDSMisc
from .ENDF_ITYPE_0_Misc import readMF3, readMF6
from .ENDF_ITYPE_3_6_Misc import MT_AtomicConfigurations
from .. import endf_endl as endf_endlModule

photonSumLabels = { 501 : 'total', 516 : 'pair product', 522 : 'photoelectric absorption' }

def readMF27( info, MT, MF27Datas, label, warningList ) :

    _class = { 502 : coherentModule.formFactor, 
               504 : None,
               505 : coherentModule.imaginaryAnomalousFactor,
               506 : coherentModule.realAnomalousFactor }[MT]

    XYs1d = coherentModule.XYs1d
    if( MT == 504 ) : XYs1d = incoherentModule.XYs1d

    MF27Data = MF27Datas[27]
    del MF27Datas[27]
    energyUnit = 'eV'
    if( MT in [ 502, 504 ] ) : energyUnit = '1/Ang'
    axes = baseModule.defaultAxes( label, energyUnit = energyUnit )
    dataLine, TAB1, MF27s = endfFileToGNDSMisc.getTAB1Regions(1, MF27Data, axes = axes, allowInterpolation6 = True, logFile = info.logs, cls = XYs1d)
    if( len( MF27s ) == 1 ) :
        if( MT == 504 ) : return( MF27s[0] )
        data = coherentModule.XYs1d( data = MF27s[0], axes = axes )
    else :
        if( MT == 504 ) :
            data = incoherentModule.regions1d( axes = axes )
        else :
            data = coherentModule.regions1d( axes = axes )
        for region in MF27s :
            if( len( region ) > 1 ) : data.append( region )
        if( MT == 504 ) : return( data )
    return( _class( data ) )

def extractMT505or506( info, MT, MTList, MTDatas, label ) :

    if( MT in MTDatas ) :
        warningList = []
        MFDatas = MTDatas[MT]
        factor = readMF27( info, MT, MFDatas, label, warningList )
        if( len( MFDatas ) != 0 ) : raise Exception( 'Unsupported MF values = %s' % MFDatas.keys( ) )
        MTList.remove( MT )
        for warning in warningList : info.logs.write( "       WARNING: %s\n" % warning, stderrWriting = True )
        return( factor )

    return( None )

def ITYPE_3( MTDatas, info, reactionSuite, singleMTOnly, parseCrossSectionOnly, verbose = 0 ) :

    MTs = sorted( MTDatas.keys( ) )
    if( 451 in MTs ) : MTs.remove( 451 )
    MTList = endfFileToGNDSMisc.niceSortOfMTs(MTs, verbose = verbose, logFile = info.logs)

    MT505 = extractMT505or506( info, 505, MTList, MTDatas, 'imaginary part of anomalous scattering factor' )
    MT506 = extractMT505or506( info, 506, MTList, MTDatas, 'real part of anomalous scattering factor' )

    iChannel = 0
    summedReactions = { 501 : None, 516 : None, 522 : None }
    for MT in MTList :
        if( MT == 500 ) : raise Exception( "MT %d not supported" % MT )
        if( ( singleMTOnly is not None ) and ( MT != singleMTOnly ) ) : continue

        if( MT in [ 523 ] ) :
            print('    Skipping MT %d as I do not know what it really is' % MT)
            continue

        warningList = []
        MTData = MTDatas[MT]

        info.logs.write( '    %3d %s' % ( MT, sorted( MTData.keys( ) ) ) )
        EPE, EFL, crossSection, LR, breakupProducts = readMF3( info, MT, MTData[23], warningList )
        if EFL != 0:
            raise NotImplementedError("Non-zero EFL flag in incident electron evaluation")
        del MTData[23]

        isTwoBody = MT in (525, 526)
        productList = []
        doubleDifferentialCrossSectionForm = None
        undefinedLevelInfo = { 'ZA' : None, 'level' : 0, 'levelIndex' : None, 'count' : 0 }

        if( 26 in MTData ) :
            isTwoBody = readMF6( MT, info, MTData[26], productList, warningList, undefinedLevelInfo, isTwoBody, crossSection, 0 )
            del MTData[26]
            if( MT == 527 ) :
                for product in productList :
                    if( isinstance( product.multiplicity[0], multiplicityModule.XYs1d ) ) :
                        multiplicity = product.multiplicity.pop( product.multiplicity[0].label )
                        multiplicity = multiplicityModule.constant1d( 1, multiplicity.domainMin, multiplicity.domainMax, 
                                axes = multiplicity.axes, label = multiplicity.label )
                        product.multiplicity.add( multiplicity )

        addTargetAsResidual = True
        productsNeeded = []
        try : 
            process = { 502 : 'coherent',                           504 : 'incoherent',                     522 : None,
                        515 : 'pair production: electron field',    517 : 'pair production: nuclear field', 525 : 'large angle electro-atomic scattering',
                        526 : 'large angle Coulomb scattering',     527 : 'bremsstrahlung',                 528 : 'excitation' }[MT]
        except KeyError:
            process = None
        if( MT in [ 502, 504, 522 ] ) :
            outputChannel = outputChannelModule.outputChannel( outputChannelModule.Genre.twoBody, process = process )
            outputChannel.Q.add( toGNDSMiscModule.returnConstantQ( info.style, 0, crossSection ) )
            if( MT in [ 502 , 504 ] ) :
                product = toGNDSMiscModule.getTypeNameGamma( info, 0 )
                product = toGNDSMiscModule.newGNDSParticle( info, product, crossSection )
                outputChannel.products.add( outputChannel.products.uniqueLabel( product ) )
                if( MT == 502 ) :
                    formFactor = readMF27( info, MT, MTData, 'coherent form factor', warningList )
                    doubleDifferentialCrossSectionForm = coherentModule.form( product.id, info.style, standardsModule.frames.labToken, 
                            formFactor, MT506, MT505 )
                    form = photonScatteringModule.coherentPhotonScattering.form( label = info.style, link = doubleDifferentialCrossSectionForm )
                else :
                    subform = readMF27( info, MT, MTData, 'incoherent scattering function', warningList )
                    doubleDifferentialCrossSectionForm = incoherentModule.form( product.id, info.style, standardsModule.frames.labToken, subform )
                    form = photonScatteringModule.incoherentPhotonScattering.form( label = info.style, link = doubleDifferentialCrossSectionForm )
                product.distribution.add( form )
        elif( MT in [ 525, 526, 527, 528 ] ) :
            if( isTwoBody ) :
                outputChannel = outputChannelModule.outputChannel( outputChannelModule.Genre.twoBody, process = process )
            else :
                outputChannel = outputChannelModule.outputChannel( outputChannelModule.Genre.NBody, process = process )
            outputChannel.Q.add( toGNDSMiscModule.returnConstantQ( info.style, 0, crossSection ) )
            productsNeeded = [ IDsPoPsModule.electron ]
            if( MT == 527 ) : productsNeeded.insert( 0, IDsPoPsModule.photon )
        else :
            if( MT in [ 515, 517 ] ) : EPE = -crossSection.domainMin
            outputChannel = outputChannelModule.outputChannel( outputChannelModule.Genre.NBody, process = process )
            outputChannel.Q.add( toGNDSMiscModule.returnConstantQ( info.style, EPE, crossSection ) )
            if( MT in MT_AtomicConfigurations ) :
                productsNeeded = [ IDsPoPsModule.electron ]
            elif( MT == 515 ) :
                productsNeeded = [ IDsPoPsModule.electron, IDsPoPsModule.electron, IDsPoPsModule.electronAnti ]
            elif( MT == 517 ) :
                productsNeeded = [ IDsPoPsModule.electron, IDsPoPsModule.electronAnti ]
            else :
                addTargetAsResidual = False

        for pid in productsNeeded :
            product = None
            for i1, product in enumerate( productList ) :
                if( pid == product.pid ) : break
                product = None
            if( product is None ) :
                product = toGNDSMiscModule.getTypeNameGamma( info, { IDsPoPsModule.photon : 0, IDsPoPsModule.electronAnti : 8, IDsPoPsModule.electron : 9 }[pid] )
                product = toGNDSMiscModule.newGNDSParticle( info, product, crossSection )
            outputChannel.products.add( outputChannel.products.uniqueLabel( product ) )
        if( ( MT >= 534 ) and ( reactionSuite.projectile == IDsPoPsModule.electron ) ) :
            product = toGNDSMiscModule.getTypeNameGamma( info, 9 )
            product = toGNDSMiscModule.newGNDSParticle( info, product, crossSection )
            outputChannel.products.add( outputChannel.products.uniqueLabel( product ) )

        if( addTargetAsResidual ) :
            elementSymbol = chemicalElementMiscPoPsModule.symbolFromZ[info.targetZA//1000]
            if( MT >= 534 ) : elementSymbol += '{%s}' % MT_AtomicConfigurations[MT]
            residual = toGNDSMiscModule.newGNDSParticle( info, elementSymbol, crossSection )
            residual.distribution.add(
                unspecifiedModule.form( info.style, productFrame = standardsModule.frames.labToken )
            )
            outputChannel.products.add( outputChannel.products.uniqueLabel( residual ) )
            info.ENDFconversionFlags.add( residual, 'implicitProduct' )

        if( len( MTData ) > 0 ) : raise Exception( 'Untranslated MF data: MFs = %s' % MTData.keys( ) )

        if( MT in summedReactions ) :
            summedReactions[MT] = [ crossSection, outputChannel ]
            info.logs.write( '\n' )
            continue

        if( MT == 525 ) :
            reaction = incompleteReactionModule.incompleteReaction( outputChannel.genre, ENDF_MT = MT )
        else :
            reaction = reactionModule.reaction( outputChannel.genre, ENDF_MT = MT )
        endf_endlModule.setReactionsOutputChannelFromOutputChannel( info, reaction, outputChannel )
        reaction.crossSection.add( crossSection )
        if( doubleDifferentialCrossSectionForm is not None ) : reaction.doubleDifferentialCrossSection.add( doubleDifferentialCrossSectionForm )

        if( MT == 525 ) :
            reactionSuite.incompleteReactions.add( reaction )
        else :
            reactionSuite.reactions.add( reaction )

        info.logs.write( '\n' )
        for warning in warningList : info.logs.write( "       WARNING: %s\n" % warning, stderrWriting = True )

    if( True ) :            # Need to test for gamma as projectile.
        for MT in sorted( summedReactions.keys( ) ) :
            if( summedReactions[MT] is None ) : continue
            crossSection, outputChannel = summedReactions[MT]
            summandMTs = list( { 501 : range( 502, 573 ), 516 : ( 515, 517 ), 522 : range( 534, 573 ) }[MT] )
            if( 525 in summandMTs ) : summandMTs.remove( 525 )

            summedCrossSection = sumsModule.crossSectionSum( label=photonSumLabels[MT], ENDF_MT=MT )
            summedCrossSection.crossSection.add( crossSection )
            for reaction in reactionSuite :
                if( reaction.ENDF_MT in summandMTs ) : summedCrossSection.summands.append( sumsModule.add( link = reaction.crossSection ) )
            Q = outputChannel.Q[info.style]
            if( Q is not None ) : summedCrossSection.Q.add( Q )
            reactionSuite.sums.crossSectionSums.add( summedCrossSection )
            iChannel += 1
