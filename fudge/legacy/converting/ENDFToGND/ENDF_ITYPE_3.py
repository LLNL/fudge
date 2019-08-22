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

from fudge.core.utilities import brb

import pqu.PQU as PQUModule
import xData.standards as standardsModule
import fudge.legacy.converting.toGNDMisc as toGNDMiscModule

import fudge.particles.nuclear as nuclearModule
import fudge.gnd.xParticle as xParticleModule

import fudge.gnd.channels as channelsModule
import fudge.gnd.sums as sumsModule

import fudge.gnd.reactions.reaction as reactionModule

import fudge.gnd.productData.distributions.photonScattering as photonScatteringModule

import endfFileToGNDMisc
from ENDF_ITYPE_0_Misc import readMF3, readMF6
from ENDF_ITYPE_3_6_Misc import MT_AtomicConfigurations

photonSumLabels = { 501 : 'total', 516 : 'pair product', 522 : 'photoelectric absorption' }

def readMF27( info, MT, MF27Datas, label, warningList ) :

    _class = { 502 : photonScatteringModule.scatteringFunction, 
               504 : None,
               505 : photonScatteringModule.imaginaryAnomalousFactor,
               506 : photonScatteringModule.realAnomalousFactor }[MT]

    MF27Data = MF27Datas[27]
    del MF27Datas[27]
    energyUnit = 'eV'
    if( MT in [ 502, 504 ] ) : energyUnit = '1/Ang'
    axes = photonScatteringModule.scatteringFactor.defaultAxes( label, energyUnit = energyUnit )
    dataLine, TAB1, MF27s = endfFileToGNDMisc.getTAB1Regions( 1, MF27Data, axes = axes, allowInterpolation6 = True, logFile = info.logs )
    if( len( MF27s ) == 1 ) :
        if( MT == 504 ) : return( photonScatteringModule.incoherent.XYs1d( data = MF27s[0], axes = axes, accuracy = 1e-3 ) )
        data = photonScatteringModule.scatteringFactor.XYs1d( data = MF27s[0], axes = axes, accuracy = 1e-3 )
    else :
        if( MT == 504 ) :
            data = photonScatteringModule.incoherent.regions1d( axes = axes )
        else :
            data = photonScatteringModule.scatteringFactor.regions1d( axes = axes )
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

def ITYPE_3( MTDatas, info, reactionSuite, singleMTOnly, parseCrossSectionOnly, verbose = False ) :

    MTs = sorted( MTDatas.keys( ) )
    if( 451 in MTs ) : MTs.remove( 451 )
    MTList = endfFileToGNDMisc.niceSortOfMTs( MTs, verbose = verbose, logFile = info.logs )

    MT505 = extractMT505or506( info, 505, MTList, MTDatas, 'imaginary part of anomalous scattering factor' )
    MT506 = extractMT505or506( info, 506, MTList, MTDatas, 'real part of anomalous scattering factor' )

    iChannel = 0
    summedReactions = { 501 : None, 516 : None, 522 : None }
    for MT in MTList :
        if( MT == 500 ) : raise Exception( "MT %d not supported" % MT )
        if( ( singleMTOnly is not None ) and ( MT != singleMTOnly ) ) : continue

        if( MT in [ 523 ] ) :
            print '    Skipping MT %d as I do not know what it really is' % MT
            continue

        warningList = []
        MTData = MTDatas[MT]

        info.logs.write( '    %3d %s' % ( MT, sorted( MTData.keys( ) ) ) )
        EPE, EFL, crossSection, breakupProducts = readMF3( info, MT, MTData[23], warningList )
        del MTData[23]

        isTwoBody, productList = MT == 526, []
        undefinedLevelInfo = { 'ZA' : None, 'level' : 0, 'levelIndex' : None, 'count' : 0 }
        if( 26 in MTData ) :
            isTwoBody = readMF6( MT, info, MTData[26], productList, warningList, undefinedLevelInfo, isTwoBody )
            del MTData[26]

        addTargetAsResidual = True
        productsNeeded = []
        try : 
            process = { 502 : 'coherent',                       504 : 'incoherent',     522 : None,
                        515 : 'electron field',                 517 : 'nuclear field',
                        526 : 'large angle Coulomb scattering', 527 : 'bremsstrahlung', 528 : 'excitation' }[MT]
        except KeyError:
            process = None
        if( MT in [ 502, 504, 522 ] ) :
            outputChannel = channelsModule.twoBodyOutputChannel( process = process )
            outputChannel.Q.add( toGNDMiscModule.returnConstantQ( info.style, 0 ) )
            if( MT in [ 502 , 504 ] ) :
                product = toGNDMiscModule.getTypeNameGamma( info, 0 )
                product = toGNDMiscModule.newGNDParticle( info, product )
                outputChannel.products.add( outputChannel.products.uniqueLabel( product ) )
                if( MT == 502 ) :
                    formFactor = readMF27( info, MT, MTData, 'coherent scattering function', warningList )
                    form = photonScatteringModule.coherent.form( info.style, standardsModule.frames.labToken, 
                            formFactor, MT506, MT505 )
                else :
                    subform = readMF27( info, MT, MTData, 'incoherent scattering function', warningList )
                    form = photonScatteringModule.incoherent.form( info.style, standardsModule.frames.labToken, subform )
                product.distribution.add( form )
        elif( MT in [ 526, 527, 528 ] ) :
            if( isTwoBody ) :
                outputChannel = channelsModule.twoBodyOutputChannel( process = process )
            else :
                outputChannel = channelsModule.NBodyOutputChannel( process = process )
            outputChannel.Q.add( toGNDMiscModule.returnConstantQ( info.style, 0 ) )
            productsNeeded = [ 'e-' ]
            if( MT == 527 ) : productsNeeded.insert( 0, 'gamma' )
        else :
            outputChannel = channelsModule.NBodyOutputChannel( process = process )
            outputChannel.Q.add( toGNDMiscModule.returnConstantQ( info.style, EPE ) )
            if( MT in MT_AtomicConfigurations ) :
                productsNeeded = [ 'e-' ]
            elif( MT in [ 515, 517 ] ) :
                productsNeeded = [ 'e-', 'e+' ]
            else :
                addTargetAsResidual = False

        if( MT in [ 526, 527, 528 ] ) : productList[0].addAttribute( 'ENDFconversionFlag', 'MF26' )
        for name in productsNeeded :
            product = None
            for i1, product in enumerate( productList ) :
                if( name == product.name ) : break
                product = None
            if( product is None ) :
                product = toGNDMiscModule.getTypeNameGamma( info, { "gamma" : 0, "e+" : 8, "e-" : 9 }[name] )
                product = toGNDMiscModule.newGNDParticle( info, product )
            outputChannel.products.add( outputChannel.products.uniqueLabel( product ) )
        if( ( MT >= 534 ) and ( reactionSuite.projectile.name == 'e-' ) ) :
            product = toGNDMiscModule.getTypeNameGamma( info, 9 )
            product = toGNDMiscModule.newGNDParticle( info, product )
            outputChannel.products.add( outputChannel.products.uniqueLabel( product ) )

        if( addTargetAsResidual ) :
            elementSymbol = nuclearModule.elementSymbolFromZ( info.targetZA / 1000 )
            if( MT >= 534 ) : elementSymbol += '{%s}' % MT_AtomicConfigurations[MT]
            element = xParticleModule.element( elementSymbol )
            residual = toGNDMiscModule.newGNDParticle( info, element )
            outputChannel.products.add( outputChannel.products.uniqueLabel( residual ) )

        if( len( MTData ) > 0 ) : raise Exception( 'Untranslated MF data: MFs = %s' % MTData.keys( ) )

        if( MT in summedReactions ) :
            summedReactions[MT] = [ crossSection, outputChannel ]
            info.logs.write( '\n' )
            continue

        EFL = PQUModule.PQU( PQUModule.pqu_float.surmiseSignificantDigits( EFL ), 'eV' )
        reaction = reactionModule.reaction( outputChannel, "%s" % iChannel, ENDF_MT = MT, date = info.Date, EFL = EFL )
        reaction.crossSection.add( crossSection )

        reactionSuite.reactions.add( reaction )
        iChannel += 1

        info.logs.write( '\n' )
        for warning in warningList : info.logs.write( "       WARNING: %s\n" % warning, stderrWriting = True )

    if( True ) :            # Need to test for gamma as projectile.
        for MT in sorted( summedReactions.keys( ) ) :
            if( summedReactions[MT] is None ) : continue
            crossSection, outputChannel = summedReactions[MT]
            summands = [ ]
            summandMTs = { 501 : range( 502, 573 ), 516 : ( 515, 517 ), 522 : range( 534, 573 ) }[MT]

            for reaction in reactionSuite :
                if( reaction.ENDF_MT in summandMTs ) :
                    summands.append( sumsModule.add( link = reaction.crossSection ) )
            summedCrossSection = sumsModule.crossSectionSum( name=photonSumLabels[MT], label=str(iChannel), ENDF_MT=MT,
                    summands = sumsModule.listOfSummands( summandList=summands ), crossSection=crossSection, date = info.Date )
            Q = outputChannel.Q[info.style]
            if( Q is not None ) : summedCrossSection.Q.add( Q )
            reactionSuite.sums.add( summedCrossSection )
            iChannel += 1
