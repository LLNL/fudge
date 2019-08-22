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

import site_packages.legacy.toENDF6.gndToENDF6 as gndToENDF6Module
import site_packages.legacy.toENDF6.endfFormats as endfFormatsModule
import xData.XYs as XYsModule
import xData.regions as regionsModule
import fudge.gnd.tokens as tokensModule
import fudge.gnd.channels as channelsModule
import fudge.gnd.reactions.reaction as reactionModule
import fudge.gnd.productData.multiplicity as multiplicityModule
import fudge.legacy.converting.endf_endl as endf_endlModule

def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ) :

    def addDecayProducts( parent, products ) :

        if( parent.decayChannel is None ) : return( False )
        doMF4AsMF6 = False
        for product in parent.decayChannel :
            if( len( product.distribution ) and ( product.name != 'gamma' ) ) :
                if( product.attributes.get( 'ENDFconversionFlag' ) == 'MF6' ) : doMF4AsMF6 = True
                products.append( product )
        return( doMF4AsMF6 )

    def divideIgnoring0DividedBy0( self, other ) :

        if isinstance( self, XYsModule.XYs ) and isinstance( other, XYsModule.XYs ):
            d = self.union( other.domainSlice( domainMin = self.domainMin( ), domainMax = self.domainMax( ) ) )
            result = []
            for p in d :
                vp = other.evaluate( p[0] )
                if( vp == 0 ) :
                    if( p[1] != 0 ) : raise Exception( 'Divide non-zero number by zero at %e' % p[0] )
                else :
                    p[1] = p[1] / vp
                result.append( [ p[0], p[1] ] )
            return( XYsModule.pointwiseXY( data = result ) )

        elif ( isinstance(self, regionsModule.regions) and isinstance(other, regionsModule.regions)
                and len(self)==len(other) ):
            regions = regionsModule.regions( dimension=1 )
            for idx in range(len(self)):
                regions.append( XYsModule.XYs( data = divideIgnoring0DividedBy0( self[idx], other[idx] ) ) )
            return regions

        else:
            raise NotImplementedError( 'Divide XYs by regions or vice-versa' )

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

    outputChannel = self.outputChannel
    if isinstance( outputChannel, channelsModule.productionChannel ) :
        print '        toENDF6 does not support writing of "%s" channel' % outputChannel.genre
        return
    MT = self.ENDF_MT
    if( flags['verbosity'] >= 10 ) : print '%s%s' % ( verbosityIndent, outputChannel.toString( simpleString = True ) )
    targetInfo['Q'] = self.getQ( 'eV', groundStateQ = True )
    targetInfo['QM'] = None

    LR, tryLR = 0, False
    if isinstance( outputChannel, channelsModule.twoBodyOutputChannel ) :
        tryLR = True
    elif isinstance( outputChannel, channelsModule.NBodyOutputChannel ) :
        if( MT == 91 ) : tryLR = True
    if( tryLR and ( len( outputChannel ) > 1 ) ) :
        secondProduct = outputChannel[1]
        primaryResidualName, decayProducts, decayChannel, = secondProduct.name.split( '_' )[0], [], secondProduct.decayChannel
        if( not( decayChannel is None ) ) :
            for decayProduct in decayChannel :
                decayProductName = decayProduct.name
                if( decayProductName not in [ primaryResidualName, 'gamma' ] ) : decayProducts.append( decayProductName )
        if( len( decayProducts ) == 1 ) :   # Kludge for Carbon breakup into 3 alphas.
            if( ( primaryResidualName == 'C' ) and ( decayProducts == [ 'He4' ] ) ) : LR = 23
        elif( len( decayProducts ) > 1 ) :                                        # This must be a breakup reaction.
            MTProducts = endf_endlModule.endfMTtoC_ProductList( 0, '' )
            MTProducts.productCounts[outputChannel[0].name] += 1
            for decayProduct in decayProducts[:-1] : MTProducts.productCounts[decayProduct] += 1
            for MT_LR in [ 22, 23, 24, 25, 28, 29, 30, 32, 33, 34, 35, 36 ] :   # 39 and 40 not allowed in ENDF6
                if( endf_endlModule.endfMTtoC_ProductLists[MT_LR].productCounts == MTProducts.productCounts ) :
                    LR = MT_LR
                    break
            if( ( LR == 32 ) and ( primaryResidualName == 'B10' ) and ( decayProducts[-1] == 'He4' ) ) : LR = 35   # Kludge for bad data.
        if( LR != 0 ) :
            QM = outputChannel.Q[targetInfo['style']].getValue( 0, 'eV' ) + decayChannel.Q[targetInfo['style']].getValue( 0, 'eV' )
            targetInfo['QM'] = QM
    targetInfo['LRs'][MT] = LR

    level = 0.
    for product in outputChannel :
        tmp = None
        if( hasattr( product, 'getLevelAsFloat' ) ) : tmp = product.getLevelAsFloat( 'eV' )
        if( tmp is not None ) : level = tmp

    targetInfo['EMin'], targetInfo['EMax'] = self.crossSection.domain( unitTo="eV" )
    if( self.EFL is not None ) : targetInfo['EFL'] = self.EFL.getValueAs( "eV" )
    self.crossSection.toENDF6( MT, endfMFList, targetInfo, level, LR )
    if( 'EFL' in targetInfo ) : del targetInfo['EFL']

    products = []
    doMF4AsMF6 = False
    for product in outputChannel :
        if( len( product.distribution ) ) : products.append( product )
        doMF4AsMF6 = doMF4AsMF6 or addDecayProducts( product, products )
    if( outputChannel.isFission( ) ) : 
        doMF4AsMF6 = False
        if( not( outputChannel.fissionEnergyReleased is None ) ) :
            outputChannel.fissionEnergyReleased.toENDF6( 1, endfMFList, flags, targetInfo )
        delayedNubar = None                                 # Special work to get delayed nubar weights.
        for product in outputChannel :
            if( 'emissionMode' in product.attributes ) :
                if( product.getAttribute( 'emissionMode' ) == tokensModule.delayedToken ) :
                    if( delayedNubar is None ) :
                        delayedNubar = product.multiplicity[targetInfo['style']]
                    else :
                        delayedNubar = delayedNubar + product.multiplicity[targetInfo['style']]
                    targetInfo['delayedRates'].append( product.getAttribute( 'decayRate' ).getValueAs('1/s') )
        if( delayedNubar is not None ) :
            targetInfo['totalDelayedNubar'] = delayedNubar
            for product in outputChannel :
                if( 'emissionMode' in product.attributes ) :
                    if( product.getAttribute( 'emissionMode' ) == tokensModule.delayedToken ) :
                        weight = divideIgnoring0DividedBy0( product.multiplicity[targetInfo['style']], delayedNubar )
                        product.ENDF6_delayedNubarWeights = thinWeights( weight )
    if( MT in [ 526, 527 ] ) : doMF4AsMF6 = True
    targetInfo['doMF4AsMF6'] = doMF4AsMF6
    targetInfo['MF6LCTs'], targetInfo['gammas'] = [], []
    for productIndex, product in enumerate( outputChannel ) :
        if( ( product.name == 'gamma' ) and not isinstance( outputChannel, channelsModule.twoBodyOutputChannel ) and ( MT != 527 ) ) :
            targetInfo['gammas'].append( product )
            continue
        targetInfo['productIndex'] = str( productIndex )
        targetInfo['productToken'] = product.name
        targetInfo['productLabel'] = product.label
        product.toENDF6( MT, endfMFList, flags, targetInfo, verbosityIndent = verbosityIndent + '    ' )
    gammas = targetInfo['gammas']
    if( len( gammas ) ) :
        gamma = gammas[0]
        targetInfo['zapID'] = gamma.name
        targetInfo['particleMass'] = gamma.getMass( 'eV/c**2' )
        multiplicity = gamma.multiplicity
        if( gamma.getAttribute( 'ENDFconversionFlag' ) == 'MF6' ) :
            targetInfo['multiplicity'] = gamma.multiplicity
            gndToENDF6Module.gammasToENDF6_MF6( MT, endfMFList, flags, targetInfo, gammas )
        elif( gamma.getAttribute( 'ENDFconversionFlag' ) == 'MF13' ) :
            targetInfo['crossSection'] = self.crossSection[targetInfo['style']]
            gndToENDF6Module.gammasToENDF6_MF12_13( MT, 13, endfMFList, flags, targetInfo, gammas )
        elif( isinstance( multiplicity, multiplicityModule.constant ) ) :      # This should probably be changed to unknown in some cases?????
            pass
        elif( isinstance( multiplicity, multiplicityModule.unknown ) ) :      # This should probably be changed to unknown in some cases?????
            pass
        else :
            gndToENDF6Module.gammasToENDF6_MF12_13( MT, 12, endfMFList, flags, targetInfo, gammas )
    if( len( targetInfo['MF6LCTs'] ) > 0 ) :
        LCT = targetInfo['MF6LCTs'][0]
        for i in targetInfo['MF6LCTs'] :
            if( not( LCT is None ) ) : break
        if( LCT is None ) : LCT = 2
        for i in targetInfo['MF6LCTs'] :
            if( i is None ) : continue
            if( i != LCT ) : LCT = 3
        if( 500 <= MT <= 572 ) : LCT = 0
        MF6or26 = { 3 : 6, 23 : 26 }[targetInfo['crossSectionMF']]
        endfMFList[MF6or26][MT].insert( 0, endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], 0, LCT, len( targetInfo['MF6LCTs'] ), 0 ) )
        endfMFList[MF6or26][MT].append( endfFormatsModule.endfSENDLineNumber( ) )

reactionModule.reaction.toENDF6 = toENDF6
