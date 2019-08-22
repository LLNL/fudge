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

import endfFormats as endfFormatsModule
import gndToENDF6 as gndToENDF6Module
import fudge.particles.nuclear as nuclear
from pqu import PQU
from fudge.legacy.converting import endf_endl
from fudge.structure import masses

import fudge

import resonances as resonancesModule
import productData.multiplicity as multiplicityModule

from fudge.gnd import styles as stylesModule
import fudge.gnd.reactionSuite as reactionSuiteModule
import fudge.gnd.covariances.covarianceSuite as covarianceSuiteModule

import fudge.processing.processingInfo as processingInfoModule

import xData.XYs as XYsModule

elementalMass = { 1000 : 9.992414e-01,  2000 : 3.968215e+00,  3000 : 6.881371e+00,  4000 : 8.934758e+00,  5000 : 1.071713e+01,
                  6000 : 1.190782e+01,  7000 : 1.388637e+01,  8000 : 1.586195e+01,  9000 : 1.883519e+01, 10000 : 2.000565e+01,
                 11000 : 2.279230e+01, 12000 : 2.409620e+01, 13000 : 2.674971e+01, 14000 : 2.784422e+01, 15000 : 3.070771e+01,
                 16000 : 3.178458e+01, 17000 : 3.514843e+01, 18000 : 3.960482e+01, 19000 : 3.876242e+01, 20000 : 3.973568e+01,
                 21000 : 4.456969e+01, 22000 : 4.748850e+01, 23000 : 5.050387e+01, 24000 : 5.154932e+01, 25000 : 5.446604e+01,
                 26000 : 5.536723e+01, 27000 : 5.842692e+01, 28000 : 5.819572e+01, 29000 : 6.300009e+01, 30000 : 6.481834e+01,
                 31000 : 6.912106e+01, 32000 : 7.196640e+01, 33000 : 7.427797e+01, 34000 : 7.828168e+01, 35000 : 7.921757e+01,
                 36000 : 8.308009e+01, 37000 : 8.473357e+01, 38000 : 8.686728e+01, 39000 : 8.814213e+01, 40000 : 9.043635e+01,
                 41000 : 9.210826e+01, 42000 : 9.511580e+01, 43000 : 9.715810e+01, 44000 : 1.002017e+02, 45000 : 1.020210e+02,
                 46000 : 1.054859e+02, 47000 : 1.069413e+02, 48000 : 1.114443e+02, 49000 : 1.138336e+02, 50000 : 1.176704e+02,
                 51000 : 1.207041e+02, 52000 : 1.265038e+02, 53000 : 1.258138e+02, 54000 : 1.301720e+02, 55000 : 1.317632e+02,
                 56000 : 1.361502e+02, 57000 : 1.377117e+02, 58000 : 1.389163e+02, 59000 : 1.396975e+02, 60000 : 1.430009e+02,
                 61000 : 1.437543e+02, 62000 : 1.491080e+02, 63000 : 1.506545e+02, 64000 : 1.558991e+02, 65000 : 1.575597e+02,
                 66000 : 1.611040e+02, 67000 : 1.635131e+02, 68000 : 1.658231e+02, 69000 : 1.674827e+02, 70000 : 1.715535e+02,
                 71000 : 1.734639e+02, 72000 : 1.769566e+02, 73000 : 1.793935e+02, 74000 : 1.822706e+02, 75000 : 1.846073e+02,
                 76000 : 1.885660e+02, 77000 : 1.905687e+02, 78000 : 1.934140e+02, 79000 : 1.952739e+02, 80000 : 1.988668e+02,
                 81000 : 2.026143e+02, 82000 : 2.054200e+02, 83000 : 2.071847e+02, 84000 : 2.072045e+02, 85000 : 2.081959e+02,
                 86000 : 2.200928e+02, 87000 : 2.210843e+02, 88000 : 2.240585e+02, 89000 : 2.250499e+02, 90000 : 2.300446e+02,
                 91000 : 2.290155e+02, 92000 : 2.359841e+02, 93000 : 2.349640e+02, 94000 : 2.419039e+02, 95000 : 2.409124e+02,
                 96000 : 2.448781e+02, 97000 : 2.448781e+02, 98000 : 2.488437e+02, 99000 : 2.498351e+02, 100000 : 2.547922e+02 }

__metaclass__ = type

def toENDF6( self, style, flags, verbosityIndent = '', covarianceSuite = None ) :

    evaluatedStyle = self.styles.getEvaluatedStyle( )
    if( evaluatedStyle is None ) : raise ValueError( 'no evaluation style found' )

    if( flags['verbosity'] >= 10 ) : print '%s%s' % ( verbosityIndent, self.inputParticlesToReactionString( suffix = " -->" ) )
    verbosityIndent2 = verbosityIndent + ' ' * ( len( self.inputParticlesToReactionString( suffix = " -->" ) ) + 1 )
    projectile, target = self.projectile, self.target
    projectileZA = projectile.getZ_A_SuffixAndZA( )[-1]
    IPART = projectileZA
    if( projectile.name == 'e-' ) : IPART = 11
    targetZA, MAT = endf_endl.ZAAndMATFromParticleName( target.name )
    targetZ, targetA = divmod( targetZA, 1000 )
    targetInfo = processingInfoModule.tempInfo( )
    targetInfo['style'] = style
    targetInfo['reactionSuite'] = self
    targetInfo['ZA'] = targetZA
    if( self.particles.hasID( 'n' ) ) :       # Need neutron mass in eV/c**2, but it may not be in the particle list.
        targetInfo['neutronMass'] = self.getParticle( 'n' ).getMass( 'eV/c**2' )
    else :
        neutronAmu = masses.getMassFromZA( 1 )
        targetInfo['neutronMass'] = PQU.PQU( neutronAmu, 'amu' ).getValueAs('eV/c**2')
    if( isinstance( target, fudge.gnd.xParticle.element ) ) :
        targetInfo['mass'] = elementalMass[targetZA]
    else :
        targetInfo['mass'] = target.getMass( 'eV/c**2' ) / targetInfo['neutronMass']

    try :
        targetInfo['LIS'] = target['levelIndex']
    except :
        targetInfo['LIS'] = 0
    targetInfo['metastables'] = []
    targetInfo['LISO'] = 0
    for key, alias in self.aliases.items( ) :
        if( alias.hasAttribute( 'nuclearMetaStable' ) ) :
            targetInfo['metastables'].append( alias.getValue() )
            if( alias.getValue() == target.name ) :
                targetInfo['LISO'] = int( alias.getAttribute( 'nuclearMetaStable' ) )
    MAT += targetInfo['LISO']
    if( self.MAT is not None ) : MAT = self.MAT

    ITYPE = 0                   # Other ITYPE sublibraries not yet supported. BRB is this still true
    for reaction in self.reactions :
        if( 500 <= reaction.ENDF_MT < 573 ) : ITYPE = 3
    targetInfo['crossSectionMF'] = { 0 : 3, 3 : 23 }[ITYPE]

    targetInfo['delayedRates'] = []
    targetInfo['totalDelayedNubar'] = None
    targetInfo['MTs'], targetInfo['MF8'], targetInfo['LRs'] = {}, {}, {}
    endfMFList = { 1 : { 451 : [] }, 2 : {}, 3 : {}, 4 : {}, 5 : {}, 6 : {}, 8 : {}, 9 : {}, 10 : {}, 12 : {}, 13 : {},
            14 : {}, 15 : {}, 23 : {}, 26 : {}, 27 : {}, 31 : {}, 32 : {}, 33 : {}, 34 : {}, 35 : {}, 40 : {} }
    if( self.resonances is not None ) :      # Add resonances, independent of reaction channels
        self.resonances.toENDF6( endfMFList, flags, targetInfo, verbosityIndent=verbosityIndent2 )

    targetInfo['production_gammas'] = {}

    for reaction in self :
        reaction.toENDF6( endfMFList, flags, targetInfo, verbosityIndent = verbosityIndent2 )
    gndToENDF6Module.upDateENDFMF8Data( endfMFList, targetInfo )
    for MT, production_gammas in targetInfo['production_gammas'].items( ) :
        MF, production_gammas = production_gammas[0], production_gammas[1:]
        for productionReaction in production_gammas :
            gammas = [ gamma for gamma in productionReaction.outputChannel ]
            targetInfo['crossSection'] = productionReaction.crossSection[targetInfo['style']]
            gndToENDF6Module.gammasToENDF6_MF12_13( MT, MF, endfMFList, flags, targetInfo, gammas )

    for particle in self.particles :              # gamma decay data.
        if( isinstance( particle, fudge.gnd.xParticle.isotope ) ) :
            for level in particle :
                if( level.gammas ) :                        # non-empty gamma information
                    for baseMT in [ 50, 600, 650, 700, 750, 800 ] :
                        residualZA = endf_endl.ENDF_MTZAEquation( projectileZA, targetZA, baseMT )[0][-1]
                        if( nuclear.nucleusNameFromZA( residualZA ) == particle.name ) : break
                    level.toENDF6( baseMT, endfMFList, flags, targetInfo )

    MFs = sorted( endfMFList.keys( ) )
    endfList = []

    totalNubar = None
    totalDelayedNubar = targetInfo['totalDelayedNubar']
    if( 455 in endfMFList[5] ) :
        MF5MT455s = endfMFList[5][455]

        endfMFList[1][455]  = [ endfFormatsModule.endfHeadLine( targetZA, targetInfo['mass'], 0, 2, 0, 0 ) ] # Currently, only LDG = 0, LNU = 2 is supported.
        endfMFList[1][455] += [ endfFormatsModule.endfHeadLine( 0, 0, 0, 0, len( targetInfo['delayedRates'] ), 0 ) ]
        endfMFList[1][455] += endfFormatsModule.endfDataList( targetInfo['delayedRates'] )

        multiplicityModule.fissionNeutronsToENDF6( 455, totalDelayedNubar, endfMFList, flags, targetInfo )

        MF5MT455List = [ endfFormatsModule.endfHeadLine( targetZA, targetInfo['mass'], 0, 0, len( MF5MT455s ), 0 ) ]
        for MF5MT455 in MF5MT455s : MF5MT455List += MF5MT455
        if( len( MF5MT455s ) == 0 ) :
            del endfMFList[5][455]
        else :
            endfMFList[5][455] = MF5MT455List + [ endfFormatsModule.endfSENDLineNumber( ) ]
    if(   'promptNubar' in targetInfo.dict ) :
        promptNubar = targetInfo['promptNubar']
        multiplicityModule.fissionNeutronsToENDF6( 456, promptNubar, endfMFList, flags, targetInfo )
        totalNubar = promptNubar
        try :
            if( not( totalDelayedNubar is None ) ) : totalNubar = totalNubar + totalDelayedNubar
        except :                                # The following is a kludge for some "bad" data.
            if( ( totalNubar.domainMax( unitTo = 'MeV' ) == 30. ) and
                ( totalDelayedNubar.domainMax( unitTo = 'MeV' ) == 20. ) ) :
                    totalDelayedNubar[-1] = [ totalNubar.domainMax( ), totalDelayedNubar.getValue( totalDelayedNubar.domainMax( ) ) ]
            totalNubar = totalNubar + totalDelayedNubar
    elif( 'totalNubar' in targetInfo.dict ) :
        totalNubar = targetInfo['totalNubar']
    if( totalNubar is not None ) :
        multiplicityModule.fissionNeutronsToENDF6( 452, totalNubar, endfMFList, flags, targetInfo )

    if( covarianceSuite ) : covarianceSuite.toENDF6( endfMFList, flags, targetInfo )

    endfDoc = self.documentation.get( 'endfDoc' )
    if( endfDoc is None ) :
        docHeader2 = [  ' %2d-%-2s-%3d LLNL       EVAL-OCT03 Unknown' % ( targetZ, fudge.particles.nuclear.elementSymbolFromZ( targetZ ), targetA ),
                        '                      DIST-DEC99                       19990101   ',
                        '----ENDL              MATERIAL %4d' % MAT,
                        '-----INCIDENT %s DATA' %
                            { 1 : 'NEUTRON', 1001 : 'PROTON', 1002 : 'DEUTERON', 1003 : 'TRITON', 2003 : 'HELION', 2004 : 'ALPHA' }[projectileZA],
                        '------ENDF-6 FORMAT' ]
        endfDoc = [ 'LLNL ENDL file translated to ENDF6 by FUDGE.', '' ' ************************ C O N T E N T S ***********************' ]
    else :
        docHeader2 = []
        endfDoc = endfDoc.getLines( )

        # update the documentation, including metadata on first 4 lines:
    try :
        self.getReaction( 'fission' )
        LFI = True
    except KeyError :
        LFI = False
    LRP = -1
    if( self.resonances is not None ) :
        if( self.resonances.scatteringRadius ) :
            LRP = 0
        elif( self.resonances.reconstructCrossSection ) :
            LRP = 1
        elif( self.resonances.unresolved and not( self.resonances.resolved )
                and self.resonances.unresolved.tabulatedWidths.forSelfShieldingOnly ) :
            LRP = 1
        else :
            LRP = 2
    EMAX = max( [ reaction.crossSection.domainMax( unitTo = 'eV' ) for reaction in self.reactions ] )

    temperature = self.styles[style].temperature.getValueAs( 'K' )
    library = evaluatedStyle.library
    version = evaluatedStyle.version
    if( library == 'ENDL' ) :           # Additional ENDF meta-data. If the library is unknown, use NLIB = -1
        NVER, LREL, NMOD = 1, 1, 1
        NLIB = -1
    else :
        NVER, LREL, NMOD = map( int, version.split( '.' ) )    # Version stored as '7.2.1'
        NLIB = { "ENDF/B" :  0,     "ENDF/A" :  1,      "JEFF"                 :  2,    "EFF"      :  3,    "ENDF/B (HE)" :  4,
                 "CENDL"  :  5,     "JENDL"  :  6,      "SG-23"                : 21,    "INDL/V"   : 31,    "INDL/A"      : 32,
                 "FENDL"  : 33,     "IRDF"   : 34,      "BROND (IAEA version)" : 35,    "INGDB-90" : 36,    "FENDL/A"     : 37,
                 "BROND"  : 41 }.get( library, -1 )

    NFOR = 6    # ENDF-6 format
    NSUB = 10 * IPART + ITYPE
    LDRV = 0
    STA = 0
    if( isinstance( self.target, fudge.gnd.xParticle.nuclearLevel ) or self.target.attributes.get( 'unstable' ) ) : STA = 1
    if( targetInfo['LISO'] ) : STA = 1
    levelIndex, level_eV = 0, 0.
    if( hasattr( self.target, 'getLevelIndex' ) ) : levelIndex, level_eV = self.target.getLevelIndex( ), self.target.getLevelAsFloat( 'eV' )
    docHeader = [ endfFormatsModule.endfHeadLine( targetZA, targetInfo['mass'], LRP, LFI, NLIB, NMOD ),
            endfFormatsModule.endfHeadLine( level_eV, STA, levelIndex, targetInfo['LISO'], 0, NFOR ),
            endfFormatsModule.endfHeadLine( self.projectile.getMass( 'eV/c**2' ) / targetInfo['neutronMass'], EMAX, LREL, 0, NSUB, NVER ),
            endfFormatsModule.endfHeadLine( temperature, 0, LDRV, 0, len( endfDoc ), -1 ) ]
    new_doc = fudge.gnd.documentation.documentation( 'endf', '\n'.join( docHeader + docHeader2 + endfDoc ) )
    endfMFList[1][451] += endfFormatsModule.toEndfStringList( new_doc )

    return( endfFormatsModule.endfMFListToFinalFile( endfMFList, MAT, lineNumbers = True ) )

reactionSuiteModule.reactionSuite.toENDF6 = toENDF6
