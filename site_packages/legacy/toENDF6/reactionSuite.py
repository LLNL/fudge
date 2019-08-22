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

import textwrap

from PoPs import IDs as IDsPoPsModule
from PoPs.quantities import halflife as halflifePoPsModule
from PoPs.families import gaugeBoson as gaugeBosonPoPsModule
from PoPs.families import lepton as leptonPoPsModule
from PoPs.families import baryon as baryonPoPsModule
from PoPs.families import nuclide as nuclidePoPsModule
from PoPs.families import nucleus as nucleusPoPsModule
from PoPs.groups import isotope as isotopePoPsModule
from PoPs.groups import chemicalElement as chemicalElementPoPsModule
from PoPs.groups import misc as chemicalElementMiscPoPsModule

from fudge.legacy.converting import endf_endl as endf_endlModule
from fudge.legacy.converting import massTracker as massTrackerModule

from fudge.gnds import documentation as documentationModule
from fudge.gnds import reactionSuite as reactionSuiteModule

from . import endfFormats as endfFormatsModule
from . import gndsToENDF6 as gndsToENDF6Module
from . import ENDFconversionFlags as ENDFconversionFlagsModule
from .productData import multiplicity as multiplicityModule

__metaclass__ = type

def toENDF6( self, style, flags, verbosityIndent = '', covarianceSuite = None, useRedsFloatFormat = False,
             lineNumbers = True, **kwargs ) :

    _useRedsFloatFormat = endfFormatsModule.useRedsFloatFormat
    endfFormatsModule.useRedsFloatFormat = useRedsFloatFormat

    evaluatedStyle = self.styles.getEvaluatedStyle( )
    if( evaluatedStyle is None ) : raise ValueError( 'no evaluation style found' )

    if( flags == {} ) : flags['verbosity'] = 0
    if( flags['verbosity'] >= 10 ) :
        print ( '%s%s' % ( verbosityIndent, self.inputParticlesToReactionString( suffix = " -->" ) ) )
    verbosityIndent2 = verbosityIndent + ' ' * ( len( self.inputParticlesToReactionString( suffix = " -->" ) ) + 1 )

    projectileZA = chemicalElementMiscPoPsModule.ZA( self.PoPs[self.projectile] )
    IPART = projectileZA
    if( self.projectile == 'e-' ) : IPART = 11

    targetID = self.target
    if( targetID in self.PoPs.aliases ) : targetID = self.PoPs[targetID].pid
    targetZA, MAT = endf_endlModule.ZAAndMATFromParticleName( targetID )
    targetZ, targetA = divmod( targetZA, 1000 )

    targetInfo = {}
    targetInfo['massTracker'] = massTrackerModule.massTracker()
    for particle in self.PoPs :
        if( isinstance( particle, ( gaugeBosonPoPsModule.particle, leptonPoPsModule.particle, baryonPoPsModule.particle, nuclidePoPsModule.particle ) ) ) : 
            ZA = chemicalElementMiscPoPsModule.ZA( particle )
            if( isinstance( particle, nuclidePoPsModule.particle ) ) :
                if( particle.index != 0 ) : continue
            if( len( particle.mass ) > 0 ) : targetInfo['massTracker'].addMassAMU( ZA, particle.getMass( 'amu' ) )

    targetInfo['ENDFconversionFlags'] = {}
    if 'LLNL' in self.applicationData:
        conversionFlags = [data for data in self.applicationData['LLNL']
                           if isinstance(data, ENDFconversionFlagsModule.ENDFconversionFlags)]
        if len(conversionFlags) > 0:
            for link in conversionFlags[0].flags:
                targetInfo['ENDFconversionFlags'][link.link] = link['flags']

    for chemicalElement in self.PoPs.chemicalElements :
        ZA = 1000 * chemicalElement.Z
        try :
            targetInfo['massTracker'].getMassAMU( ZA )
        except :                                                # If not present, i.e., a raise, add.
            targetInfo['massTracker'].addMassAMU( ZA, targetInfo['massTracker'].getElementalMassAMU( ZA ) )

    targetInfo['style'] = style
    targetInfo['reactionSuite'] = self
    targetInfo['ZA'] = targetZA

    try :
        target = self.PoPs[targetID]
    except :
        target = self.PoPs.chemicalElements.getSymbol( targetID )

    levelIndex = 0
    levelEnergy_eV = 0
    if( isinstance( target, nuclidePoPsModule.particle ) ) :      # isomer target
        levelIndex = target.index
        levelEnergy_eV = target.energy[0].float( 'eV' )
    targetInfo['mass'] = targetInfo['massTracker'].getMassAWR( targetZA, levelEnergyInEv = levelEnergy_eV )

    targetInfo['LIS'] = levelIndex
    targetInfo['metastables'] = {}
    targetInfo['LISO'] = 0
    if( levelIndex > 0 ) : targetInfo['LISO'] = 1
# BRBBRB
    for alias in self.PoPs.aliases :
        if( hasattr( alias, 'metaStableIndex' ) ) :
            targetInfo['metastables'][alias.pid] = alias
    MAT += targetInfo['LISO']
    if( self.MAT is not None ) : MAT = self.MAT

    ITYPE = 0
    for reaction in self.reactions :
        if( 500 <= reaction.ENDF_MT < 573 ) : ITYPE = 3
    targetInfo['crossSectionMF'] = { 0 : 3, 3 : 23 }[ITYPE]

    if ITYPE == 3:
        targetInfo['EFL'] = 0

    targetInfo['delayedRates'] = []
    targetInfo['MTs'], targetInfo['MF8'], targetInfo['LRs'] = {}, {}, {}
    endfMFList = { 1 : { 451 : [] }, 2 : {}, 3 : {}, 4 : {}, 5 : {}, 6 : {}, 8 : {}, 9 : {}, 10 : {}, 12 : {}, 13 : {},
            14 : {}, 15 : {}, 23 : {}, 26 : {}, 27 : {}, 31 : {}, 32 : {}, 33 : {}, 34 : {}, 35 : {}, 40 : {} }
    if( self.resonances is not None ) :      # Add resonances, independent of reaction channels
        self.resonances.toENDF6( endfMFList, flags, targetInfo, verbosityIndent=verbosityIndent2 )

    targetInfo['production_gammas'] = {}

    for multiplicitySum in self.sums.multiplicities:
        if multiplicitySum.ENDF_MT == 452:
            targetInfo['totalNubar'] = multiplicitySum.multiplicity.evaluated
        elif multiplicitySum.ENDF_MT == 455:
            targetInfo['totalDelayedNubar'] = multiplicitySum.multiplicity.evaluated

    for reaction in self :
        reaction.toENDF6( endfMFList, flags, targetInfo, verbosityIndent = verbosityIndent2 )
    gndsToENDF6Module.upDateENDFMF8Data( endfMFList, targetInfo )
    for MT, production_gammas in targetInfo['production_gammas'].items( ) :
        MF, production_gammas = production_gammas[0], production_gammas[1:]
        for productionReaction in production_gammas :
            gammas = [ gamma for gamma in productionReaction.outputChannel ]
            targetInfo['crossSection'] = productionReaction.crossSection[targetInfo['style']]
            gndsToENDF6Module.gammasToENDF6_MF12_13( MT, MF, endfMFList, flags, targetInfo, gammas )

    for particle in self.PoPs :
        if( isinstance( particle, ( nucleusPoPsModule.particle, nuclidePoPsModule.particle ) ) ) :
            if( len( particle.decayData.decayModes ) > 0 ) :
                for baseMT in [ 50, 600, 650, 700, 750, 800, 1000 ] :   # 1000 causes raise in endf_endlModule.ENDF_MTZAEquation.
                    residualZA = endf_endlModule.ENDF_MTZAEquation( projectileZA, targetZA, baseMT )[0][-1]
                    if( chemicalElementMiscPoPsModule.ZA( particle ) == residualZA ) : break
                addDecayGamma( self.PoPs, particle, baseMT, endfMFList, flags, targetInfo )

    if 'totalNubar' in targetInfo:
        multiplicityModule.fissionNeutronsToENDF6( 452, targetInfo['totalNubar'], endfMFList, flags, targetInfo )
    if 'promptNubar' in targetInfo:
        multiplicityModule.fissionNeutronsToENDF6( 456, targetInfo['promptNubar'], endfMFList, flags, targetInfo )
    if 'totalDelayedNubar' in targetInfo:
        MF5MT455s = endfMFList[5][455]

        endfMFList[1][455]  = [ endfFormatsModule.endfHeadLine( targetZA, targetInfo['mass'], 0, 2, 0, 0 ) ] # Currently, only LDG = 0, LNU = 2 is supported.
        endfMFList[1][455] += [ endfFormatsModule.endfHeadLine( 0, 0, 0, 0, len( targetInfo['delayedRates'] ), 0 ) ]
        endfMFList[1][455] += endfFormatsModule.endfDataList( targetInfo['delayedRates'] )

        multiplicityModule.fissionNeutronsToENDF6( 455, targetInfo['totalDelayedNubar'], endfMFList, flags, targetInfo )

        MF5MT455List = [ endfFormatsModule.endfHeadLine( targetZA, targetInfo['mass'], 0, 0, len( MF5MT455s ), 0 ) ]
        for MF5MT455 in MF5MT455s : MF5MT455List += MF5MT455
        if( len( MF5MT455s ) == 0 ) :
            del endfMFList[5][455]
        else :
            endfMFList[5][455] = MF5MT455List + [ endfFormatsModule.endfSENDLineNumber( ) ]

    if( covarianceSuite ) : covarianceSuite.toENDF6( endfMFList, flags, targetInfo )

    endfDoc = self.documentation.get( 'endfDoc' )
    if( endfDoc is None ) :
        docHeader2 = [  ' %2d-%-2s-%3d LLNL       EVAL-OCT03 Unknown' % ( targetZ, chemicalElementMiscPoPsModule.symbolFromZ[targetZ], targetA ),
                        '                      DIST-DEC99                       19990101   ',
                        '----ENDL              MATERIAL %4d' % MAT,
                        '-----INCIDENT %s DATA' %
                            { 1 : 'NEUTRON', 1001 : 'PROTON', 1002 : 'DEUTERON', 1003 : 'TRITON', 2003 : 'HELION', 2004 : 'ALPHA' }[projectileZA],
                        '------ENDF-6 FORMAT' ]
        endfDoc = [ 'LLNL ENDL file translated to ENDF6 by FUDGE.', '' ' ************************ C O N T E N T S ***********************' ]
        endlDoc = self.documentation.get( 'ENDL' ).getLines() 
        endlDoc2 = [ ]
        if endlDoc != None: 
            for line in endlDoc:
                if 'endep' in line:    # remove text from the line before the ones containing 'endep'
                    del endlDoc2[-1]
                    break
                newline = textwrap.wrap(line,width=66,drop_whitespace=False,subsequent_indent='    ')
                if len(newline)==0: newline = [' ']   # do not let blank lines disappear altogether
                endlDoc2 += newline
            endfDoc = endlDoc2  + endfDoc
    else :
        docHeader2 = []
        endfDoc = endfDoc.getLines( )

        # update the documentation, including metadata on first 4 lines:
    if len( [reac for reac in self.reactions if reac.outputChannel.isFission()] ) > 0:
        LFI = True
    else:
        LFI = False
    LRP = -1
    if( self.resonances is not None ) :
        LRP = 0
        if( self.resonances.reconstructCrossSection ) :
            LRP = 1
        elif( self.resonances.unresolved is not None and self.resonances.resolved is None ) : # self-shielding only
            LRP = 1
        elif( self.resonances.resolved is not None or self.resonances.unresolved is not None ) :
            LRP = 2

    crossSectionScale = self.reactions[0].domainUnitConversionFactor( 'eV' )
    EMAX = max( [ crossSectionScale * reaction.crossSection.domainMax for reaction in self.reactions ] )

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

    NLIB = kwargs.get( 'NLIB', NLIB )

    STA = 0
    try :
        target = self.PoPs[targetID]
    except :
        target = self.PoPs.chemicalElements.getSymbol( targetID )
    if( isinstance( target, nuclidePoPsModule.particle ) ) :
        if( len( target.nucleus.halflife ) > 0 ) :
            if( target.nucleus.halflife[0].value == halflifePoPsModule.UNSTABLE ) : STA = 1
    if( levelIndex > 0 ) : STA = 1

    projectileMass = targetInfo['massTracker'].getMassAWR( projectileZA, asTarget = False )
    docHeader = [ endfFormatsModule.endfHeadLine( targetZA, targetInfo['mass'], LRP, LFI, NLIB, NMOD ),
            endfFormatsModule.endfHeadLine( levelEnergy_eV, STA, levelIndex, targetInfo['LISO'], 0, NFOR ),
            endfFormatsModule.endfHeadLine( projectileMass, EMAX, LREL, 0, NSUB, NVER ),
            endfFormatsModule.endfHeadLine( temperature, 0, LDRV, 0, len( docHeader2 + endfDoc ), -1 ) ]
    new_doc = documentationModule.documentation( 'endf', '\n'.join( docHeader + docHeader2 + endfDoc ) )
    endfMFList[1][451] += endfFormatsModule.toEndfStringList( new_doc )

    endfFormatsModule.useRedsFloatFormat = _useRedsFloatFormat

    return( endfFormatsModule.endfMFListToFinalFile( endfMFList, MAT, lineNumbers = lineNumbers ) )

reactionSuiteModule.reactionSuite.toENDF6 = toENDF6

def addDecayGamma( PoPs, particle, baseMT, endfMFList, flags, targetInfo ) :

    MF = 12
    LP = 0
    MT = baseMT + particle.index
    gammaData = []
    levelEnergy_eV = particle.energy[0].float( 'eV' )
    for decayMode in particle.decayData.decayModes :
        IDs = [ product.pid for decay in decayMode.decayPath for product in decay.products ]
        IDs.remove( IDsPoPsModule.photon )
        if( len( IDs ) != 1 ) : raise Exception( 'Do not know how to handle this.' )
        probability = decayMode.probability[0].value
        finalEnergy_eV = PoPs[IDs[0]].energy[0].float( 'eV' )
        _data = [finalEnergy_eV, probability]
        if decayMode.photonEmissionProbabilities:
            _data.append( decayMode.photonEmissionProbabilities[0].value )
        gammaData.append( _data )

    gammaData.sort( reverse = True )
    nGammas = len( gammaData )
    LGp = len( gammaData[0] )
    endfMFList[MF][MT] = [ endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], 2, LGp - 1, MT - baseMT, 0 ),
        endfFormatsModule.endfHeadLine( levelEnergy_eV, 0., LP, 0, LGp * nGammas, nGammas ) ]

    endfMFList[MF][MT] += endfFormatsModule.endfNdDataList( gammaData )
    endfMFList[MF][MT].append( endfFormatsModule.endfSENDLineNumber( ) )

        # Currently, assume all distributions are isotropic
    endfMFList[14][MT] = [ endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'], 1, 0, nGammas, 0 ) ]
    endfMFList[14][MT].append( endfFormatsModule.endfSENDLineNumber( ) )
