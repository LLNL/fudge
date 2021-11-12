# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from brownies.legacy.toENDF6 import endfFormats as endfFormatsModule
from brownies.legacy.toENDF6 import gndsToENDF6 as gndsToENDF6Module

from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import \
    CoulombPlusNuclearElastic as CPNElasticModule, nuclearAmplitudeExpansion as nuclearAmplitudeExpansionModule, \
    RutherfordScattering as RutherfordScatteringModule
from fudge.reactionData.doubleDifferentialCrossSection.chargedParticleElastic import \
    nuclearPlusInterference as nuclearPlusInterferenceModule
from fudge.productData.distributions import angular as angularModule

#
# CoulombPlusNuclearElastic.form
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    targetInfo['LIDP'] = self.identicalParticles
    targetInfo['productFrame'] = self.productFrame
    self.data.toENDF6( MT, endfMFList, flags, targetInfo )
    del targetInfo['LIDP'], targetInfo['productFrame']

CPNElasticModule.form.toENDF6 = toENDF6

#
# nuclearAmplitudeExpansion
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    def LTP_oneSubParsing( LTP, LIDP, nuclear, interferenceReal, interferenceImaginary, lineData ) :

        legendreDat = nuclear.coefficients.copy( )
        if LIDP:
            legendreDat = legendreDat[::2]  # ENDF doesn't store odd-L coefficients
            NL = len( legendreDat ) - 1
            NW = 3 * NL + 3
        else :
            NL = ( len( legendreDat ) - 1 ) // 2
            NW = 4 * NL + 3
        lineData.append( endfFormatsModule.endfContLine( 0, nuclear.outerDomainValue, LTP, 0, NW, NL ) )

        for j, r in enumerate( interferenceReal.coefficients ) :
            legendreDat.append( r )
            legendreDat.append( interferenceImaginary[j] )
        lineData += endfFormatsModule.endfDataList( legendreDat )

    counts, interpolationFlagsList, lineData = 0, [], []
    LTP = 1                     # indicates this is a nuclear + interference section

    reactionSuite = targetInfo['reactionSuite']
    projectile = reactionSuite.PoPs[reactionSuite.projectile]
    if hasattr(projectile, 'nucleus'): projectile = projectile.nucleus
    projectileSpin = 0
    if( len( projectile.spin ) > 0 ) : projectileSpin = projectile.spin[0].value
    LIDP = targetInfo['LIDP']
    if( isinstance( self.nuclearTerm.data, angularModule.XYs2d ) ) :
        for ridx in range( len( self.nuclearTerm.data ) ) :
            counts += 1
            nuclear, interferenceReal, interferenceImaginary = (self.nuclearTerm.data[ridx],
                                                                self.realInterferenceTerm.data[ridx],
                                                                self.imaginaryInterferenceTerm.data[ridx] )
            LTP_oneSubParsing( LTP, LIDP, nuclear, interferenceReal, interferenceImaginary, lineData )
        interpolationFlagsList += [ counts, gndsToENDF6Module.gndsToENDFInterpolationFlag( self.nuclearTerm.data.interpolation ) ]
    elif( isinstance( self.nuclearTerm.data, angularModule.regions2d ) ) :
        for regionIndex, region in enumerate( self.nuclearTerm.data ) :
            interferenceReal, interferenceImaginary = ( self.realInterferenceTerm.data[regionIndex],
                                                        self.imaginaryInterferenceTerm.data[regionIndex] )
            for energyIndex, nuclear in enumerate( region ) :
                if( ( regionIndex != 0 ) and ( energyIndex == 0 ) ) : continue
                counts += 1
                LTP_oneSubParsing( LTP, LIDP, nuclear, interferenceReal[energyIndex], interferenceImaginary[energyIndex], lineData )
            interpolationFlagsList += [ counts, gndsToENDF6Module.gndsToENDFInterpolationFlag( region.interpolation ) ]
    else :
        raise NotImplementedError( "Unknown data storage inside nuclearAmplitudeExpansion: %s" % type(self.nuclearTerm.data) )
    interpolationFlags = endfFormatsModule.endfInterpolationList( interpolationFlagsList )
    ENDFDataList = [ endfFormatsModule.endfContLine( projectileSpin, 0, LIDP, 0,
            len( interpolationFlagsList ) / 2, counts ) ] + interpolationFlags + lineData
    LAW = 5
    gndsToENDF6Module.toENDF6_MF6(MT, endfMFList, flags, targetInfo, LAW, targetInfo['productFrame'], ENDFDataList)

nuclearAmplitudeExpansionModule.nuclearAmplitudeExpansion.toENDF6 = toENDF6

#
# nuclearPlusInterference
#
def toENDF6( self, MT, endfMFList, flags, targetInfo ) :

    LI, LTT, MF6 = self.distribution.data.toENDF6( flags, { 'doMF4AsMF6' : True } )
    reactionSuite = targetInfo['reactionSuite']
    projectile = reactionSuite.PoPs[reactionSuite.projectile]
    if hasattr(projectile, 'nucleus'): projectile = projectile.nucleus
    projectileSpin = 0
    if( len( projectile.spin ) > 0 ) : projectileSpin = projectile.spin[0].value
    LIDP = targetInfo['LIDP']
    NR, NE = [ int(a) for a in MF6[0].split()[-2:] ]
    MF6[0] = endfFormatsModule.endfContLine( projectileSpin, 0, LIDP, 0, NR, NE )
    LAW = 5
    gndsToENDF6Module.toENDF6_MF6( MT, endfMFList, flags, targetInfo, LAW, targetInfo['productFrame'], MF6 )

nuclearPlusInterferenceModule.nuclearPlusInterference.toENDF6 = toENDF6

#
# RutherfordScattering (write back as nuclearPlusInterference with 0 cross section / isotropic distribution data)
#
def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent ) :

    MT=2
    NR, NE = 1, 2
    LAW = 5

    reactionSuite = targetInfo['reactionSuite']
    projectile = reactionSuite.PoPs[reactionSuite.projectile]
    if hasattr(projectile, 'nucleus'): projectile = projectile.nucleus
    ZAP = 1000 * projectile.Z + projectile.A
    AWP = targetInfo['massTracker'].getMassAWR(ZAP)
    LIDP = reactionSuite.projectile == reactionSuite.target
    projectileSpin = 0
    if len( projectile.spin ) > 0: projectileSpin = projectile.spin[0].value
    domain = reactionSuite.styles.getEvaluatedStyle().projectileEnergyDomain

    # add 0-cross section data to MF3
    endfMFList[3][MT] = [
        endfFormatsModule.endfHeadLine(targetInfo['ZA'], targetInfo['mass'], 0, 0, 0, 0),
        endfFormatsModule.endfContLine(0, 0, 0, 0, 1, 2),
        endfFormatsModule.endfInterpolationLine([2,2])
    ]
    endfMFList[3][MT] += endfFormatsModule.endfDataList([domain.min, 0, domain.max, 0])
    endfMFList[3][MT].append(endfFormatsModule.endfSENDLineNumber())

    # add isotropic distribution to MF=6, LAW=5, LTP=12
    endfMFList[6][MT] = [
        endfFormatsModule.endfHeadLine(targetInfo['ZA'], targetInfo['mass'], 0, 2, 1, 0),       # 2: COM frame, 1: 1 subsection
        endfFormatsModule.endfContLine(ZAP, AWP, 0, LAW, NR, NE),
        endfFormatsModule.endfInterpolationLine([2,2])
    ]
    endfMFList[6][MT] += endfFormatsModule.endfDataList([domain.min, 1, domain.max, 1])         # constant multiplicity
    endfMFList[6][MT] += [
        endfFormatsModule.endfContLine(projectileSpin, 0, LIDP, 0, NR, NE),
        endfFormatsModule.endfInterpolationLine([2,2])
    ]
    LTP = 12
    for energy in (domain.min, domain.max):
        endfMFList[6][MT].append( endfFormatsModule.endfContLine(0, energy, LTP, 0, 4, 2) )
        endfMFList[6][MT] += endfFormatsModule.endfDataList([-1, 0.5, 1, 0.5])
    endfMFList[6][MT].append(endfFormatsModule.endfSENDLineNumber())

RutherfordScatteringModule.RutherfordScattering.toENDF6 = toENDF6
