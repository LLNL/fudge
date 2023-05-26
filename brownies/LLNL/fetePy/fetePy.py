#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os, sys
from argparse import ArgumentParser

from xData import enums as xDataEnumsModule
from xData import XYs1d as XYs1dModule

from PoPs import IDs as IDsPoPsModule
from PoPs.chemicalElements import misc as chemicalElementMiscPoPsModule
from PoPs.families import nuclide as nuclideModule

from fudge import outputChannel as outputChannelModule
from brownies.legacy.endl import fudgeParameters as fudgeParametersModule

fudgeParametersModule.VerboseMode = 0

from brownies.legacy.converting import endfFileToGNDS as endfFileToGNDSModule, endf_endl as endf_endlModule
from brownies.legacy.endl import endlZA as endlZAClass
from brownies.legacy.endl import bdfls as bdflsModule, endlmisc as endlmiscModule

from fudge import enums as enumsModule
from fudge.productData import multiplicity as multiplicityModule
from fudge.productData.distributions import angular as angularModule
from fudge.productData.distributions import energy as energyModule
from fudge.productData.distributions import uncorrelated as uncorrelatedModule
from fudge.productData.distributions import unspecified as unspecifiedModule

outputDefault = 'ascii'
EMaxDefault = 20
temperatureDefault = 2.58522e-08

energyEps = 2e-7
fractionDefault = 1e-6
style = 'eval'

CSQList = {}

usage = \
"""
Translates gamma data in an ENDF file into ENDL yo07* and *c55* files. If MT is given, only translate that specified MT.
"""

parser = ArgumentParser( description = usage )
parser.add_argument( 'inputFile', type = str, 
        help = 'ENDF-6 file whose gamma data are translated to ENDL equivalent data/files.' )
parser.add_argument( '-MT', type = int, nargs = '?', default = None,
        help = 'Only process reaction gammas for specified MT.' )
parser.add_argument( '-d', '--date', default = None, action = 'store', type = int,
        help = '''The date to use for the ENDL date. Must be of the format yymmdd. If not present, today's date is used.''' )
parser.add_argument( '-e', '--ENDF', default = False, action = 'store_true',
        help = '''Write the original ENDF file to the target's directory.''' )
parser.add_argument( '-g', '--GNDS', default = False, action = 'store_true',
        help = '''Write an GNDS/XML file representing the ENDF file to the target's directory.''' )
parser.add_argument( '-o', '--output', type = str, default = outputDefault, action = 'store',
        help = '''Write the target's files (i.e., (yi0?/za??????/*) to the specified directory. Default = "%s".''' % outputDefault )
parser.add_argument( '-s', '--skipBadData', default = False, action = 'store_true',
        help = 'Skips format errors if possible when translating the ENDF file to GNDS.' )
parser.add_argument( '--MT2Skip', type = int, action = 'append',
        help = 'When translating from ENDF to GNDS, add MT to list of MTs to skip.' )
parser.add_argument( '-c', '--continuumSpectraFix', action = 'store_true',
        help = 'When translating from ENDF to GNDS, set the "continuumSpectraFix" flag to True.' )
parser.add_argument( '-b', '--bdfls', action = 'store', default = None,
        help = 'The next argument specifies the bdfls file to use.' )
parser.add_argument( '-t', '--temperature', action = 'store', type = float, default = temperatureDefault,
        help = 'Sets the temperature in each ENDL dataset to temperature. Default = %s.' % temperatureDefault )
parser.add_argument( '--EMax', action = 'store', type = float, default = EMaxDefault,
        help = 'Only data for gammas with non-zero multiplicity below EMax are written to an ENDL file. Default = %s.' % EMaxDefault )
parser.add_argument( '-f', '--fraction', action = 'store', type = float, default = fractionDefault,
        help = 'If at incident energy E_i the total multiplicity is zero, all gamma multiplcities are evaluated E_i + dE' +
                ' where dE = faction * DE. DE is either ( E_{i+1} - E_i ) or if that multiplicity is also zero, it is ( E_i - E_{i-1} ).' +
                ' Default = %s.' % fractionDefault )
parser.add_argument( '-v', '--verbose', default = 0, action = 'count',
        help = 'Enable verbose output.' )

args = parser.parse_args( )

if( args.MT2Skip is None ) : args.MT2Skip = []

bdflsFile = args.bdfls
if( bdflsFile is None ) :
    bdflsFile = os.path.join( args.inputFile.split( '/endf/' )[0], 'bdfls' )
    if( not( os.path.exists( bdflsFile ) ) ) : bdflsFile = None
if( bdflsFile is not None ) : bdflsFile = bdflsModule.bdfls( template = bdflsFile )

rce = endfFileToGNDSModule.endfFileToGNDS( args.inputFile, singleMTOnly = args.MT, MTs2Skip = args.MT2Skip,
        toStdOut = args.verbose > 2, skipBadData = args.skipBadData, doCovariances = False, verbose = 0,
        printBadNK14 = False, continuumSpectraFix = args.continuumSpectraFix, acceptBadMF10FissionZAP = True )
reactionSuite = rce['reactionSuite']
reactionSuite.convertUnits( { 'eV' : 'MeV' } )

yi = endlmiscModule.incidentParticleTags( str( reactionSuite.projectile ) )[1]

suffix = ''
ELevel = 0
if( reactionSuite.target in reactionSuite.PoPs.aliases ) :
    target = reactionSuite.PoPs.aliases[reactionSuite.target]
    suffix = 'm%s' % target.metaStableIndex
    target = reactionSuite.PoPs[target.pid]
else :
    target = reactionSuite.PoPs[reactionSuite.target]
ZA = chemicalElementMiscPoPsModule.ZA( target )
if( isinstance( target, nuclideModule.Particle ) ) : ELevel = target.nucleus.energy[0].float( 'MeV' )

endlZA = endlZAClass( ZA, yi, suffix = suffix, workDir = args.output, bdflsFile = bdflsFile )
if( args.verbose > 0 ) : print( 'Using bdfls files %s' % endlZA.bdflsFile.template )

workDir = endlZA.workDir
os.system( 'rm -rf %s/*' % endlZA.workDir )

if( args.ENDF ) : os.system( 'cp %s %s' % ( args.inputFile, workDir ) )
if( args.GNDS ) : reactionSuite.saveToFile( os.path.join( workDir, 'gnds.xml' ) )

def getMultiplicityForDistributionSum( self, energy, energies, numberOfGammas, zeroTotal ) :

    if( numberOfGammas == 0 ) : return( 1 )
    multiplicity = self.multiplicity.evaluate( energy )
    if( ( multiplicity == 0 ) and zeroTotal ) :
        found = False
        for i2, e2 in enumerate( energies ) :
            if( found ) : break
            if( e2 == energy ) : found = True
        ePlus = energy + args.fraction * ( e2 - energy )
        multiplicity = self.multiplicity.evaluate( ePlus )
        if( ( multiplicity is None ) and ( i2 > 1 ) ) :
            ePlus = energy - args.fraction * ( energy - energies[i2-2] )
            multiplicity = self.multiplicity.evaluate( ePlus )
    return( multiplicity )

def gammaDeltaDistribution( energy ) :

    deltaEnergy = float( "1e%s" % ("%.0e" % ( energy * energyEps )).split( 'e' )[1] )
    eMin, eMax = energy - deltaEnergy, energy + deltaEnergy
    if( eMin < 0 ) : return( XYs1dModule.XYs1d( [ [ energy, 1. ], [ eMax, 0. ] ] ).normalize( ) )
    return( XYs1dModule.XYs1d( [ [ eMin, 0. ], [ energy, 1 ], [ eMax, 0. ] ] ).normalize( ) )

def finalDistributionCheck( fullDistribution ) :

    for i1, ( Ein, PofEout ) in enumerate( fullDistribution  ) :
        for i2, PEout2 in enumerate( PofEout ) :
            if( i2 > 0 ) :
                if( ( PEout2[0] - PEout1[0] ) < 0.1 * energyEps * PEout2[0] ) : PEout2[0] = ( 1 + 0.1 * energyEps ) * PEout1[0]
            PEout1 = PEout2
        PEout = XYs1dModule.XYs1d( fullDistribution[0][1] )
        if( abs( PEout.integrate( ) - 1 ) > 1e-5 ) : fullDistribution[i1][1] = PEout.normalize( ).copyDataToXYs( )

class PrimaryGamma :

    def __init__( self, multiplicity, energy, massRatio, angularForm ) :

        self.multiplicity = multiplicity
        self.energy = energy
        self.massRatio = massRatio
        self.angularForm = angularForm

    def getDistribution( self, energy, energies, numberOfGammas, zeroTotal ) :
        """Returns delta function gamma spectrum."""

        if( self.multiplicity.domainMin <= energy <= self.multiplicity.domainMax ) :
            multiplicity = getMultiplicityForDistributionSum( self, energy, energies, numberOfGammas, zeroTotal )
            return( multiplicity * gammaDeltaDistribution( self.energy + self.massRatio * energy ) )
        else :
            return( None )

class DiscreteGamma :

    def __init__( self, multiplicity, energy, angularForm ) :

        self.multiplicity = multiplicity
        self.energy = energy
        self.angularForm = angularForm

    def getDistribution( self, energy, energies, numberOfGammas, zeroTotal ) :
        """Returns delta function gamma spectrum."""

        if( self.multiplicity.domainMin <= energy <= self.multiplicity.domainMax ) :
            multiplicity = getMultiplicityForDistributionSum( self, energy, energies, numberOfGammas, zeroTotal )
            return( multiplicity * gammaDeltaDistribution( self.energy ) )
        else :
            return( None )

class ContinuumGammaIsotropic :

    def __init__( self, multiplicity, energyForm ) :

        self.multiplicity = multiplicity
        self.energyForm = energyForm

    def getDistribution( self, energy, energies, numberOfGammas, zeroTotal ) :

        if( self.multiplicity.domainMin <= energy <= self.multiplicity.domainMax ) :
            domainMin, domainMax = self.energyForm.domainMin, self.energyForm.domainMax
            dE = energy - domainMin
            if( dE < 0 ) :
                if( abs( dE ) < 1e-15 * energy ) :
                    energy = domainMin
                else :
                    return( None )
            dE = domainMax - energy
            if( dE < 0 ) :
                if( abs( dE ) < 1e-15 * energy ) :
                    energy = domainMin
                else :
                    return( None )

            energyForm = self.energyForm
            if( isinstance( energyForm, ( energyModule.Regions2d, ) ) ) :
                for region in energyForm :
                    if( energy < region.domainMax ) : break
                energyForm = region

            PofEp1 = PofEp2 = energyForm[0]
            for PofEp2 in energyForm :
                if( PofEp2.value >= energy * ( 1 - 1e-15 ) ) : break
                PofEp1 = PofEp2
            if(   abs( energy - PofEp1.value ) < ( energy * 1e-3 ) ) :            # If energy is close to PofEp1.value use PofEp1
                energyForm = PofEp1.toPointwise_withLinearXYs( accuracy = 1e-4, lowerEps = energyEps, upperEps = energyEps )
            elif( abs( energy - PofEp2.value ) < ( energy * 1e-3 ) ) :            # If energy is close to PofEp2.value use PofEp2
                energyForm = PofEp2.toPointwise_withLinearXYs( accuracy = 1e-4, lowerEps = energyEps, upperEps = energyEps )
            else :
                fraction2 = ( energy - PofEp1.value ) / ( PofEp2.value - PofEp1.value )
                fraction2 = min( 1.0, max( 0.0, fraction2 ) )
                energyForm1 = PofEp1.toPointwise_withLinearXYs( accuracy = 1e-4, lowerEps = energyEps, upperEps = energyEps )
                energyForm2 = PofEp2.toPointwise_withLinearXYs( accuracy = 1e-4, lowerEps = energyEps, upperEps = energyEps )
                try :
                    energyForm = XYs1dModule.pointwiseXY_C.unitbaseInterpolate( energy, PofEp1.value, energyForm1, PofEp2.value, energyForm2, 1 )
                except :
                    if( fraction2 > 0.5 ) :
                        energyForm = energyForm2
                    else :
                        energyForm = energyForm1
            lowerEps = energyEps
            if( energyForm.domainMin == 0 ) : lowerEps = 0
            energyForm = energyForm.dullEdges( lowerEps = lowerEps, upperEps = energyEps )
            energyForm = XYs1dModule.XYs1d( energyForm )
            multiplicity = getMultiplicityForDistributionSum( self, energy, energies, numberOfGammas, zeroTotal )
            distribution = multiplicity * energyForm
            return( distribution )
        else :
            return( None )

def processGamma( gamma, gammas, crossSection ) :

    distribution = gamma.distribution[style]
    if( isinstance( distribution, unspecifiedModule.Form ) ) :
        if( args.verbose > 1 ) : print( '            No distribution data' )
        return

    MF13 = False
    if( 'ENDFconversionFlag' in gamma.attributes ) : MF13 = gamma.attributes['ENDFconversionFlag'] == 'MF13'
    multiplicity = gamma.multiplicity[style]
    if(   isinstance( multiplicity, ( multiplicityModule.Constant1d ) ) ) :
        domainMin, domainMax = multiplicity.domainMin, multiplicity.domainMax
        multiplicity = XYs1dModule.XYs1d( [ [ float( domainMin ), multiplicity.value ], [ float( domainMax ), multiplicity.value ] ],
                axes = XYs1dModule.XYs1d.defaultAxes( labelsUnits = { 1 : ( 'energy_in', 'MeV' ) } ) )
    elif( isinstance( multiplicity, ( multiplicityModule.XYs1d, multiplicityModule.Regions1d ) ) ) :
        if( isinstance( multiplicity, multiplicityModule.Regions1d ) and MF13 ) :
            for region in multiplicity :
                if region.interpolation == xDataEnumsModule.Interpolation.flat:
                    print( '    WARNING: This may need to be fixed.' )
        if( MF13 and isinstance( multiplicity, multiplicityModule.XYs1d ) ) :
            pass
        else :
            multiplicity = multiplicity.toPointwise_withLinearXYs( lowerEps = energyEps, upperEps = energyEps )
    else :
        raise Exception( 'Unsupported multiplicity form "%s"' % multiplicity.moniker )
    multiplicity = multiplicity.domainSlice( domainMax = args.EMax )
    if( ( len( multiplicity ) < 2 ) or ( multiplicity.rangeMax == 0 ) ) :
        if( args.verbose ) : print( '''INFO: 1) not writing reaction's gamma data as gamma multiplicity is 0. below EMax = %s''' % args.EMax )
        return

    if( crossSection is not None ) :
        if( MF13 ) :
            values = [ [ E1, m1 * crossSection.evaluate( E1 ) ] for E1, m1 in multiplicity ]
            multiplicity = XYs1dModule.XYs1d( values, axes = multiplicity.axes, interpolation = multiplicity.interpolation )
            multiplicity = multiplicity.toPointwise_withLinearXYs( lowerEps = energyEps, upperEps = energyEps )
            multiplicity.axes[0].unit = "b"
        else :
            crossSection = crossSection.domainSlice( domainMax = args.EMax )
            crossSection.setInfill( False )
            multiplicity.setInfill( False )
            try :
                multiplicity = multiplicity * crossSection
            except :
                print( 'multiplicity' )
                print( multiplicity.toString( ) )
                print( 'crossSection' )
                print( crossSection.toString( ) )
                print( multiplicity.domainMin, multiplicity.domainMax, crossSection.domainMin, crossSection.domainMax )
                raise
            multiplicity = multiplicity.trim( )

    if( isinstance( distribution, uncorrelatedModule.Form ) ) :
        angularForm = distribution.angularSubform.data
        energyForm = distribution.energySubform.data
        if( isinstance( angularForm, angularModule.XYs2d ) ) :
            if( not( angularForm.isIsotropic( ) ) ) : print( '    WARNING: treating angular XYs2d as isotropic' )
            isIsotropic = True
        else :
            isIsotropic = angularForm.isIsotropic( )
            if( not( isIsotropic ) ) : raise Exception( 'unsupported angular distribution' )
        if( isIsotropic ) :
            if( isinstance( energyForm, energyModule.PrimaryGamma ) ) :
                gammas['primaries'].append( PrimaryGamma( multiplicity, energyForm.value,
                        float( energyForm.massRatio ), angularForm ) )
            elif( isinstance( energyForm, energyModule.DiscreteGamma ) ) :
                gammas['discretes'].append( DiscreteGamma( multiplicity, energyForm.value, angularForm ) )
            else :
                if( isinstance( energyForm, ( energyModule.XYs2d, energyModule.Regions2d ) ) ) :
                    pass
                else :
                    raise Exception( 'Unsupported energy form "%s"' % energyForm.moniker )
                gammas['continuum'].append( ContinuumGammaIsotropic( multiplicity, energyForm ) )
        else :
            raise Exception( 'Unsupported angular form "%s"' % angularForm.moniker )
    else :
        raise Exception( 'Unsupported distribution type "%s"' % distribution.moniker )

def processChannel( channel, gammas, branchingGammas, crossSection = None ) :

    numberOfGammasProcessed = 0
    for product in channel :
        if( product.pid in branchingGammas ) :
            gammas['branchingGammas'].append( ( product.pid, branchingGammas[product.pid] ) )
        else :
            if( product.pid == IDsPoPsModule.photon ) :
                energy = ''
                distribution = product.distribution[style]
                if( isinstance( distribution, uncorrelatedModule.Form ) ) :
                    energyForm = distribution.energySubform.data
                    if( isinstance( energyForm, energyModule.PrimaryGamma ) ) :
                        energy = ": primary energy = %s" % ( energyForm.value )
                    elif( isinstance( energyForm, energyModule.DiscreteGamma ) ) :
                        energy = ": discrete energy = %s" % ( energyForm.value )
                if( args.verbose > 1 ) : print( '        %-40s: %-12s%s' % ( product, product.label, energy ) )
                if channel.genre == enumsModule.Genre.twoBody:
                    print('    WARNING: skipping two body gamma')
                else:
                    processGamma( product, gammas, crossSection )
                    numberOfGammasProcessed += 1
            if( product.outputChannel is not None ) :
                numberOfGammasProcessed += processChannel( product.outputChannel, gammas, branchingGammas, crossSection = crossSection )
    return( numberOfGammasProcessed )

def getDiscretePrimaryMultiplicity( gammas, multiplicity, gammaList ) :

    for gamma in gammas['discretes'] + gammas['primaries'] :
        gammaList.append( gamma )
        if( multiplicity is None ) :
            multiplicity = gamma.multiplicity.copy( )
        else :
            try :
                m1, m2 = multiplicity.mutualify(    lowerEps1 = energyEps, upperEps1 = energyEps, positiveXOnly1 = 1, 
                        other = gamma.multiplicity, lowerEps2 = energyEps, upperEps2 = energyEps, positiveXOnly2 = 1 )
                multiplicity = m1 + m2
            except :
                print( multiplicity.domainMin, multiplicity.domainMax )
                print( gamma.multiplicity.domainMin, gamma.multiplicity.domainMax )
                raise

    return( multiplicity )

def sumDistributions( multiplicity, gammaList ) :

    fullDistribution = []
    energies = [ energy for energy, m1 in multiplicity ]
    for energy, m1 in multiplicity :
        distributionSum = None
        zeroTotal = True
        for gamma in gammaList :
            value = gamma.multiplicity.evaluate( energy )
            if( value is None ) : continue
            if( value > 0 ) : zeroTotal = False
        for gamma in gammaList :
            distribution = gamma.getDistribution( energy, energies, len( gammaList ), zeroTotal )
            if( distribution is None ) : continue
            if( distribution.integrate( ) == 0 ) : continue
            if( distributionSum is None ) :
                distributionSum = distribution
            else :
                try :
                    distributionSum = distributionSum + distribution
                except :
                    print( distributionSum.toString( ) )
                    print( distribution.toString( ) )
                    raise
# FIXME, how do we sum distributions when multiplicity is 0.
        if( distributionSum is None ) : continue
        norm = distributionSum.normalize( )
        fullDistribution.append( [ energy, norm.copyDataToXYs( ) ] )
    return( fullDistribution )

def writeGammas( ENDLFiles, C, S, Q, X1, multiplicity, fullDistribution ) :

    trimmed = multiplicity.trim( )
    if( isinstance( trimmed, XYs1dModule.XYs1d ) ) :
        domainMin, rangeMax = trimmed.domainMin, trimmed.rangeMax
    else :
        domainMin, rangeMax = trimmed.domainMin( ), trimmed.rangeMax( )
    if( ( domainMin >= args.EMax ) or ( rangeMax == 0 ) ) :
        if( args.verbose ) : print( '''INFO: 2) not writing reaction's gamma data as gamma multiplicity is 0. below EMax = %s''' % args.EMax )
        return
    
    yo, I, ENDLFileName = 7, 9, "yo%2.2dc%2.2di%3.3ds%3.3d" % ( 7, C, 9, S )
    if( C == 55 ) : yo, I, ENDLFileName = 0, 0, "yo%2.2dc%2.2di%3.3ds%3.3d" % ( 0, C, 0, S )
    if( ENDLFileName not in ENDLFiles ) : ENDLFiles[ENDLFileName] = endlZA.addFile( yo, C, I, S )
    ENDLFile = ENDLFiles[ENDLFileName]
    ENDLFile.addData( multiplicity.copyDataToXYs( ), Q = Q, X1 = X1, temperature = args.temperature )

    if( fullDistribution is not None ) :
        finalDistributionCheck( fullDistribution )
        ENDLFileName = "yo%2.2dc%2.2di%3.3ds%3.3d" % ( 7, C, 4, S )
        if( ENDLFileName not in ENDLFiles ) : ENDLFiles[ENDLFileName] = endlZA.addFile( 7, C, 4, S )
        ENDLFile = ENDLFiles[ENDLFileName]
        ENDLFile.addData( [ [ 0, fullDistribution ] ], Q = Q, X1 = X1, temperature = args.temperature )

def addGammas(multiplicity, ID, photonBranchingData, branchingData):

    for probability, energy, residualID, photonEmissionProbability in photonBranchingData[ID]['photons'] :
        branchingData.append([multiplicity * probability * photonEmissionProbability, energy])
        addGammas(multiplicity * probability, residualID, photonBranchingData, branchingData)

branchingData = {}
photonBranchingData = reactionSuite.photonBranchingData( )
for isotopeID in photonBranchingData :
    for residualID in photonBranchingData[isotopeID] :
        branchingData[residualID] = []
        addGammas( 1.0, residualID, photonBranchingData[isotopeID], branchingData[residualID] )

branchingGammas = {}
for residualID in branchingData :
        multiplicity = 0
        spectra = None
        for probability, gammaEnergy in branchingData[residualID] :
            spectrum = probability * gammaDeltaDistribution( float( gammaEnergy ) )
            data = [ [ float( "%.14e" % x ), y ] for x, y in spectrum ]
            spectrum = XYs1dModule.XYs1d( data, axes = spectrum.axes )
            multiplicity += probability
            if( spectra is None ) :
                spectra = spectrum
            else :
                spectra = spectra + spectrum
        if( spectra is not None ) :
            spectra = spectra.normalize( )
            branchingGammas[residualID] = ( multiplicity, spectra )

if( args.verbose ) : print( '%s ->' % reactionSuite.inputParticlesToReactionString( ) )

multiplicitySums = {}
for _sum in reactionSuite.sums.multiplicitySums :
    multiplicity = _sum.multiplicity.toPointwise_withLinearXYs( lowerEps = energyEps, upperEps = energyEps )
    multiplicitySums[_sum.ENDF_MT] = multiplicity.domainSlice( domainMax = args.EMax )
ENDLFiles = {}
for reaction in reactionSuite.reactions :
    MT = reaction.ENDF_MT
    C, S = endf_endlModule.getCSFromMT( MT )
    if( C < 0 ) :
        crossSection = reaction.crossSection.toPointwise_withLinearXYs( lowerEps = energyEps, upperEps = energyEps )
        print( '    WARNING: skipping MT=%s as it has no ENDL C equivalent. Domain range is %s to %s.' % ( MT, crossSection.domainMin, crossSection.domainMax ) )
        continue
    if( C == 10 ) : continue                    # There are no gammas for elastic scattering, but logic below may produce gammas for 
                                                # meta-stables with branching data so let's skip it now.
    Q = reaction.getQ( unit = 'MeV' )           # BRB, this is the wrong Q.
    X1 = 0.
    if( len( reaction.outputChannel ) == 2 ) :
        particle = reaction.outputChannel[1].particle
        if( isinstance( particle, nuclideModule.Particle ) ) :
            X1 = particle.energy[0].value
    Q += X1                                     # BRB, is this the right Q?
    if( C not in CSQList ) : CSQList[C] = {}
    if( S not in CSQList[C] ) : CSQList[C][S] = Q
    Q2 = CSQList[C][S]
    if( abs( Q2 - Q ) < 1e-7 * abs( Q ) ) : Q = Q2
    if( S not in [ 1, 3 ] ) : X1 = 0
        
    if( args.verbose ) :
        print( '    %-40s: MT = %3s C = %2s S = %s Q = %-12s X1 = %-12s' % ( reaction, MT, C, S, Q, X1 ) )
        sys.stdout.flush( )
    gammas = { 'continuum' : [], 'discretes' : [], 'primaries' : [], 'branchingGammas' : [] }
    processChannel( reaction.outputChannel, gammas, branchingGammas )

    multiplicity = None
    gammaList = []
    if( len( gammas['continuum'] ) > 0 ) :
        if( len( gammas['continuum'] ) > 1 ) : raise Exception( 'multiple continuum gammas no currently supported' )
        continuum = gammas['continuum'][0]
        multiplicity = continuum.multiplicity
        gammaList.append( continuum )

    multiplicity = getDiscretePrimaryMultiplicity( gammas, multiplicity, gammaList )

    if( multiplicity is not None ) :
        if( len( gammas['branchingGammas'] ) != 0 ) : raise Exception( 'Currently, branching gammas are not support with other gammas.' )
        fullDistribution = sumDistributions( multiplicity, gammaList )
    else :
        if( len( gammas['branchingGammas'] ) not in ( 0, 1 ) ) : raise Exception( 'need to support this' )
        for levelName, ( gammaMultiplicity, spectrum ) in gammas['branchingGammas'] :
            domainMin, domainMax = reaction.domainMin, reaction.domainMax
            multiplicity = XYs1dModule.XYs1d( [ [ domainMin, gammaMultiplicity ], [ domainMax, gammaMultiplicity ] ] )
            fullDistribution = [ [ domainMin, spectrum ], [ domainMax, spectrum ] ]

    if( multiplicity is not None ) :
        if( MT in multiplicitySums ) : multiplicity = multiplicitySums[MT]
        writeGammas( ENDLFiles, C, S, Q, X1, multiplicity, fullDistribution )

gammas = { 'continuum' : [], 'discretes' : [], 'primaries' : [] }

productionMT = None
for production in reactionSuite.orphanProducts :
    MT = production.ENDF_MT
    if( args.verbose ) : print( '    %-40s: MT = %3s C = 55' % ( production, MT ) )
    if( productionMT is None ) : MT = productionMT
    C55CrossSection = production.crossSection.toPointwise_withLinearXYs( lowerEps = energyEps, upperEps = energyEps )
    numberOfGammasProcessed = processChannel( production.outputChannel, gammas, {}, crossSection = C55CrossSection )
    if( ( numberOfGammasProcessed > 0 ) and ( MT != productionMT ) ) :
        raise Exception( 'Currently only support production gammas for one MT: %d %d' % ( productionMT, MT ) )

if( len( gammas['discretes'] ) > 0 ) :
    discretes = sorted( [ [ discrete.energy, discrete ] for discrete in gammas['discretes'] ] )
    for energy, discrete in discretes :
        domainMin, domainMax = discrete.multiplicity.domainMin, discrete.multiplicity.domainMax
        angularDistribution = [ [ domainMin, [ [ -1, 0.5 ], [ 1, 0.5 ] ] ], [ domainMax, [ [ -1, 0.5 ], [ 1, 0.5 ] ] ] ]
        writeGammas( ENDLFiles, 55, 3, 0, discrete.energy, discrete.multiplicity, None )
        ENDLFileName = "yo%2.2dc%2.2di%3.3ds%3.3d" % ( 7, 55, 1, 3 )
        if( ENDLFileName not in ENDLFiles ) : ENDLFiles[ENDLFileName] = endlZA.addFile( 7, 55, 1, 3 )
        ENDLFiles[ENDLFileName].addData( angularDistribution, Q = 0, X1 = discrete.energy, temperature = args.temperature )
    gammas['discretes'] = None

multiplicity = None
gammaList = []
for continuumGamma in gammas['continuum'] :
    if( multiplicity is None ) :
        multiplicity = continuumGamma.multiplicity
    else :
        multiplicity = multiplicity + continuumGamma.multiplicity
    gammaList.append( continuumGamma )

if( len( gammas['primaries'] ) != 0 ) : raise Exception( 'need to add primary gammas to continuum' )
if( multiplicity is not None ) :
    fullDistribution = sumDistributions( multiplicity, gammaList )
    writeGammas( ENDLFiles, 55, 0, 0, 0, multiplicity, fullDistribution )

for data in endlZA.findDatas( ) :
    data.setFormat( 15 )
    data.setELevel( ELevel )
    if( args.date is not None ) : data.setDate( args.date )

endlZA.save( )
