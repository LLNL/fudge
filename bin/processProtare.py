#! /usr/bin/env python
 
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

from __future__ import print_function
import sys
import os
import shutil

from argparse import ArgumentParser

from PoPs import IDs as IDsPoPsModule

from xData import axes as axesModule
from xData import values as valuesModule
from xData import XYs as XYsModule
from xData import physicalQuantity as physicalQuantityModule

from pqu import PQU as PQUModule

import fudge.gnds.physicalQuantity as temperatureModule
import fudge.gnds.styles as stylesModule
import fudge.gnds.reactionSuite as reactionSuiteModule
from fudge.legacy.endl import bdfls as bdflsModule

from fudge.processing import flux as fluxModule
from fudge.processing import transportables as transportablesModule
from fudge.processing import group as groupModule

dateTimeStr = None  ### should handle date automagically

lMaxDefault = 9
fidDefault = 'LLNL_fid_1'
bdflsDefault = '/usr/gapps/data/nuclear/bdfls.archive/bdfls.Audi_etal.2003.12.22'
reconAccuracyDefault = 1e-6
CoulombPlusNuclearMuCutOffDefault = 0.94                                            # This is the mu cutoff used for ENDL data.
tempsDefault = '2.586e-8'
tempUnitDefault = 'MeV/k'
energyUnitDefault = 'MeV'
outputDefault = 'SRC'

EVprefixDefault = 'eval'
RCprefixDefault = 'recon'
muCutoffPrefixDefault = 'muCutoff'
AEPprefixDefault = 'apd'
HTprefixDefault = 'heated'
MGprefixDefault = 'MultiGroup'
MCprefixDefault = 'MonteCarlo'
UPprefixDefault = 'UpScatter'

    ### Default group structures
gids = {'n'      : 'LLNL_gid_7',
        'H1'     : 'LLNL_gid_71', 
        'H2'     : 'LLNL_gid_71', 
        'H3'     : 'LLNL_gid_71', 
        'He3'    : 'LLNL_gid_71', 
        'He4'    : 'LLNL_gid_71', 
        'photon' : 'LLNL_gid_70' }
gidDefaults = ""
for gid in gids : gidDefaults += '%s="%s", ' % ( gid, gids[gid] )
gidDefaults = gidDefaults[:-2]
 
description = """Processes all data in a GNDS file."""

parser = ArgumentParser( description = description )
parser.add_argument( "gnds", type = str,                                                                help = 'input gnds file to process')
parser.add_argument( "outputFile", nargs = "?", default = None,                                         help = 'output file' )

parser.add_argument( "--tag", type = str, default = '.proc',                                            help = 'tag to indicate data has been processed')
parser.add_argument( '-o', '--output', default = outputDefault,                                         help = 'directory to write output file. If "SRC", writes to source directory. Default is %s' % outputDefault )
parser.add_argument( "--writeConvertedUnits", action = 'store_true', default = False,                   help = 'write data to file after units converted' )

parser.add_argument( "--energyUnit", type = str, default = energyUnitDefault,                           help = 'energy unit to convert to. Default is %s ' % energyUnitDefault)

parser.add_argument( "--reconstruct", type = str, choices = ['all','crossSection','angular'], default='crossSection', 
                                                                                                        help = "What kind of reconstructing should be done ('all','crossSection','angular')" )
parser.add_argument( "--CoulombPlusNuclearMuCutOff", type = float, default = CoulombPlusNuclearMuCutOffDefault,
                                                                                                        help = 'For Coulomb + nuclear elastic scattering mu is limited to [ -1, muCutOff ] to make the cross section finite. For identical particles mu is limited to [ -muCutOff, muCutOff ]. muCutOff must be in the range ( -1 to 1 ). Default is %s.' % CoulombPlusNuclearMuCutOffDefault )

parser.add_argument( "-t", "--temperatures", type = float, action='append', default = None,             help = "temperatures for heating, use successive -t arguments" )
parser.add_argument( "--temperatureUnit", type = str, default = tempUnitDefault,                        help = 'temperature unit to convert to, default is %s ' % tempUnitDefault )

parser.add_argument( "--bdfls", type = str, default = bdflsDefault,                                     help = "bdfls file to use. Default is %s " % bdflsDefault )
parser.add_argument( "--fluxID", type = str, default = fidDefault,                                      help = 'Flux ID from bdfls. Default is %s ' % fidDefault )
parser.add_argument( "-g", "--gid", type = str, action='append', default = None,                        help = "particle grouping schemes : <particle>=LLNL_gid_<integer>. Default is %s." % gidDefaults )
parser.add_argument( "--legendreMax", type = int, default = lMaxDefault,                                help = "Maximum Legendre order for Sn prcessed data. Default is %s." % lMaxDefault )
parser.add_argument( "-mg", "--MultiGroup", action = 'store_true',                                      help = "flag for MultiGroup processing" )
parser.add_argument( "-up", "--UpScatter", action = 'store_true',                                       help = "flag for Upscatter processing" )

parser.add_argument( "-mc", "--MonteCarlo", action = 'store_true',                                      help = "flag for Monte Carlo processing" )

parser.add_argument( "--threads",type = int, default = 1,                                               help = "number of threads to use for temperatures loop" )
parser.add_argument( "-v", "--verbose", action = "count", default = 0,                                  help = "enable verbose output" )
parser.add_argument( "--reconAccuracy", type = float, default = reconAccuracyDefault,                   help = "Accuracy for reconstructing resonances. Default is %.1e" % reconAccuracyDefault )

parser.add_argument( "--prefixEval", type = str, default=EVprefixDefault,                               help = "prefix for Evaluation styles. Default is %s." % EVprefixDefault )
parser.add_argument( "--prefixRecon", type = str, default=RCprefixDefault,                              help = "prefix for Resonance Reconstruction styles. Default is %s." % RCprefixDefault )
parser.add_argument( '--prefixMuCutoff', type = str, default = muCutoffPrefixDefault,                   help = 'prefix for Coulomb + nuclear elastic scattering mu cutoff style. Default is "%s".' % muCutoffPrefixDefault )
parser.add_argument( "--prefixAep", type = str, default=AEPprefixDefault,                               help = "prefix for Average Energy to Product styles. Default is %s." % AEPprefixDefault )
parser.add_argument( "--prefixHeated", type = str, default=HTprefixDefault,                             help = "prefix for heated styles. Default is %s." % HTprefixDefault )
parser.add_argument( "--prefixMultiGroup", type = str, default=MGprefixDefault,                         help = "prefix for MultiGroup styles. Default is %s." % MGprefixDefault )
parser.add_argument( "--prefixMonteCarlo", type = str, default=MCprefixDefault,                         help = "prefix for MonteCarlo styles. Default is %s." % MCprefixDefault )
parser.add_argument( "--prefixUpscatter", type = str, default=UPprefixDefault,                          help = "prefix for UpScatter styles. Default is %s." % UPprefixDefault )

args = parser.parse_args( )

logFile = open( 'logFile' + args.tag, 'w' )

bdflsFile = args.bdfls 
bdfls = bdflsModule.getDefaultBdfls( template = bdflsFile )

### get flux data
bdflsFlux = bdfls.flux( args.fluxID )
if( len( bdflsFlux ) != 1 ) : raise Exception( 'flux order greater than 0 currently not supported.' )
axes = axesModule.axes( rank = 3 )
axes[2] = axesModule.axis( 'energy_in', 2, args.energyUnit )
axes[1] = axesModule.axis( 'mu', 1, '' )
axes[0] = axesModule.axis( 'flux(energy_in,mu)', 0, '1/s' )
fluxData = fluxModule.XYs2d( axes = axes )
for energy_MeV, C0 in bdflsFlux.EF_l[0] :
    newE = physicalQuantityModule.physicalQuantity( energy_MeV, 'MeV' ).getValueAs(args.energyUnit)
    fluxData.append( fluxModule.LegendreSeries( [ C0 ], value = newE, axes = axes ) )
flux = fluxModule.flux( args.fluxID, fluxData )

if not ( args.gid == None ) :
    for gidVal in args.gid:
        key, value = gidVal.split('=')
        if key not in gids: raise UserWarning("particle grouping specified for unknown particle name : %s " % (key) )
        gids[key] = value


transportables = []
for particle in gids :
    gbs = valuesModule.values( [ PQUModule.PQU( boundary, 'MeV' ).getValueAs( args.energyUnit ) for boundary in bdfls.group(gids[particle]).gb ] )
    grid = axesModule.grid( 'energy_in', 0, args.energyUnit, axesModule.boundariesGridToken, gbs )
    group = groupModule.group( gids[particle], grid )
    transportables.append( transportablesModule.transportable( particle, transportablesModule.conserve.number, group ) )
        
#### read in GNDS file
reactionSuite = reactionSuiteModule.readXML( args.gnds )

if( reactionSuite.projectile != IDsPoPsModule.neutron ) : tempsDefault = 0
if( args.temperatures == None ) : args.temperatures = [ tempsDefault ]
args.temperatures.sort( )

### fail on detection of existing processed data
for style in reactionSuite.styles :
    if( isinstance( style, ( stylesModule.heated, stylesModule.griddedCrossSection, stylesModule.multiGroup, stylesModule.griddedCrossSection ) ) ) : 
        raise Exception( "File already contains processed data. Please use a clean file!" )

### convert units if necessary
reactionSuite.convertUnits( {'MeV': args.energyUnit, 'eV': args.energyUnit } )
if( args.writeConvertedUnits ) :
    convertedUnitsFileName = args.outputFile
    if( convertedUnitsFileName is None ) : convertedUnitsFileName = args.gnds
    convertedUnitsFileName = convertedUnitsFileName + '.convertedUnits'
    if( outputDefault != args.output ) : convertedUnitsFileName = os.path.join( args.output, os.path.basename( convertedUnitsFileName ) )
    reactionSuite.saveToFile( convertedUnitsFileName, xs_pdf_cdf1d_singleLine = True )

try :
    evalStyle = reactionSuite.styles[args.prefixEval]
except :
    raise Exception( "Cannot find the requested evaluation label! : %s != %s" % ( evalStyle.label, args.prefixEval ) )
baseStyle = evalStyle

### reconstruct resonances
if( reactionSuite.supportsResonanceReconstruction() ) :
### FIXME: need to add logic about the choice of reconstruction types (different style for reconstructed angular vs crossSection)
    if( args.verbose > 0 ) : print('Processing resonances')
    reconStyle = reactionSuite.styles.getStyleOfClass(stylesModule.crossSectionReconstructed)
    if reconStyle == None: 
        reconStyle = stylesModule.crossSectionReconstructed( args.prefixRecon, evalStyle.label, date = dateTimeStr )
        reactionSuite.styles.add( reconStyle )
        reactionSuite.reconstructResonances( reconStyle )   ### FIXME need logfile arguments in this process call
    baseStyle = reconStyle

CoulombPlusNuclearMuCutoffs = reactionSuite.CoulombPlusNuclearMuCutoffs( )
if( CoulombPlusNuclearMuCutoffs is not None ) :
    if( args.CoulombPlusNuclearMuCutOff not in CoulombPlusNuclearMuCutoffs ) :
        muCutoffStyle = stylesModule.CoulombPlusNuclearElasticMuCutoff( args.prefixMuCutoff, baseStyle.label, args.CoulombPlusNuclearMuCutOff, date = dateTimeStr )
        reactionSuite.styles.add( muCutoffStyle )
        reactionSuite.processCoulombPlusNuclearMuCutoff( muCutoffStyle )
        baseStyle = muCutoffStyle

### Calculate Average Product Data  ( Energy and Momenta )
apdStyle = reactionSuite.styles.getStyleOfClass( stylesModule.averageProductData )
if apdStyle == None: 
    if( args.verbose > 0 ) : print('Processing average product data')
    apdStyle = stylesModule.averageProductData( args.prefixAep, baseStyle.label, date = dateTimeStr )
    reactionSuite.styles.add( apdStyle )
    reactionSuite.calculateAverageProductData( apdStyle, indent = '  ', verbosity = args.verbose - 2 )  ### FIXME need logfile arguments in this process call
preLoopStyle = apdStyle

if( args.MonteCarlo ) :
    if( args.verbose > 0 ) : print('Processing Monte Carlo')
    MonteCarloStyle = stylesModule.MonteCarlo_cdf( '%s' % ( args.prefixMonteCarlo ), apdStyle.label, date = dateTimeStr )
    reactionSuite.styles.add( MonteCarloStyle )
    reactionSuite.processMC_cdf( MonteCarloStyle, verbosity = args.verbose - 1 )
    preLoopStyle = MonteCarloStyle
        
if( args.MultiGroup ) :
    multiGroupStyle = stylesModule.multiGroup( '%s' % ( args.prefixMultiGroup ), args.legendreMax, date = dateTimeStr )
    for transportable in transportables :
        multiGroupStyle.transportables.add( transportable )
    reactionSuite.styles.add( multiGroupStyle )

### heat the cross sections and AEPs with a style for each temperature
for temperatureIndex, temperatureValue in enumerate( args.temperatures ) :

    suffix = '_%03d' % temperatureIndex

    temperature = temperatureModule.temperature( temperatureValue, args.temperatureUnit )
    if( reactionSuite.projectile != IDsPoPsModule.neutron ) :
        if( float( temperature ) > 0 ) : raise ValueError( 'Can only heat neutron as projectile cross sections.' )
        temperature = reactionSuite.styles[0].temperature
    if( float( temperature ) >= 0 ) :
        if( args.verbose > 0 ) : print( 'Heating to %s %s' % ( temperature.value, temperature.unit ) )

        heatStyle = stylesModule.heated( args.prefixHeated + suffix, preLoopStyle.label, temperature, date = dateTimeStr )
        reactionSuite.styles.add( heatStyle )
        reactionSuite.heatCrossSections( heatStyle, setThresholdToZero = True, heatBelowThreshold = False )  ### FIXME need logfile arguments in this process call
        loopStyle = heatStyle
    else :
        loopStyle = preLoopStyle

    if( args.MonteCarlo ) :
        if( args.verbose > 1 ) : print('  Processing gridded cross sections')
        griddedCrossSectionStyle = stylesModule.griddedCrossSection( args.prefixMonteCarlo + suffix, loopStyle.label, date = dateTimeStr )
        reactionSuite.styles.add( griddedCrossSectionStyle )
        reactionSuite.processGriddedCrossSections( griddedCrossSectionStyle, verbosity = args.verbose - 2, indent = '    ' )

    if( args.MultiGroup ) :
        if( args.verbose > 1 ) : print('  Processing multi group')
        heatedMultiGroupStyle = stylesModule.heatedMultiGroup( args.prefixMultiGroup + suffix, loopStyle.label, multiGroupStyle.label, flux, date = dateTimeStr )
        reactionSuite.styles.add( heatedMultiGroupStyle )
        workDir = os.path.join( 'Merced.work', reactionSuite.projectile, reactionSuite.target, "T_%s_%s" % 
                ( temperatureValue, args.temperatureUnit.replace( " ", "" ).replace( os.sep, "_" ) ) )
        if( os.path.exists( workDir ) ) : shutil.rmtree( workDir )
        reactionSuite.processMultiGroup( heatedMultiGroupStyle, verbosity = args.verbose - 2, logFile = logFile, indent = '    ', workDir = workDir )
        if args.UpScatter: 
            heatUpscatterStyle = stylesModule.UpScatter( args.prefixUpScatter + suffix, heatedMultiGroupStyle.label, date = dateTimeStr ) 
            reactionSuite.styles.add( heatUpscatterStyle )
            reactionSuite.processUP( heatUpscatterStyle, verbosity = args.verbose - 2 )  ### not implemented yet?  ### FIXME need logfile arguments in this process call

processedFileName = args.outputFile
if( processedFileName is None ) : processedFileName = args.gnds.replace( '.xml', '%s.xml' % args.tag )
if( outputDefault != args.output ) : processedFileName = os.path.join( args.output, os.path.basename( processedFileName ) )
reactionSuite.saveToFile( processedFileName, xs_pdf_cdf1d_singleLine = True )

logFile.close()
