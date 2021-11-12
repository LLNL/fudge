#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import copy
from argparse import ArgumentParser

from xData import XYs as XYsModule

from brownies.legacy.endl import fudgeParameters

fudgeParameters.VerboseMode = 0

from brownies.legacy.endl import endlProject as endlProjectClass
from fudge import reactionSuite as reactionSuiteModule

from fudge.productData.distributions import \
        angularEnergy as angularEnergyModule, \
        LLNL_angularEnergy as LLNL_angularEnergyModule

from fudge.reactionData.doubleDifferentialCrossSection.thermalNeutronScatteringLaw import \
        coherentElastic as coherentElasticModule, \
        incoherentElastic as incoherentElasticModule, \
        incoherentInelastic as incoherentInelasticModule

description1 = """
This module converts processed data from a GNDS theramal scattering law file into ENDL and adds it to a base endl ZA directory.
Thermal scattering law parameters should already be expanded into double-differential form, e.g. using processProtare with -mc option.
Sample usage:

    processProtare.py /path/to/tsl-HinCH2.gnds.xml tsl-HinCH2.proc.xml -mc -t 2.586e-8
    TSL2ENDL.py tsl-HinCH2.proc.xml heated_000 /usr/gapps/data/nuclear/endl_official/endl2009.2 1001 1901
"""

# Jurgenson recommended names:
# Name                  Endl    ENDF 7.1       ENDF 8.0 b4   endlName
#
# H in H2O               x          y               z          1801
# H in CH2               x          x               y          1901
# D in D2O               x          y               z          1902
# Be-metal               x          x               y          4809
# Be in BeO              x          y               z          4909
# Graphite               x          x               y          6912
# O in BeO               x          y               z          8916
# Zr in ZrH                         x               x          40900
# H in ZrH                          x               x          1701
# Al                                x               x          13927
# Fe                                x               x          26956
# U in UO                           x               y          92938
# O in UO                           x               y          8816
#
# Benzene                           x               x
# Ortho-D                           x               x
# Ortho-H                           x               x
# Para-D                            x               x
# Para-H                            x               x
# s-methane                         x               x
# l-methane                         x               x
#
# SiO2 alpha                        x (as just SiO) y          14928
# SiO2 beta                                         x          14828
# C in SiC                                          x          6712
# H in C5O2H8  (Lucite)                             x          1401
# H in Ice                                          x          1601
# H in YH2                                          x          1501
# N in UN                                           x          7914
# O in D2O                                          x          8716
# O in Ice                                          x          8616
# Si in SiC                                         x          14728
# U in UN                                           x          92838
# Y in YH2                                          x          39989
# Reactor graphite - 10P                            x          6812
# Reactor graphite - 30P                            x          6612
#
# Letters to indicate successive versions of each evaluation in ENDF/B-VIII-beta4

__doc__ = description1

parser = ArgumentParser( description = description1 )
parser.add_argument( 'TNSLFile', type = str,                                            help = 'Path to the GNDS thermal scattering law file.' )
parser.add_argument( 'style', type = str,                                               help = 'Style of data to extract from GNDS data file.' )
parser.add_argument( 'endlDatabase', type = str,                                        help = 'Path to the endl database containing the ZA data to use as the base endl ZA (e.g., "/usr/gapps/data/nuclear/endl_official/endl2009.2").' )
parser.add_argument( 'ZA', type = int,                                                  help = 'ZA in the endl database to update with the thermal scattering law data (e.g., 4009).' )
parser.add_argument( 'tslZA', type = int,                                               help = 'ZA to use for thermal scattering law data (e.g., 4909).' )
parser.add_argument( '-o', '--output', type = str, action = 'store', default = './',    help = 'The path where the endlZA with the thermal scattering law data will be written. Default is "./".' )
parser.add_argument( '--onlyTNSL', action = 'store_true',                               help = 'If present, only the TNSL data are written, all other data are removed.' )

args = parser.parse_args( )

TNSLFile = reactionSuiteModule.readXML( args.TNSLFile )
#TNSLFile.convertUnits( {'eV':'MeV'} )   # not working right now, onus on user to ensure data are processed in MeV

endlProject = endlProjectClass( database = args.endlDatabase, projectile = 1, readOnly = True )
endlZA = endlProject.readZA( args.ZA )
endlZA.read( )

I0CrossSection = endlZA.findData( C = 10, I = 0 )
temperature_MeV = I0CrossSection.getTemperature( )

def getEndlFile( endlZA, yo, C, I, S ) :

    endlFiles = endlZA.findFiles( yo, C, I, S )
    if( len( endlFiles ) == 0 ) : return( endlZA.addFile( yo, C, I, S ) )
    return( endlFiles[0] )

def addData( endlZA, X1, temperature, energyMin_MeV, energyMax_MeV, I0Data, I1Data, I3Data, halflife = None ) :

    endlFile = getEndlFile( endlZA, 0, 11, 0, 1 )
    data = [ [ energy, xSec ] for energy, xSec in I0Data ]
    if( data[0][0] > energyMin_MeV ) :
        if( data[0][1] != 0 ) : data.insert( 0, [ data[0][0], 0 ] )
        data.insert( 0, [ energyMin_MeV, 0 ] )
    if not args.onlyTNSL:
        if( data[-1][0] < energyMax_MeV ) :
            if( data[-1][1] != 0 ) : data.append( [ data[-1][0], 0 ] )
            data.append( [ energyMax_MeV, 0 ] )
    endlData = endlFile.addData( data, Q = 0, X1 = X1, temperature = temperature, halflife = halflife )
    endlData.setFormat( 15 )

    if( I1Data[0][0] > energyMin_MeV ) : I1Data.insert( 0, [ energyMin_MeV, copy.copy( I1Data[0][1] ) ] )
    if not args.onlyTNSL:
        if( I1Data[-1][0] < energyMax_MeV ) : I1Data.append( [ energyMax_MeV, copy.copy( I1Data[-1][1] ) ] )
    _I1Data = []
    for energy, PofMu in I1Data :
        PofMu = XYsModule.XYs1d( PofMu ).normalize( )
        _I1Data.append( [ energy, [ [ P, Mu ] for P, Mu in PofMu ] ] )

    endlFile = getEndlFile( endlZA, 1, 11, 1, 1 )
    endlData = endlFile.addData( _I1Data, Q = 0, X1 = X1, temperature = temperature, halflife = halflife )
    endlData.setFormat( 15 )

    if I3Data is not None:
        if( I3Data[0][0] > energyMin_MeV ) : I3Data.insert( 0, [ energyMin_MeV, copy.copy( I3Data[0][1] ) ] )
        if not args.onlyTNSL:
            if( I3Data[-1][0] < energyMax_MeV ) : I3Data.append( [ energyMax_MeV, copy.copy( I3Data[-1][1] ) ] )

        _I3Data = []
        for energy, PofEprimeMu in I3Data:
            tmp = []
            for mu, PofEprime in PofEprimeMu:
                PofEprime = XYsModule.XYs1d( PofEprime ).normalize( )
                tmp.append( [ mu, [ [ Ep, P ] for Ep, P in PofEprime ] ] )
            _I3Data.append( [ energy, tmp ] )

        endlFile = getEndlFile( endlZA, 1, 11, 3, 1 )
        endlData = endlFile.addData( _I3Data, Q = 0, X1 = X1, temperature = temperature, halflife = halflife )
        endlData.setFormat( 15 )

energyMin_MeV = min( 1e-11, I0CrossSection.data[0][0] )
energyMax_MeV = I0CrossSection.data[-1][0]
halflife = I0CrossSection.getHalflife()

def findMaxEnergy( xsc ):
    """ find energy of maximum point with non-zero cross section """
    for x,y in xsc.copyDataToXYs()[::-1]:
        if y > 0:
            break
    return x

TNSLEnergyMax_MeV = []
coherentElastic = incoherentElastic = incoherentInelastic = None
for reaction in TNSLFile.reactions:
    if not reaction.doubleDifferentialCrossSection:
        continue
    evaluated = reaction.doubleDifferentialCrossSection.evaluated
    TNSLEnergyMax_MeV.append( findMaxEnergy( reaction.crossSection[args.style] ) )
    if isinstance(evaluated, coherentElasticModule.form):
        coherentElastic = reaction
    elif isinstance(evaluated, incoherentElasticModule.form):
        incoherentElastic = reaction
    elif isinstance(evaluated, incoherentInelasticModule.form):
        incoherentInelastic = reaction
    else:
        print("WARNING: encountered unknown double-differential cross section form '%s'" % evaluated)

if len(set(TNSLEnergyMax_MeV)) > 1:
    print("WARNING: inconsistent upper energy limits for TNSL cross sections!")
    print(TNSLEnergyMax_MeV)
    print("Using minimum cutoff for transition")
    TNSLEnergyMax_MeV = min(TNSLEnergyMax_MeV)
else:
    TNSLEnergyMax_MeV = TNSLEnergyMax_MeV[0]

if incoherentInelastic is None:
    raise Exception( 'Missing incoherent inelastic thermal scattering law data.' )

if args.onlyTNSL: endlZA.removeFile( )

if coherentElastic is not None:
    I0Data = coherentElastic.crossSection[args.style]
    I1Data = [[xys1d.outerDomainValue, xys1d.copyDataToXYs()] for xys1d in
            coherentElastic.outputChannel.getProductWithName("n").distribution[args.style].angularSubform.data]
    I3Data = None
    addData( endlZA, 2e-12, temperature_MeV, energyMin_MeV, energyMax_MeV, I0Data, I1Data, I3Data, halflife=halflife )

if incoherentElastic is not None:
    I0Data = incoherentElastic.crossSection[args.style]
    I1Data = [[xys1d.outerDomainValue, xys1d.copyDataToXYs()] for xys1d in
            incoherentElastic.outputChannel.getProductWithName("n").distribution[args.style].angularSubform.data]
    I3Data = None
    addData( endlZA, 2e-12, temperature_MeV, energyMin_MeV, energyMax_MeV, I0Data, I1Data, I3Data, halflife=halflife )

if incoherentInelastic is not None:
    I0Data = incoherentInelastic.crossSection[args.style]
    I0Data = I0Data.domainSlice(energyMin_MeV, TNSLEnergyMax_MeV)

    I1Data = []
    I3Data = []
    distribution = incoherentInelastic.outputChannel.getProductWithName("n").distribution[args.style]
    if isinstance(distribution, LLNL_angularEnergyModule.LLNLAngularEnergyForm):
        # older versions of processProtare used LLNL-specific distribution form
        angular = distribution.angularSubform.data.domainSlice(energyMin_MeV, TNSLEnergyMax_MeV)
        angularEnergy = distribution.angularEnergySubform.data.domainSlice(energyMin_MeV, TNSLEnergyMax_MeV)
        for xys1d in angular:
            I1Data.append( [xys1d.outerDomainValue, xys1d.copyDataToXYs()] )
        for xys2d in angularEnergy:
            I3Data.append( [xys2d.outerDomainValue, [[xys1d.outerDomainValue, xys1d.copyDataToXYs()] for xys1d in xys2d]] )
    elif isinstance(distribution, angularEnergyModule.form):
        # produced by newer version of processProtare. Need to renormalize
        angularEnergy = distribution.angularEnergySubform
        for xys2d in angularEnergy:
            I1Data.append( [xys2d.outerDomainValue,
                [[xys1d.outerDomainValue, float(xys1d.integrate())] for xys1d in xys2d]] )
            I3Data.append( [xys2d.outerDomainValue,
                [[xys1d.outerDomainValue, xys1d.normalize().copyDataToXYs()] for xys1d in xys2d]] )

    addData( endlZA, 1e-12, temperature_MeV, energyMin_MeV, energyMax_MeV, I0Data, I1Data, I3Data, halflife=halflife )

for data in endlZA.findDatas( ) : data.setZA( args.tslZA )

if( not args.onlyTNSL ) :
    I0CrossSection.setValue( TNSLEnergyMax_MeV, I0CrossSection.getValue( TNSLEnergyMax_MeV ) )
    for index, energyCrossSection in enumerate( I0CrossSection.data ) :
        if( energyCrossSection[0] == TNSLEnergyMax_MeV ) : break
    I0CrossSection.data = I0CrossSection.data[index:]
    I0CrossSection.data.insert( 0, [ TNSLEnergyMax_MeV, 0 ] )
    I0CrossSection.data.insert( 0, [ energyMin_MeV, 0 ] )

output = os.path.join( args.output, "za%.6d" % args.tslZA )
if( os.path.exists( output ) ) : os.system( 'rm -rf %s' % output )
endlZA.save( path = output )
