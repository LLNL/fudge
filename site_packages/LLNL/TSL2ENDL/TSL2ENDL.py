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

import os
import copy
from argparse import ArgumentParser

from xData import XYs as XYsModule

from fudge import fudgeParameters
fudgeParameters.VerboseMode = 0

from fudge.legacy.endl import endlProject as endlProjectClass
from fudge.gnd import thermalScattering as thermalScatteringModule

import coherentElastic as coherentElasticModule
import incoherentElastic as incoherentElasticModule
import incoherentInelastic as incoherentInelasticModule

description1 = """
This module converts data in a GND theramal scattering law file into ENDL data and adds it to a base endl ZA directory.

    TSL2ENDL.py /usr/workspace/wsa/ndg/Sab/tsl-HinCH2.gnd.xml /usr/gapps/data/nuclear/endl_official/endl2009.2 1001 1901
"""

__doc__ = description1

parser = ArgumentParser( description = description1 )
parser.add_argument( 'TSLFile', type = str,
        help = 'Path to the GND thermal scattering law file.' )
parser.add_argument( 'endlDatabase', type = str,
        help = 'Path to the endl database containing the ZA data to use as the base endl ZA' +
        '(e.g., "/usr/gapps/data/nuclear/endl_official/endl2009.2").' )
parser.add_argument( 'ZA', type = int,
        help = 'ZA in the endl database to update with the thermal scattering law data (e.g., 4009).' )
parser.add_argument( 'tslZA', type = int,
        help = 'ZA to use for thermal scattering law data (e.g., 4909).' )
parser.add_argument( '-o', '--output', type = str, action = 'store', default = './',
        help = 'The path where the endlZA with the thermal scattering law data will be written. Default is "./".' )

args = parser.parse_args( )

TSLFile = thermalScatteringModule.readXML( args.TSLFile )

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
    if( data[-1][0] < energyMax_MeV ) :
        if( data[-1][1] != 0 ) : data.append( [ data[-1][0], 0 ] )
        data.append( [ energyMax_MeV, 0 ] )
    endlData = endlFile.addData( data, Q = 0, X1 = X1, temperature = temperature, halflife = halflife )
    endlData.setFormat( 15 )

    if( I1Data[0][0] > energyMin_MeV ) : I1Data.insert( 0, [ energyMin_MeV, copy.copy( I1Data[0][1] ) ] )
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
if( TSLFile.incoherentInelastic is None ) : raise Exception( 'Missing incoherent inelastic thermal scattering law data.' )
TSLEnergyMax_MeV = TSLFile.incoherentInelastic.scatteringAtoms[0].e_max.copyToUnit( 'MeV' ).value

if( TSLFile.coherentElastic is not None ) :
    I0Data, I1Data, I3Data = coherentElasticModule.toENDL( TSLFile.coherentElastic, energyMin_MeV, TSLEnergyMax_MeV, temperature_MeV )
    addData( endlZA, 2e-12, temperature_MeV, energyMin_MeV, energyMax_MeV, I0Data, I1Data, I3Data, halflife=halflife )

if( TSLFile.incoherentElastic is not None ) :
    I0Data, I1Data, I3Data = incoherentElasticModule.toENDL( TSLFile.incoherentElastic, energyMin_MeV, TSLEnergyMax_MeV, temperature_MeV )
    addData( endlZA, 2e-12, temperature_MeV, energyMin_MeV, energyMax_MeV, I0Data, I1Data, I3Data, halflife=halflife )

if( TSLFile.incoherentInelastic is not None ) :
    I0Data, I1Data, I3Data = incoherentInelasticModule.toENDL( TSLFile.incoherentInelastic, energyMin_MeV, TSLEnergyMax_MeV, temperature_MeV )
    addData( endlZA, 1e-12, temperature_MeV, energyMin_MeV, energyMax_MeV, I0Data, I1Data, I3Data, halflife=halflife )

for data in endlZA.findDatas( ) : data.setZA( args.tslZA )

I0CrossSection.setValue( TSLEnergyMax_MeV, I0CrossSection.getValue( TSLEnergyMax_MeV ) )
for index, energyCrossSection in enumerate( I0CrossSection.data ) :
    if( energyCrossSection[0] == TSLEnergyMax_MeV ) : break
I0CrossSection.data = I0CrossSection.data[index:]
I0CrossSection.data.insert( 0, [ TSLEnergyMax_MeV, 0 ] )
I0CrossSection.data.insert( 0, [ energyMin_MeV, 0 ] )

output = os.path.join( args.output, "za%.6d" % args.tslZA )
if( os.path.exists( output ) ) : os.system( 'rm -rf %s' % output )
endlZA.save( path = output )
