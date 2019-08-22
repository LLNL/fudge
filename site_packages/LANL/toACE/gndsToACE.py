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

"""
Issues:
    -) When delayed neutrons are present, their multiplicities are added to prompt, and prompt and total nu_bar
        are written. However, their decay constants, etc are currently not written.
    -) Only neutrons distribution data are outputted.
"""

import time

from fudge.core.utilities import brb

from pqu import PQU as PQUModule

from xData import standards as standardsModule
from xData import axes as axesModule
from xData import values as valuesModule
from xData import Ys1d as Ys1dModule
from xData import XYs as XYsModule

from PoPs import IDs as IDsPoPsModule
from PoPs.groups import misc as miscPoPsModule

from fudge.gnds import tokens as tokensModule
from fudge.gnds import styles as stylesModule
from fudge.gnds.productData import multiplicity as multiplicityModule
from fudge.gnds.productData.distributions import angular as angularModule
from fudge.gnds.productData.distributions import energy as energyModule
from fudge.legacy.converting import endf_endl as endf_endlModule

from . import angularEnergy as angularEnergyModule
from . import specialMF6 as specialMF6Module

floatFormat = '%20.12E'
floatFormat = ' %19.11E'

def toACE( self, styleLabel, fileName, evaluationId, productData, addAnnotation ) :

    style = self.styles[styleLabel]
    if( not( isinstance( style, stylesModule.griddedCrossSection ) ) ) :
        raise TypeError( 'style labelled "%s" not a "griddedCrossSection" style' % styleLabel )

    massUnit = 'eV/c**2'
    neutron = self.PoPs[IDsPoPsModule.neutron]
    target = self.PoPs[self.target]
    if( target.id in self.PoPs.aliases ) : target = self.PoPs[target.pid]
    targetZ, targetA, targetZA, level = miscPoPsModule.ZAInfo( target )

    neutronMass = neutron.getMass( massUnit )
    targetMass = target.getMass( massUnit )
    target2NeutronMass = targetMass / neutronMass

    temperature = style.temperature
    temperature_K = temperature.getValueAs( 'K' )
    temperature_MeV = temperature.getValueAs( 'MeV/k' )

    processingTime = time.localtime( )
    processingTime = [ 2014, 5, 10 ]
    strRecords = [ "%6d.%.2dc%12.7f %11.5E %.2d/%.2d/%.2d" % ( targetZA, evaluationId, target2NeutronMass, 
            temperature_MeV, processingTime[1], processingTime[2], processingTime[0] ) ]

    nonFissionEnergyDependentNeutronMultiplicities = []
    nonFissionEnergyDependentNeutronMultiplicitiesIndex = 100

    MAT = endf_endlModule.ZAAndMATFromParticleName( target.id )[1]
    MAT_ID = '  mat %d' % MAT
    evaluated = self.styles.getEvaluatedStyle( )
    HK = '%2d-%s at %.1fK from %s-%s Fudge.toACE' % \
            ( targetZ, target.id, temperature_K, evaluated.library, evaluated.version )
    HK = "%-70s" % HK[:70]
    strRecords.append( HK + MAT_ID )
    for i in xrange( 4 ) : strRecords.append( '      0   0.000000      0   0.000000      0   0.000000      0   0.000000' )

    NXS = 16 * [ 0 ]
    JXS = 32 * [ 0 ]
    MTR = []
    NU_prompt = None
    NU_total = None
    LQR = []
    TYP = []
    SigData = {}
    neutronAngular = []
    neutronEnergies = []

    NXS[2-1] = targetZA

    energyGrid = style.grid.values

    totalXSec = Ys1dModule.Ys1d( valuesModule.values( [] ) )
    absorptionXSec = Ys1dModule.Ys1d( valuesModule.values( [] ) )

    sortedMTs = sorted( [ [ MT_MTData[0], i1 ] for i1, MT_MTData in enumerate( productData ) ] ) # Sort MTs like NJOY.
    delayedNeutronDatas = []

    for MT, i1 in sortedMTs :
        _MT, MTData = productData[i1]
        XSec = MTData['ESZ']
        neutronDatas = MTData[IDsPoPsModule.neutron]
        multiplicity = 0
        totalXSec = totalXSec + XSec
        if( MT == 2 ) :
            elasticXSec = XSec
            if( len( neutronDatas ) > 1 ) : raise Exception( 'Only one type of neutron data are supported for the elastic reaction.' )
            if( len( neutronDatas ) > 0 ) : neutronAngular.append( ( MT, neutronDatas[0]['angularData'] ) )
        else :
            MTR.append( MT )
            LQR.append( MTData['Q_initial'] )
            SigData[MT] = XSec
            if( len( neutronDatas ) == 0 ) :
                absorptionXSec = absorptionXSec + XSec
                neutronMultiplicity = 0
            else :
                NXS[5-1] += 1
                product = neutronDatas[0]['product']
                frame = neutronDatas[0]['frame']
                angularData = neutronDatas[0]['angularData']
                energyData = neutronDatas[0]['energyData']
                if( MTData['isFission'] ) :
                    neutronMultiplicity = 19
                    totalPresent = False
                    _neutronDatas = []
                    for neutronData in neutronDatas :
                        product = neutronData['product']
                        if( ( product.attributes['emissionMode'] == tokensModule.promptToken ) or ( product.attributes['emissionMode'] == 'total' ) ) :
                            totalPresent = product.attributes['emissionMode'] == 'total'
                            NU_prompt = neutronData['multiplicity']
                            _neutronDatas = []
                        else :
                            delayedNeutronDatas.append( neutronData )
                        neutronDatas = _neutronDatas
                    if( NU_prompt is None ) : raise Exception( 'Missing prompt neutron' )
                    if( len( delayedNeutronDatas ) > 0 ) :
                        if( totalPresent ) : raise Exception( 'Total nu_bar present and delayed neutron data.' )
                        NU_total = NU_prompt
                        for neutronData in delayedNeutronDatas :
                            product = neutronData['product']
                            if( product.attributes['emissionMode'] == 'delayed' ) :
                                NU_delayed = neutronData['multiplicity']
                                if( NU_delayed.domainMax < NU_total.domainMax ) : NU_delayed = NU_delayed.dullEdges( upperEps = 1e-8 )
                                NU_total = NU_total + NU_delayed
                else :
                    if( len( neutronDatas ) == 1 ) :
                        neutronMultiplicity = neutronDatas[0]['multiplicity']
                    else :
                        neutronMultiplicity, angularData, energyData = specialMF6Module.neutrons( neutronDatas )
                    if( isinstance( neutronMultiplicity, float ) ) :
                        if( neutronMultiplicity != int( neutronMultiplicity ) ) :
                            raise Exception( 'Bad neutronMultiplicity = %e' % neutronMultiplicity )
                        neutronMultiplicity = int( neutronMultiplicity )
                    if( not( isinstance( neutronMultiplicity, int ) ) ) :
                        nonFissionEnergyDependentNeutronMultiplicitiesIndex += 1
                        nonFissionEnergyDependentNeutronMultiplicities.append( [ MT, nonFissionEnergyDependentNeutronMultiplicitiesIndex, neutronMultiplicity ] )
                        neutronMultiplicity = nonFissionEnergyDependentNeutronMultiplicitiesIndex
 
                    if( frame == standardsModule.frames.centerOfMassToken ) : neutronMultiplicity *= -1  # Negative for COM frame.
                    if( energyData is None ) :
                        if( 50 <= MT <= 91 ) :
                            energyData = n_nPrimeEnergyData( MT, target2NeutronMass, MTData['Q_initial'] )
                        else :
                            raise 'hell: but first implement energy for MT =' % MT
                    else :
                        if( isinstance( energyData, energyModule.NBodyPhaseSpace ) ) : angularData = None
                neutronAngular.append( ( MT, angularData ) )
                neutronEnergies.append( [ MT, energyGrid[XSec.Ys.start], energyGrid[-1], energyData ] )
            TYP.append( neutronMultiplicity )

    annotates = []
    XSS = []
    
# 1) Add the ESZ block.
    NXS[3-1] = len( energyGrid )
    updateXSSInfo( 'energyGrid', annotates, XSS, energyGrid.values )
    updateXSSInfo( 'totalXSec', annotates, XSS, totalXSec )
    if( len( absorptionXSec ) == 0 ) :
        updateXSSInfo( 'absorption cross section', annotates, XSS, len( energyGrid ) * [ 0. ] )
    else :
        updateXSSInfo( 'absorption cross section', annotates, XSS, mapCrossSectionToGrid( absorptionXSec ) )
    updateXSSInfo( 'elastic cross section', annotates, XSS, elasticXSec )
    averageHeating = len( energyGrid ) * [ 0. ]
    updateXSSInfo( 'average heating', annotates, XSS, averageHeating )

# 2) Add the NU block.
    if( NU_prompt is not None ) :
        JXS[2-1] = len( XSS ) + 1
        NU_prompt = NU_prompt.toACE( )
        if( NU_total is not None ) : NU_prompt.insert( 0, -len( NU_prompt ) )
        updateXSSInfo( 'NU', annotates, XSS, NU_prompt )
        if( NU_total is not None ) : 
            NU_total = NU_total.toACE( )
            updateXSSInfo( 'NU(total)', annotates, XSS, NU_total )

# 3) Add the MTR block.
    NXS[4-1] = len( MTR )
    JXS[3-1] = len( XSS ) + 1
    updateXSSInfo( 'MTR', annotates, XSS, MTR )

# 4) Add the LQR block.
    JXS[4-1] = len( XSS ) + 1
    updateXSSInfo( 'LQR', annotates, XSS, LQR )

# 5) Add the TYP block.
    JXS[5-1] = len( XSS ) + 1
    updateXSSInfo( 'TYP', annotates, XSS, TYP )

# 6 and 7) Add the LSIG and SIG blocks.
    SIG = []
    for MT, i1 in sortedMTs :
        _MT, MTData = productData[i1]
        if( MT not in [ 2 ] ) : addSigData( MT, SIG, SigData[MT] )
    JXS[6-1] = len( XSS ) + 1
    LSIG = [ 1 ]
    for MT, firstNonZero, reactionSIG in SIG[:-1] : LSIG.append( LSIG[-1] + len( reactionSIG ) + 2 )
    updateXSSInfo( 'LSIG', annotates, XSS, LSIG )
    JXS[7-1] = len( XSS ) + 1
    for MT, firstNonZero, reactionSIG in SIG :
        reactionSIG_ = [ firstNonZero + 1, len( reactionSIG ) ] + reactionSIG
        updateXSSInfo( 'SIG(MT=%s)' % MT, annotates, XSS, reactionSIG_ )

# 8 and 9) Add the LAND and AND blocks.
    LAND, MT_AND, length = [], [], 0
    for MT, angular in neutronAngular :
        if( angular is None ) :
            LAND.append( -1 )
        elif( angular.isIsotropic( ) ) :
            LAND.append( 0 )
        elif( isinstance( angular, ( angularModule.XYs2d, angularModule.regions2d ) ) ) :
            if( isinstance( angular, angularModule.regions2d ) ) :
                regions = angular
                angular = []
                for region in regions :
                    for xs_pdf_cdf in region : angular.append( xs_pdf_cdf )
            LAND.append( length + 1 )
            length += 2 * len( angular ) + 1
            energies_in = [ len( angular ) ]
            LCs = []
            Ps = []
            for xs_pdf_cdf in angular :
                energies_in.append( xs_pdf_cdf.value )
                LCs.append( -length - 1 )
                Ps += [ 2, len( xs_pdf_cdf ) ] + xs_pdf_cdf.xs.values.values + xs_pdf_cdf.pdf.values.values + xs_pdf_cdf.cdf.values.values
                length += 3 * len( xs_pdf_cdf ) + 2
            MT_AND.append( ( MT, energies_in + LCs + Ps ) )
        else :
            raise Exception( 'Unsupport neutron angular distribution type = %s' % type( angular ) )
    JXS[8-1] = len( XSS ) + 1
    updateXSSInfo( 'LAND', annotates, XSS, LAND )
    JXS[9-1] = len( XSS ) + 1
    for MT, AND in MT_AND : updateXSSInfo( 'AND(MT=%s)' % MT, annotates, XSS, AND )

# 10 and 11) Add the LDLW and DLW blocks.
    LDLW, MT_DLW = processEnergyData( massUnit, neutronMass, neutronEnergies )
    JXS[10-1] = len( XSS ) + 1
    updateXSSInfo( 'LDLW', annotates, XSS, LDLW )
    JXS[11-1] = len( XSS ) + 1
    for MT, DLW in MT_DLW : updateXSSInfo( 'DLW(MT=%s)' % MT, annotates, XSS, DLW )

# 24, 25, 26 and 27) Added the DNU, BDD, DNEDL and DNED blocks

    if( len( delayedNeutronDatas ) > 0 ) :
        skipDelayed = False
        for delayedNeutronData in delayedNeutronDatas :
            if( delayedNeutronData['energyData'] is None ) : skipDelayed = True
        if( skipDelayed ) : delayedNeutronDatas = []

    if( len( delayedNeutronDatas ) > 0 ) :
        NXS[8-1] = len( delayedNeutronDatas )
        delayedFissionNeutrons = []
        axes = axesModule.axes( labelsUnits = { 1 : ( '', 'MeV' ), 0 : ( '', '' ) } )
        totolDelayedMultiplicity = XYsModule.XYs1d( data = [], axes = axes )
        for delayedNeutronData in delayedNeutronDatas :
            decayRate = 1e-8 * delayedNeutronData['product'].attributes['decayRate'].getValueAs( '1/s' )
            multiplicity = delayedNeutronData['multiplicity']
            if( not( isinstance( multiplicity, multiplicityModule.XYs1d ) ) ) : raise Exception( 'fix me' )
            if( multiplicity.interpolation != standardsModule.interpolation.linlinToken ) : raise Exception( 'fix me' )
            delayedFissionNeutrons.append( [ decayRate, multiplicity ] )
            totolDelayedMultiplicity += multiplicity
            
        energies, multiplicities = [], []
        for energy, multiplicity in totolDelayedMultiplicity :
            energies.append( energy )
            multiplicities.append( multiplicity )
        DNU = [ 2, 0, len( totolDelayedMultiplicity ) ] + energies + multiplicities
        JXS[24-1] = len( XSS ) + 1
        updateXSSInfo( 'DNU', annotates, XSS, DNU )

        BDD = []
        for decayRate, multiplicity in delayedFissionNeutrons :
            probability = multiplicity / totolDelayedMultiplicity
            probability = probability.thin( 1e-3 )
            energies, probabilities = [], []
            for energy, _probability in probability :
                energies.append( energy )
                probabilities.append( _probability )

            BDD += [ decayRate, 0, len( probability ) ] + energies + probabilities
        JXS[25-1] = len( XSS ) + 1
        updateXSSInfo( 'BDD', annotates, XSS, BDD )

        delayedNeutronEnergies = []
        for delayedNeutronData in delayedNeutronDatas :
            multiplicity = delayedNeutronData['multiplicity']
            delayedNeutronEnergies.append( [ 18, multiplicity.domainMin, multiplicity.domainMax, delayedNeutronData['energyData'] ] )
        LDLW, MT_DLW = processEnergyData( massUnit, neutronMass, delayedNeutronEnergies )
        JXS[26-1] = len( XSS ) + 1
        updateXSSInfo( 'DNEDL', annotates, XSS, LDLW )
        JXS[27-1] = len( XSS ) + 1
        for MT, DLW in MT_DLW : updateXSSInfo( 'DNED(MT=%s)' % MT, annotates, XSS, DLW )

# 12) GPD block. Not needed.

# 13) MTRP block (gammas).

# 14) LSIGP block (gammas).

# 15) SIGP block (gammas).

# 16) LANDP block (gammas).

# 17) ANDP block (gammas).

# 18) LDLWP block (gammas).

# 19) DLWP block (gammas).

# 20) YP block (gammas).

# 21) FIS block. Not needed.

# 22) UNR block.

# Now fixup the TYP data whose abs( value ) is greater than 100.
    for MT, index, multiplicity in nonFissionEnergyDependentNeutronMultiplicities :
        i1 = JXS[5-1]
        found = False
        for i2 in range( NXS[4-1] ) :
            if( abs( XSS[i1 + i2 - 1] ) ==  index ) :
                n_multiplicity = multiplicity.toACE( )
                LNU, n_multiplicity = n_multiplicity[0], n_multiplicity[1:]
                if( LNU != 2 ) : raise Exception( 'Only tabular neutron multiplicity is allowed, not LNU = %d' % LNU )
                offset = len( XSS ) - JXS[11-1] + 101 + 1
                if( XSS[i1+i2-1] < 0 ) : offset *= -1
                XSS[i1+i2-1] = offset
                updateXSSInfo( 'Neutron Yields(MT=%s)' % MT, annotates, XSS, n_multiplicity )
                found = True
                break
        if( not( found ) ) : raise Exception( 'Neutron multiplicity for index %s not found in TYP table' % index )

# Time to wrap it all up.
    NXS[1-1] = len( XSS )
    JXS[1-1] = 1
    JXS[22-1] = len( XSS ) + 1      # Location of the last word of XSS (i.e., its length).

    strRecords += intArrayToRecords( NXS )
    strRecords += intArrayToRecords( JXS )
    strRecords += XSSToStrings( annotates, XSS, addAnnotation )

    strRecords.append( '' )
    fOut = open( fileName, 'w' )
    fOut.write( '\n'.join( strRecords ) )
    fOut.close( )        

def updateXSSInfo( label, annotates, XSS, data ) :

    annotates.append( ( len( XSS ), label ) )
    XSS += data

def intArrayToRecords( _array ) :

    s = [ "%9d" % d for d in _array ]
    records = []
    while( len( s ) ) :
        record = []
        for i1 in xrange( 8 ) : record.append( s.pop( 0 ) )
        records.append( ''.join( record ) )
    return( records )

def XSSToStrings( annotates, XSS, addAnnotation ) :

    strData, record, i3 = [], [], 0
    annotates.append( ( len( XSS ), '' ) )
    i2, label = annotates[i3]
    for i1, datum in enumerate( XSS ) :
        if( type( datum ) == int ) :
            record.append( "%20d" % datum )
        else :
            try :
                record.append( floatFormat % datum )
            except :
                print i1, len( XSS ), datum, annotates[i3]
                raise
        if( len( record ) == 4 ) :
            annotate = ''
            if( i2 <= i1 ) :
                annotate = ' !'
                sep = ' '
                while( i2 <= i1 ) :
                    annotate += sep + label + ' (%s)' % ( i2 - i1 + 3 )
                    sep = ', '
                    i3 += 1
                    i2, label = annotates[i3]
            if( not addAnnotation ) :   annotate = ''
            strData.append( ''.join( record ) + annotate )
            record = []
    if( record ) : strData.append( ''.join( record ) )
    return( strData )

def mapCrossSectionToGrid( XSec ) :

    Ys = XSec.Ys
    return( Ys.start * [ 0.0 ] + Ys.values )

def addSigData( MT, SIG, SigData ) :

    SIG.append( ( MT, SigData.Ys.start, SigData.Ys.values ) )

class n_nPrimeEnergyData :

    def __init__( self, MT, A, Q ) :

        self.LAW = 3
        self.MT = MT
        self.A = float( A )
        self.Q = float( Q )

    def toACE( self, label, offset, weight, **kwargs ) :

        r1 = self.A / ( 1. + self.A )
        return( [ 0, self.LAW, offset + len( weight ) + 4 ] + weight + [ abs( self.Q ) / r1, r1 * r1 ] )

def processEnergyData( massUnit, neutronMass, energyDatas ) :

    LDLW, MT_DLW, length = [], [], 0
    for MT, EMin, EMax, energyData in energyDatas :
        LDLW.append( length + 1 )
        defaultWeight = [ 0, 2, EMin, EMax, 1.0, 1.0 ]
        label = 'MT = %s for product = "neutron"' % MT
        _DLW = energyData.toACE( label, length, defaultWeight, massUnit = massUnit, neutronMass = neutronMass )
        MT_DLW.append( ( MT, _DLW ) )
        length += len( _DLW )
    return( LDLW, MT_DLW )
