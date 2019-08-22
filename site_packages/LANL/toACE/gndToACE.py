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

"""
Issues:
    -) When delayed neutrons are present, their multiplicities are added to prompt, and prompt and total nu_bar
        are written. However, their decay constants, etc are currently not written.
    -) Only neutrons distribution data are outputted.
"""

import time
from fudge.core.math.xData import axes, XYs
from fudge.gnd import tokens
from fudge.gnd.productData import distributions
from fudge.legacy.converting import endf_endl
from pqu import PQU

import angularEnergy, specialMF6

floatFormat = '%20.12E'
floatFormat = ' %19.11E'

def toACE( self, fileName, evaluationId, temperature, productData, addAnnotation ) :

    massUnit = 'eV/c**2'
    projectile, target = self.projectile, self.target
    targetZ, targetA, targetSuffix, targetZA = target.getZ_A_SuffixAndZA( )
    neutronMass = projectile.getMass( massUnit )
    targetMass = target.getMass( massUnit )
    target2NeutronMass = targetMass / neutronMass

    if( not( isinstance( temperature, PQU.PQU ) ) ) :
        temperature = PQU.PQU( temperature )

    if( temperature.isTemperature( ) ) :
        temperature_MeV = temperature.getValueAs( 'MeV / k' )
        temperature_K = temperature.getValueAs( 'K' )
    else :
        temperature_MeV = temperature.getValueAs( 'MeV' )
        temperature_K = temperature.getValueAs( 'k * K' )

    processingTime = time.localtime( )
    strRecords = [ "%6d.%.2dc%12.7f %11.5E %.2d/%.2d/%.2d" % ( targetZA, evaluationId, target2NeutronMass, temperature_MeV, processingTime[1],
        processingTime[2], processingTime[0] ) ]

    nonFissionEnergyDependentNeutronMultiplicities = []
    nonFissionEnergyDependentNeutronMultiplicitiesIndex = 100

    MAT = endf_endl.ZAAndMATFromParticleName( target.name )[1]
    MAT_ID = '  mat %d' % MAT
    HK = '%2d-%s at %.1fK from %s-%s Fudge.toACE' % \
        ( targetZ, target.name, temperature_K, self.styles['evaluated'].attributes['library'], self.styles['evaluated'].attributes['version'] )
    HK = "%-70s" % HK[:70]
    strRecords.append( HK + MAT_ID )
    for i in xrange( 4 ) : strRecords.append( '      0   0.000000      0   0.000000      0   0.000000      0   0.000000' )

    NXS = 16 * [ 0 ]
    JXS = 32 * [ 0 ]
    MRT, NU_prompt, NU_total, LQR, TYP, SigData, neutronAngular, neutronEnergies = [], None, None, [], [], {}, [], []
    NXS[2-1] = targetZA

    for MT, MTData in productData :
        if( MT == 2 ) :                 # Elastic (MT = 2) must always be present in ACE file.
            EMin = MTData['ESZ'].domainMin( unitTo = 'MeV' )
            break
    totalXSec = XYs.XYs( axes.defaultAxes( labelsUnits = { 0 : ( '', 'MeV' ), 1 : ( '', 'b' ) } ), [], 1e-3 )
    absorptionXSec = XYs.XYs( axes.defaultAxes( labelsUnits = { 0 : ( '', 'MeV' ), 1 : ( '', 'b' ) } ), [], 1e-3 )
    sortedMTs = sorted( [ [ MTData[0], i1 ] for i1, MTData in enumerate( productData ) ] ) # Sort MTs like NJOY.
    for MT, i1 in sortedMTs :
        MT_, MTData = productData[i1]
        XSec = MTData['ESZ']
        XSec = XSec.convertAxisToUnit( 0, 'MeV' )
        XSec = XSec.convertAxisToUnit( 1, 'b' )
        neutronDatas = MTData['n']
        multiplicity = 0
        if( MT == 2 ) :
            elasticXSec = XSec
            if( len( neutronDatas ) > 1 ) : raise Exception( 'Only one type of neutron data is support for the elastic reaction.' )
            if( len( neutronDatas ) > 0 ) : neutronAngular.append( ( MT, neutronDatas[0]['angularData'] ) )
        else :
            MRT.append( MT )
            LQR.append( MTData['Q'] )
            if( XSec[0][0] != 0 ) :
                if( XSec[0][0] > EMin ) : XSec = XSec.dullEdges( lowerEps = 1e-8 )
            SigData[MT] = XSec
            totalXSec = totalXSec + XSec
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
                    for neutronData in neutronDatas :
                        product = neutronData['product']
                        if( ( product.attributes['emissionMode'] == tokens.promptToken ) or ( product.attributes['emissionMode'] == 'total' ) ) :
                            totalPresent = product.attributes['emissionMode'] == 'total'
                            NU_prompt = neutronData['multiplicity']
                            break
                    if( NU_prompt is None ) : raise Exception( 'Missing prompt neutron' )
                    if( len( neutronDatas ) > 1 ) :
                        if( totalPresent ) : raise Exception( 'Total nu_bar present and delayed neutron data.' )
                        NU_total = NU_prompt
                        for neutronData in neutronDatas :
                            product = neutronData['product']
                            if( product.attributes['emissionMode'] == 'delayed' ) :
                                NU_delayed = neutronData['multiplicity']
                                if( NU_delayed.domainMax( ) < NU_total.domainMax( ) ) : NU_delayed = NU_delayed.dullEdges( upperEps = 1e-8 )
                                NU_total = NU_total + NU_delayed
                else :
                    if( len( neutronDatas ) == 1 ) :
                        neutronMultiplicity = neutronDatas[0]['multiplicity']
                    else :
                        neutronMultiplicity, angularData, energyData = specialMF6.neutrons( neutronDatas )
                    if( not ( isinstance( neutronMultiplicity, int ) ) ) :
                        nonFissionEnergyDependentNeutronMultiplicitiesIndex += 1
                        nonFissionEnergyDependentNeutronMultiplicities.append( [ MT, nonFissionEnergyDependentNeutronMultiplicitiesIndex, neutronMultiplicity ] )
                        neutronMultiplicity = nonFissionEnergyDependentNeutronMultiplicitiesIndex
 
                    if( frame == axes.centerOfMassToken ) : neutronMultiplicity *= -1  # Negative for COM frame.
                    if( energyData is None ) :
                        if( 50 <= MT <= 91 ) :
                            energyData = n_nPrimeEnergyData( MT, target2NeutronMass, MTData['Q'] )
                        else :
                            raise 'hell: but first implement energy for MT =' % MT
                    else :
                        if( isinstance( energyData, distributions.energy.NBodyPhaseSpace ) ) : angularData = None
                neutronAngular.append( ( MT, angularData ) )
                neutronEnergies.append( [ MT, XSec.xMin( ), XSec.xMax( ), energyData ] )
            TYP.append( neutronMultiplicity )
    totalXSec = elasticXSec + totalXSec
    annotates, XSS, energyGrid, totalSigma = [], [], [], []
    for E, sigma in totalXSec :
        energyGrid.append( E )
        totalSigma.append( sigma )

# 1) Add the ESZ block.
    NXS[3-1] = len( energyGrid )
    updateXSSInfo( 'energyGrid', annotates, XSS, energyGrid )
    updateXSSInfo( 'totalSigma', annotates, XSS, totalSigma )
    if( len( absorptionXSec ) == 0 ) :
        updateXSSInfo( 'absorption cross section', annotates, XSS, len( energyGrid ) * [ 0. ] )
    else :
        updateXSSInfo( 'absorption cross section', annotates, XSS, mapEnergyToTotal( totalXSec, absorptionXSec ))
    updateXSSInfo( 'elastic cross section', annotates, XSS, mapEnergyToTotal( totalXSec, elasticXSec ) )
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

# 3) Add the MRT block.
    NXS[4-1] = len( MRT )
    JXS[3-1] = len( XSS ) + 1
    updateXSSInfo( 'MRT', annotates, XSS, MRT )

# 4) Add the LQR block.
    JXS[4-1] = len( XSS ) + 1
    updateXSSInfo( 'LQR', annotates, XSS, LQR )

# 5) Add the TYP block.
    JXS[5-1] = len( XSS ) + 1
    updateXSSInfo( 'TYP', annotates, XSS, TYP )

# 6 and 7) Add the LSIG and SIG blocks.
    SIG = []
    for MT, MTData in productData :
        if( MT not in [ 2 ] ) : addSigData( MT, SIG, energyGrid, SigData[MT] )
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
        elif( isinstance( angular, ( distributions.angular.linear, angularEnergy.angularFor_angularEnergy ) ) ) :
            LAND.append( length + 1 )
            length += 2 * len( angular ) + 1
            energies_in = [ len( angular ) ]
            LCs, Ps = [], []
            for xys in angular :
                energy_in = PQU.PQU( xys.value, angular.axes[0].getUnit( ) ).getValueAs( 'MeV' )
                energies_in.append( energy_in )
                LCs.append( -length - 1 )
                xys = xys.normalize( )
                mus, pdf = [], []
                for x, y in xys :
                    mus.append( x )
                    pdf.append( y )
                cdf = xys.runningIntegral( )
                cdf[-1] = 1.                    # Make sure last point is really 1.
                Ps += [ 2, len( mus ) ] + mus + pdf + cdf
                length += 3 * len( xys ) + 2
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

# Now fixup the TYP data whose abs( value ) is greater than 100.
    for MT, index, multiplicity in nonFissionEnergyDependentNeutronMultiplicities :
        i1, found = JXS[5-1], False
        for i2 in range( NXS[4-1] ) :
            if( abs( XSS[i1 + i2 - 1] ) ==  index ) :
                n_multiplicity = multiplicity.toACE( )
                LNU, n_multiplicity = n_multiplicity[0], n_multiplicity[1:]
                if( LNU != 2 ) : raise Exception( 'Only tabular neutron multiplicity is allowed, not LNU = %d' % LNU )
                offset = len( XSS ) - JXS[11-1] + 101 + 1
                if( XSS[i1 + i2 - 1] < 0 ) : offset *= -1
                XSS[i1 + i2 - 1] = offset
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

def mapEnergyToTotal( totalXSec, XSec ) :

    return( [ XSec.getValue( E ) for E, sigma in totalXSec ] )

def addSigData( MT, SIG, energyGrid, SigData ) :

    lastNonZero, firstNonZero, xSec = 2, -1, []
    for i1, energy in enumerate( energyGrid ) :
        value = SigData.getValue( energy )
        if( value is None ) : value = 0.
        xSec.append( value )
        if( value != 0. ) :
            if( firstNonZero == -1 ) : firstNonZero = i1
            lastNonZero = i1 + 2
    if( firstNonZero > 0 ) : firstNonZero -= 1
    SIG.append( ( MT, firstNonZero, xSec[firstNonZero:lastNonZero] ) )

class n_nPrimeEnergyData :

    def __init__( self, MT, A, Q ) :

        self.LAW = 3
        self.MT = MT
        self.A = float( A )
        self.Q = float( Q )

    def toACE( self, label, offset, weight, **kwargs ) :

        r1 = self.A / ( 1. + self.A )
        return( [ 0, self.LAW, len( weight ) + 4 ] + weight + [ abs( self.Q ) / r1, r1 * r1 ] )

def processEnergyData( massUnit, neutronMass, energyDatas ) :

    LDLW, MT_DLW, length = [], [], 0
    for MT, EMin, EMax, energyData in energyDatas :
        LDLW.append( length + 1 )
        defaultWeight = [ 0, 2, EMin, EMax, 1.0, 1.0 ]
        label = 'MT = %s for product = "neutron"' % MT
        DLW_ = energyData.toACE( label, length, defaultWeight, massUnit = massUnit, neutronMass = neutronMass )
        MT_DLW.append( ( MT, DLW_ ) )
        length += len( DLW_ )
    return( LDLW, MT_DLW )
