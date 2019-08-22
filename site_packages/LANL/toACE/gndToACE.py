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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

import time
from fudge.core.utilities import brb
from fudge.core.math.xData import axes, XYs
from fudge.gnd.productData import distributions
from fudge.legacy.converting import endf_endl
from pqu import physicalQuantityWithUncertainty

floatFormat = '%20.12E'
floatFormat = ' %19.11E'

def toACE( self, fileName, version, temperature, data ) :

    projectile, target = self.projectile, self.target
    targetZ, targetA, targetSuffix, targetZA = target.getZ_A_SuffixAndZA( )
    neutronMass = projectile.getMass( 'eV/c**2' )
    targetMass = target.getMass( 'eV/c**2' )
    if( temperature.isTemperature( temperature ) ) :
        temperature_MeV = temperature.getValueAs( 'MeV / k' )
        temperature_K = temperature.getValueAs( 'K' )
    else :
        temperature_MeV = temperature.getValueAs( 'MeV' )
        temperature_K = temperature.getValueAs( 'k * K' )
    processingTime = time.localtime( )
    strRecords = [ "%6d.%.2dc%12.7f %11.5E %.2d/%.2d/%.2d" % ( targetZA, version, targetMass / neutronMass, temperature_MeV, processingTime[1],
        processingTime[2], processingTime[0] ) ]

    MAT = endf_endl.ZAAndMATFromParticleName( target.getName( ) )[1]
    MAT_ID = '  mat %d' % MAT
    HK = '%2d-%s at %.1fK from %s-%s Fudge.toACE' % \
        ( targetZ, target.getName( ), temperature_K, self.styles['evaluated'].attributes['library'], self.styles['evaluated'].attributes['version'] )
    HK = "%-70s" % HK[:70]
    strRecords.append( HK + MAT_ID )
    for i in xrange( 4 ) : strRecords.append( '      0   0.000000      0   0.000000      0   0.000000      0   0.000000' )

    NXS = 16 * [ 0 ]
    JXS = 32 * [ 0 ]
    MRT, NU_prompt, NU_total, LQR, TYP, SigData, neutronAngular, neutronEnergy = [], [], [], [], [], {}, [], []
    NXS[2-1] = targetZA

    for MT, MTData in data :
        if( MT == 2 ) :                 # Elastic (MT = 2) must always be present in ACE file.
            EMin = MTData['ESZ'].domainMin( unitTo = 'MeV' )
            break
    totalXSec = XYs.XYs( axes.defaultAxes( labelsUnits = { 0 : ( '', 'MeV' ), 1 : ( '', 'b' ) } ), [], 1e-3 )
    absorptionXSec = XYs.XYs( axes.defaultAxes( labelsUnits = { 0 : ( '', 'MeV' ), 1 : ( '', 'b' ) } ), [], 1e-3 )
    for MT, MTData in data :
        XSec = MTData['ESZ']
        XSec = XSec.convertAxisToUnit( 0, 'MeV' )
        XSec = XSec.convertAxisToUnit( 1, 'b' )
        neutronData = MTData['n']
        if( MT == 2 ) :
            elasticXSec = XSec
            if( len( neutronData ) > 1 ) : raise Exception( 'Only one type of neutron data is support except for fission.' )
            if( len( neutronData ) > 0 ) : neutronAngular.append( ( MT, neutronData[0][2] ) )
        else :
            MRT.append( MT )
            LQR.append( MTData['Q'] )
            if( XSec[0][0] != 0 ) :
                if( XSec[0][0] > EMin ) : XSec = XSec.dullEdges( lowerEps = 1e-8 )
            SigData[MT] = XSec
            totalXSec = totalXSec + XSec
            if( len( neutronData ) == 0 ) :
                absorptionXSec = absorptionXSec + XSec
                TYP.append( 0 )
            else :
                NXS[5-1] = 1
                if( MTData['fission'] ) :
                    raise 'hell: need to implement fission'
                else :
                    if( len( neutronData ) > 1 ) : raise Exception( 'Only one type of neutron data is support except for fission.' )
                    frame, multiplicity, angularData, energyData = neutronData[0]
                    if( energyData is None ) :
                        if( frame == axes.centerOfMassToken ) : multiplicity *= -1
                    else :
                        if( isinstance( energyData, distributions.energy.NBodyPhaseSpace ) ) :
                            angularData = None
                        else :
                            raise Exception( 'energyData needs to be implemented' )
                    TYP.append( multiplicity )
                    neutronAngular.append( ( MT, angularData ) )
                    if( energyData is not None ) : neutronEnergy.append( energyData )
    totalXSec = elasticXSec + totalXSec
    annotates, XSS, energyGrid, totalSigma = [], [], [], []
    for E, sigma in totalXSec :
        energyGrid.append( E )
        totalSigma.append( sigma )

# 1) Add the ESZ block.
    NXS[3-1] = len( energyGrid )
    updateXSSInfo( 'energyGrid', annotates, XSS, energyGrid )
    updateXSSInfo( 'totalSigma', annotates, XSS, totalSigma )
    updateXSSInfo( 'absorption cross section', annotates, XSS, mapEnergyToTotal( totalXSec, absorptionXSec ))
    updateXSSInfo( 'elastic cross section', annotates, XSS, mapEnergyToTotal( totalXSec, elasticXSec ) )
    averageHeating = len( energyGrid ) * [ 0. ]
    updateXSSInfo( 'average heating', annotates, XSS, averageHeating )

# 2) Add the NU block.
    if( NU_prompt != [] ) :
        if( NU_total != [] ) :
            pass
    elif( NU_total != [] ) :
        pass
    if( ( len( NU_prompt ) + len( NU_total ) ) > 0 ) : JXS[2-1] = len( XSS ) + 1
    updateXSSInfo( 'NU', annotates, XSS, NU_prompt + NU_total )

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
    for MT, MTData in data :
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
    LAND, AND, length = [], [], 0
    for MT, angular in neutronAngular :
        AND_MT = []
        if( angular is None ) :
            LAND.append( -1 )
        elif( angular.isIsotropic( ) ) :
            LAND.append( 0 )
        elif( isinstance( angular, ( distributions.angular.linear ) ) ) :
            LAND.append( length + 1 )
            length += 2 * len( angular ) + 1
            energies_in = [ len( angular ) ]
            LCs, Ps = [], []
            for xys in angular :
                energy_in = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( xys.value, angular.axes[0].getUnit( ) ).getValueAs( 'MeV' )
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
            AND.append( ( MT, energies_in + LCs + Ps ) )
        else :
            raise Exception( 'Unsupport neutron angular distribution type = %s' % type( angular ) )
    JXS[8-1] = len( XSS ) + 1
    updateXSSInfo( 'LAND', annotates, XSS, LAND )
    JXS[9-1] = len( XSS ) + 1
    for MT, AND_MT in AND : updateXSSInfo( 'AND(MT=%s)' % MT, annotates, XSS, AND_MT )

# Time to wrap it all up.
    NXS[1-1] = len( XSS )
    JXS[1-1] = 1

    strRecords += intArrayToRecords( NXS )
    strRecords += intArrayToRecords( JXS )
    strRecords += XSSToStrings( annotates, XSS )

    strRecords.append( '' )
    f = open( fileName, 'w' )
    f.write( '\n'.join( strRecords ) )
    f.close( )        

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

def XSSToStrings( annotates, XSS ) :

    strData, record, i3 = [], [], 0
    annotates.append( ( len( XSS ), '' ) )
    i2, label = annotates[i3]
    for i1, datum in enumerate( XSS ) :
        if( type( datum ) == int ) :
            record.append( "%20d" % datum )
        else :
            record.append( floatFormat % datum )
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
