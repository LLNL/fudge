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

import math

from xData import axes as axesModule
from xData import XYs as XYsModule

def calculatePofMu( energy, index, SofE, epsilon ) :

    def addDeltaFunction( PofMu, energy, energy_i, S_i, S_im1, epsilon ) :

        dS = S_i - S_im1
        if( dS == 0 ) : return
        mu_i = 1 - 2 * energy_i / energy
        delta = abs( epsilon * mu_i )
        if( delta == 0 ) : delta = epsilon * 1e-3
        mu_lower = mu_i - delta
        mu_upper = mu_i + delta
        if( mu_lower < -1 ) :
            mu_lower = -1
            mu_i = mu_lower + delta
            mu_upper = mu_i + delta
        elif( mu_upper > 1 ) :
            mu_upper = 1
            mu_i = mu_upper - delta
            mu_lower = mu_i - delta

        if( len( PofMu ) > 0 ) :
            if( PofMu[-1][0] <= mu_lower ) :
                raise ValueError( 'Delta function overlap at energy = %s and energy_i = %s: %s and %s' % 
                        ( energy, energy_i, PofMu[-1][0], mu_i ) )
        PofMu.insert( 0, [ mu_upper, 0 ] )
        PofMu.insert( 0, [ mu_i,     dS / delta ] )
        PofMu.insert( 0, [ mu_lower, 0 ] )

    if( index == 0 ) : index = 1

    PofMu = []
    S_im1 = 0
    for i1 in range( index ) :
        energy_i, S_i = SofE[i1]
        addDeltaFunction( PofMu, energy, energy_i, S_i, S_im1, epsilon )
        S_im1 = S_i

    if( len( PofMu ) > 0 ) :
        if( PofMu[0][0]  > -1 ) : PofMu.insert( 0, [ -1, 0 ] )
        if( PofMu[-1][0] <  1 ) : PofMu.append(    [  1, 0 ] )
    return( PofMu )

def bisectionTest( I0Fine, x1, y1, x2, y2, accuracy, level = 0 ) :

    y2p = y1 * x1 / x2
    xMid = 0.5 * ( x1 + x2 )
    yMid = y1 * x1 / xMid
    yEstimate = 0.5 * ( y1 + y2p )
    if( abs( yEstimate - yMid ) <= accuracy * yMid ) : return
    bisectionTest( I0Fine, x1,   y1,   xMid, yMid, accuracy, level + 1 )
    I0Fine.append( [ xMid, yMid ] )
    bisectionTest( I0Fine, xMid, yMid, x2,   y2,   accuracy, level + 1 )

def toENDL( coherentElastic, energyMin_MeV, energyMax_MeV, temperature_MeV, accuracy = 0.001, epsilon = 1e-6 ) :

    S_table = coherentElastic.S_table
    gridded2d = S_table.gridded2d

    array = gridded2d.array.constructArray( )

    energyGrid = gridded2d.axes[1].copy( [] )

    temperatureGrid = gridded2d.axes[2].copy( [] )
    temperatureGrid.convertToUnit( 'MeV/k' )
    for index2, temperature2 in enumerate( temperatureGrid.values ) :
        if( temperature2 >= temperature_MeV ) : break
    if( temperature_MeV > temperature2 ) :
        SofE = XYsModule.XYs1d( ( energyGrid.values, array[-1] ), dataForm = 'XsAndYs', interpolation = energyGrid.interpolation )
    elif( ( index2 == 0 ) or ( temperature2 == temperature_MeV ) ) :
        SofE = XYsModule.XYs1d( ( energyGrid.values, array[index2] ), dataForm = 'XsAndYs', interpolation = energyGrid.interpolation )
    else :
        index1 = index2 - 1
        temperature1 = temperatureGrid.values[index1]
        fraction = ( temperature_MeV - temperature1 ) / ( temperature2 - temperature1 )

        SofE1 = XYsModule.XYs1d( ( energyGrid.values, array[index1] ), dataForm = 'XsAndYs', interpolation = energyGrid.interpolation )
        SofE2 = XYsModule.XYs1d( ( energyGrid.values, array[index2] ), dataForm = 'XsAndYs', interpolation = energyGrid.interpolation )
        SofE = ( 1 - fraction ) * SofE1 + fraction * SofE2

    axes = axesModule.axes( )
    axes[0].label = gridded2d.axes[0].label
    axes[0].unit = gridded2d.axes[0].unit
    axes[1].label = energyGrid.label
    axes[1].unit = energyGrid.unit
    SofE.axes = axes
    SofE.convertUnits( { "eV" : "MeV" } )

    _SofE = SofE
    SofE = []
    for energy, S2 in _SofE :
        if( energy > energyMax_MeV ) : break
        SofE.append( [ energy, S2 ] )
    if( ( SofE[-1][0] < energyMax_MeV ) and ( SofE[-1][0] < _SofE[-1][0] ) ) :
        eMax = _SofE[-1][0]
        if( energyMax_MeV < eMax ) : eMax = energyMax_MeV
        SofE.append( [ eMax, SofE[-1][1] ] )

    S1 = 0
    I0 = []
    I1 = []
    n_minus1 = len( SofE ) - 1
    for index, energyS in enumerate( SofE ) :
        energy, S2 = energyS

        energy1 = ( 1 - epsilon ) * energy
        if( index == 0 ) :
            S1 = S2
            energy1 = energy

        energy2 = ( 1 + epsilon ) * energy
        if( index == n_minus1 ) : energy2 = energy

        if( energy >= energyMax_MeV ) :
            energy = energyMax_MeV
            energy2 = energyMax_MeV
        if( energy2 > energyMax_MeV ) :
            energy2 = energyMax_MeV

        I0.append( [ energy1, S1 / energy ] )
        if( energy != energy1 ) : I0.append( [ energy2, S2 / energy ] )

        PofMu = calculatePofMu( energy1, index,     SofE, epsilon )
        if( len( PofMu ) > 0 ) : I1.append( [ energy1, PofMu ] )
        if( energy != energy1 ) :
            PofMu = calculatePofMu( energy2, index + 1, SofE, epsilon )
            if( len( PofMu ) > 0 ) : I1.append( [ energy2, PofMu ] )

        S1 = S2
        if( energy == energyMax_MeV ) : break

    accuracy = min( 0.1, max( 1e-5, accuracy ) )
    I0Fine = [ I0.pop( 0 ) ]
    x1, y1 = I0Fine[0]
    for x2, y2 in I0 :
        if( ( x2 - x1 ) > ( 2 * epsilon * x1 ) ) : bisectionTest( I0Fine, x1, y1, x2, y2, accuracy )
        I0Fine.append( [ x2, y2 ] )
        x1, y1 = x2, y2

    return( I0Fine, I1, None )
