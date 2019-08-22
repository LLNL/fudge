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

import sys
import copy
import math

import fudge
from fudge.core.math import *
from fudge.core.utilities import brb
from xData import axes as axesModule
from xData import XYs as XYsModule

from fudge.gnd.productData import energyDeposition as energyDepositionModule
from fudge.gnd.productData import momentumDeposition as momentumDepositionModule

__metaclass__ = type

def calculateDepositionEnergyFromEpP( E, EpP ) :
    "EpP is a list of [ E', P(E') ]"

    axes = axesModule.axes( )
    axes[0] = axesModule.axis( 'a', 0, EpP.axes[0].getUnit( ) )
    axes[1] = axesModule.axis( 'b', 1, EpP.axes[0].getUnit( ) )
    Ep = XYsModule.XYs( data = [ [ EpP[0][0], EpP[0][0] ], [ EpP[-1][0], EpP[-1][0] ] ], axes = axes, accuracy = 1e-12 )
    return( float( EpP.integrateTwoFunctions( Ep ) ) )

def calculateDepositionEnergyFromAngular_angularEnergy( processInfo, derivedStyles, angular, energy, multiplicity, 
        doingGammaMomentum = False, accuracy = 1e-6 ) :

    energyUnit = energy.axes[0].getUnit( )
    energyPrimeUnit = energy.axes[2].getUnit( )
    momentumDepositionUnit = energyPrimeUnit + '/c'

    sqrtP6 = 0.77459666924148337704 / 2.      # sqrt( 0.6 ) / 2
    depEnergy = []
    for indexE, muEpPs in enumerate( energy ) :
        E = muEpPs.value
        I1MuP = angular[indexE]
        sum = 0.
        for indexMu, muEpP in enumerate( muEpPs ) :
            mu2, P2 = I1MuP[indexMu]
            mu2 = muEpP.value
            Ep2 = muEpP.integrateWithWeight_x( )
            if( indexMu != 0 ) :
                muMid = 0.5 * ( mu1 + mu2 )
                EpMid = muEpPs.interpolateAtW( muMid, unitBase = True ).integrateWithWeight_x( )
                if( doingGammaMomentum ) :
                    dMu = sqrtP6 * ( mu2 - mu1 )
                    muG1 = muMid - dMu
                    muEpG1 = muG1 * muEpPs.interpolateAtW( muG1, unitBase = True ).integrateWithWeight_x( )
                    muG2 = muMid + dMu
                    muEpG2 = muG2 * muEpPs.interpolateAtW( muG2, unitBase = True ).integrateWithWeight_x( )
                    sum += ( mu2 - mu1 ) * ( 5. * ( muEpG1 + muEpG2 ) + 8. * muMid * EpMid )
                else :
                    sum += 3. * ( mu2 - mu1 ) * ( P1 * Ep1 + 2 * ( P1 + P2 ) * EpMid + P2 * Ep2 )   # 3 due to 18 instead of 6 below.
            P1 = P2
            mu1 = mu2
            Ep1 = Ep2
        depEnergy.append( [ E, multiplicity.getValue( E ) * sum / 18. ] )

    if( doingGammaMomentum ) : return( depEnergy )

    axes = energyDepositionModule.pointwise.defaultAxes( energyUnit = energyUnit, energyDepositionUnit = energyUnit )
    return( energyDepositionModule.pointwise( data = depEnergy, axes = axes, label = processInfo.style.label, accuracy = accuracy ) )

def GaussQuadrature2( function, parameters, a, b ) :

    if( a == b ) : return( 0. )
    xp, m, width = 0.57735026918962576451, 0.5 * ( a + b ), b - a                 # sqrt( 1. / 3. ), mid-point, width
    x1, x2 = m - 0.5 * width * xp, m + 0.5 * width * xp
    return( 0.5 * width * ( function( x1, parameters ) + function( x2, parameters ) ) )

def GnG_adaptiveQuadrature( function, a, b, parameters, quadrature, tolerance, maxEvaluations = 1000 ) :

    class quadratureInfo :

        def __init__( self, function, parameters, quadrature, estimate, maxEvaluations ) :

            self.function = function
            self.parameters = parameters
            self.quadrature = quadrature
            self.estimate = estimate
            self.evaluations = 0
            self.totalEevaluations = 0
            self.maxEvaluations = maxEvaluations
            self.maxLevelReached = 0

        def GnG_adaptiveQuadrature_2( self, course, a, b, level ) :

            if( a == b ) : return( 0. )
            self.evaluations += 1
            if( self.evaluations > self.maxEvaluations ) : return( 0. )
            level += 1
            if( level > self.maxLevelReached ) : self.maxLevelReached = level
            m = 0.5 * ( a + b )
            l, r = self.quadrature( self.function, self.parameters, a, m ), self.quadrature( self.function, self.parameters, m, b )
            fine = l + r
            extrapolate = ( 16. * fine - course ) / 15.
            if( self.estimate + ( extrapolate - fine ) == self.estimate ) : return( fine )
            return( self.GnG_adaptiveQuadrature_2( l, a, m, level ) + self.GnG_adaptiveQuadrature_2( r, m, b, level ) )

        def __repr__( self ) :

            return( "evaluations = %d, totalEevaluations = %d, maxLevelReached = %d" % ( self.evaluations, self.totalEevaluations,
                self.maxLevelReached ) )

    if( a == b ) : return( 0., None )
    estimate = 0.
    for r in [ 0.9501, 0.2311, 0.6068, 0.4860, 0.8913 ] : estimate += function( a + ( b - a ) * r, parameters )
    estimate = 0.5 * ( estimate * 0.2 * ( b - a ) + quadrature( function, parameters, a, b ) )
    if( estimate == 0 ) : estimate = b - a
    if( tolerance < sys.float_info.epsilon ) : tolerance = sys.float_info.epsilon
    quadInfo = quadratureInfo( function, parameters, quadrature, tolerance * estimate / sys.float_info.epsilon, maxEvaluations )
    course = quadrature( function, parameters, a, b )
    value = quadInfo.GnG_adaptiveQuadrature_2( course, a, b, 0 )
    r = value / estimate
    if( ( r < .1 ) or ( r > 10. ) ) :
        quadInfo.estimate = tolerance * value / sys.float_info.epsilon
        quadInfo.totalEevaluations = quadInfo.evaluations
        quadInfo.evaluations = 0
        quadInfo.maxLevelReached = 0
        value = quadInfo.GnG_adaptiveQuadrature_2( course, a, b, 0 )
    quadInfo.totalEevaluations += quadInfo.evaluations
    return( value, quadInfo )
