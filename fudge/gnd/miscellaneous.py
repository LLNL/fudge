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

from fudge.core.math import *
from xData import standards as standardsModule
from xData import axes as axesModule

__metaclass__ = type

def setStringMember( self, memberName, value ) :

    if( not( isinstance( value, str ) ) ) : TypeError( 'value for member %s is not a string' % memberName )
    exec( 'self.%s = "%s"' % ( memberName, value ) )

def find_gndPath( reactionSuite, gndPath ) :

    def find_gndPath2( parent, paths, gndPath ) :

        if( len( paths ) == 0 ) : return( parent )
        current, others = paths[0], paths[1:]
        if( current[:6] == 'label:' ) :
            label = current[6:]
            for element in parent :
                if( element.getLabel( ) == label ) : return( find_gndPath2( element, others, gndPath ) )
            raise Exception( 'Could not find path fragment = %s in gndPath = %s' % ( '/'.join( paths ), gndPath ) )
        else :
            raise Exception( 'Finding element currently not supported for gndPath = %s' % gndPath )

    return( find_gndPath2( reactionSuite, gndPath.split( '/' )[1:], gndPath ) )

def floatToString( f, precision = 12 ) :
    """Returns the shortest string for the float from %f and %e conversion at precision."""

    s = ( '%%.%df' % precision ) % f
    s = s.rstrip( '0' )
    t = ( '%%.%de' % precision )  % f
    t2 = t.split( 'e' )
    p = t2[0].rstrip( '0' )
    t3 = p + 'e' + t2[1]
    e = int( t2[1] )
    if( ( e < 0 ) and ( ( len( p ) - 2 - e ) > precision ) ) : s = t3
    if( len( t3 ) <= len( s ) ) : s = t3
    return( s )

def TMs2Form( processInfo, tempInfo, newComponents, TM_1, TM_E, axes_p ) :

    from productData import distributions
    from fudge.core.math import matrix
    crossSectionUnit = 'b'                          # ?????? 'b' should not be hardwired.
    axes = axesModule.axes( rank = 3 )
    axes[0] = axesModule.axis( 'C_l(energy_in,energy_out)', 0, crossSectionUnit )
    axes[1] = axesModule.axis( 'energy_out',          1, axes_p[0].unit )
    axes[2] = axesModule.axis( axes_p[-1].getLabel( ), 2, axes_p[-1].unit )
    # convert TM_1 and TM_E from dicts into matrix objects:
    if( not ( TM_1 is None ) ) :
        TM_1_new = []
        for i1 in range(len(TM_1)):
            TM_1_new.append( matrix.matrix( [TM_1[i1][i2] for i2 in range(len(TM_1[i1]))], form="sparse_asymmetric" ) )
        component = distributions.Legendre.component( )
        component.addForm( distributions.Legendre.grouped( axes, TM_1_new, standardsModule.frames.labToken ) )
        newComponents.append( component )
    if( not ( TM_E is None ) ) :
        TM_E_new = []
        for i1 in range(len(TM_E)):
            TM_E_new.append( matrix.matrix( [TM_E[i1][i2] for i2 in range(len(TM_E[i1]))], form="sparse_asymmetric" ) )
        component = distributions.Legendre.energyConservationComponent( )
        component.addForm( distributions.Legendre.energyConservationGrouped( axes, TM_E_new, standardsModule.frames.labToken ) )
        newComponents.append( component )

def makeGrouped( self, processInfo, tempInfo, data, normType = 'groupedCrossSectionNorm' ) :

    from fudge.processing import miscellaneous

    projectile, target = processInfo.getProjectileName( ), processInfo.getTargetName( )
    norm = tempInfo[normType]
    crossSection = tempInfo['crossSection']
    label, unit = data.axes[1].label, data.axes[1].unit
    axes_ = crossSection.axes.copy( )
    if( normType == 'groupedFlux' ) :
        axes_[1].label = "%s %s" % ( label,  axes_[1].label )
        axes_[1].unit = "%s %s" % ( unit,  axes_[1].unit )
    else :
        axes_[1].label = "%s" % label
        axes_[1].unit = "%s" % unit
    independent, dependent, qualifier = axes_[0].interpolation.getInterpolationTokens( )
    axes_[0].interpolation = axes.interpolationXY( independent, standardsModule.interpolation.flatToken )
    grouped = miscellaneous.groupTwoFunctionsAndFlux( projectile, processInfo, crossSection, data, norm = norm )
    return( axes_, grouped )
