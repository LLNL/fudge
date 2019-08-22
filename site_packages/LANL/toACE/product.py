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

"""
This module adds the method toACE to the product class.
"""

from fudge.gnd import product

def toACE( self, MTData ) :

    from fudge.gnd.productData import distributions
    from fudge.gnd.productData.distributions import angular, energy

    name = self.getName( )
    if( name in [ 'n', 'gamma' ] ) :
        if( self.multiplicity.isConstant ) :
            multiplicity = self.multiplicity.getValue( 0 )
        else :
            multiplicity = self.multiplicity.getNativeData( )
    if( self.getName( ) == 'n' ) :
        nativeData = self.distributions.getNativeData( )
        energyData = None
        if( isinstance( nativeData, distributions.uncorrelated.component ) ) :
            angularForm = nativeData.angularComponent.getNativeData( )
            energyForm = nativeData.energyComponent.getNativeData( )
            energyData = energyForm
        else :
            angularForm = nativeData
        if( isinstance( angularForm, ( angular.isotropic ) ) ) :
            angularData = angularForm
        elif( isinstance( angularForm, ( angular.pointwise, angular.linear, angular.piecewise, angular.LegendrePointwise,
                angular.LegendrePiecewise, angular.mixedRanges ) ) ) :
            angularData = angularForm.toPointwise_withLinearXYs( accuracy = 1e-4 )
        else :
            raise Exception( "angular form '%s' not supported" % angularForm.moniker )
        frame = angularData.getProductFrame( )
        MTData['n'].append( ( frame, multiplicity, angularData, energyData ) )

product.product.toACE = toACE
