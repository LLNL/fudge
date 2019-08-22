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
This module adds the method toACE to the classes in the fudge.gnd.product module.
"""

from fudge.gnd import product
from fudge.gnd.productData import distributions
from xData import axes
import angularEnergy

def toACE( self, MTData, MT ) :

    from fudge.gnd.productData import distributions
    from fudge.gnd.productData.distributions import angular

    name = self.name
    if( name in [ 'n', 'gamma' ] ) :
        if( self.multiplicity.isConstant( ) ) :
            multiplicity = self.multiplicity.getValue( 0 )
        else :
            multiplicity = self.multiplicity.getNativeData( )
    if( self.name == 'n' ) :
        nativeDataToken = self.distributions.getNativeDataToken( )
        angularForm, energyData = None, None
        if( ( MT == 18 ) and ( nativeDataToken == distributions.base.unspecifiedComponentToken ) ) :
            pass
        else :
            nativeData = self.distributions.getNativeData( )
            if( nativeDataToken == distributions.base.angularComponentToken ) :
                angularForm = nativeData
            elif( nativeDataToken == distributions.base.uncorrelatedComponentToken ) :
                angularForm = nativeData.angularComponent.getNativeData( )
                energyData = nativeData.energyComponent.getNativeData( )
            elif( nativeDataToken == distributions.base.energyAngularComponentToken ) :
                if( isinstance( nativeData, distributions.energyAngular.KalbachMann ) ) :
                    energyData = nativeData
                else :
                    raise Exception( '%s is not implemented' % nativeDataToken )
            elif( nativeDataToken == distributions.base.LegendreComponentToken ) :
                accuracy = 1e-3
                if( nativeData.axes[1].isFlat( qualifierOk = True ) ) : accuracy = 1e-6
                energyData = nativeData.toPointwise_withLinearXYs( accuracy )
            elif( nativeDataToken == distributions.base.angularEnergyComponentToken ) :
                angularForm = angularEnergy.angularFor_angularEnergy( nativeData )
                energyData = nativeData
            else :
                print 'token', self.distributions.getNativeDataToken( )
                raise Exception( 'Unsupported nativeData = %s' % nativeData )

        if( angularForm is None ) :
            angularData = angularForm
        elif( isinstance( angularForm, distributions.angularEnergy.pointwise ) ) :
            angularData = angularForm
        elif( isinstance( angularForm, angular.isotropic ) ) :
            angularData = angularForm
        elif( isinstance( angularForm, angularEnergy.angularFor_angularEnergy ) ) :
            angularData = angularForm
        elif( isinstance( angularForm, ( angular.pointwise, angular.linear, angular.piecewise, angular.LegendrePointwise,
                angular.LegendrePiecewise, angular.mixedRanges ) ) ) :
            angularData = angularForm.toPointwise_withLinearXYs( accuracy = 1e-4 )
        else :
            from fudge.core.utilities import brb
            brb.objectoutline( angularForm )
            raise Exception( "angular form '%s' not supported" % angularForm.moniker )

        frame = axes.labToken
        if( angularData is not None ) :
            frame = angularData.getProductFrame( )
        else :
            if( energyData is not None ) : frame = energyData.getProductFrame( )
        MTData['n'].append( { 'product' : self, 'frame' : frame, 'multiplicity' : multiplicity, 
            'angularData' : angularData, 'energyData' : energyData } )

product.product.toACE = toACE
