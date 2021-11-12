# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData.distributions import uncorrelated as uncorrelatedModule

def toENDL( self ) :

    return( self.data.toENDL( ) )

uncorrelatedModule.angularSubform.toENDL = toENDL

def toENDL( self ) :

    return( self.data.toENDL( ) )

uncorrelatedModule.energySubform.toENDL = toENDL

#
# form
#
def toENDL( self ) :

    I1And4 = { 1 : self.angularSubform.toENDL( ) }
    I1And4[4] = [ [ 0, self.energySubform.toENDL( ) ] ]
    return( I1And4 )

uncorrelatedModule.form.toENDL = toENDL
