# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData.distributions import distribution as distributionModule

#
# component
#
def toENDL( self ) :

    try :
        return( self[0].toENDL( ) )
    except :
        print( 'Distribution "%s" is missing "toENDL" function' % self[0].moniker )
        return( None )

distributionModule.component.toENDL = toENDL
