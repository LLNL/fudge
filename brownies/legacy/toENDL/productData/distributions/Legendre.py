# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData.distributions import Legendre as LegendreModule

def toENDL( self ) :

    data = []
    for order in self.data :
        orderData = []
        for xyz in order :
            orderData.append( [ xyz.outerDomainValue, [ [ x, y ] for x, y in xyz ] ] )
        data.append( [ order.outerDomainValue, orderData ] )
    return( { 4 : data } )

LegendreModule.LLNLPointwise.toENDL = toENDL

#
# form
#
def toENDL( self ) :

    return( self.Legendre.toENDL( ) )

LegendreModule.Form.toENDL = toENDL
