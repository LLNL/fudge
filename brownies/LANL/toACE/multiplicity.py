# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module adds the method toACE to the classes in the fudge.productData.multiplicity module.
"""

from xData import enums as xDataEnumsModule

from fudge.productData import multiplicity as multiplicityModule


#
#   Polynomial1d multiplicity.
#
def toACE( self ) :

    factor = 1
    coefficients = []
    for coefficient in self :
        coefficients.append( coefficient / factor )
    return( [ 1, len( self ) ] + coefficients )

multiplicityModule.Polynomial1d.toACE = toACE

#
#   XYs1d multiplicity.
#
def toACE( self ) :

    interpolation = 0
    if self.interpolation == xDataEnumsModule.Interpolation.flat:
        interpolation = 1
    elif self.interpolation == xDataEnumsModule.Interpolation.linlin:
        interpolation = 2
    if( interpolation == 0 ) : raise Exception( 'Interpolation "%s" not supported' % self.interpolation )

    return( [ 2, 0, len( self ) ] + [ E1 for E1, m1 in self ] + [ m1 for E1, m1 in self ] )

multiplicityModule.XYs1d.toACE = toACE
