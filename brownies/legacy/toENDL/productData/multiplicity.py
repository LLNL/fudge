# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.productData import multiplicity as multiplicityModule


#
# unspecified
#
def toENDL( self ) :

    return( None )

multiplicityModule.unspecified.toENDL = toENDL

#
# constant1d
#
def toENDL( self ) :

    return( None )

multiplicityModule.constant1d.toENDL = toENDL

#
# XYs1d
#
def toENDL( self ) :

    return( [ [ x, y ] for x, y in self ] )

multiplicityModule.XYs1d.toENDL = toENDL

#
# branching1d
#
def toENDL( self ) :

    return( self )

multiplicityModule.branching1d.toENDL = toENDL

#
# component
#
def toENDL( self ) :

    return( self[0].toENDL( ) )

multiplicityModule.component.toENDL = toENDL
