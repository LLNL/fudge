# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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

multiplicityModule.Unspecified.toENDL = toENDL

#
# Constant1d
#
def toENDL( self ) :

    return( None )

multiplicityModule.Constant1d.toENDL = toENDL

#
# XYs1d
#
def toENDL( self ) :

    return( [ [ x, y ] for x, y in self ] )

multiplicityModule.XYs1d.toENDL = toENDL

#
# Branching1d
#
def toENDL( self ) :

    return( self )

multiplicityModule.Branching1d.toENDL = toENDL

#
# component
#
def toENDL( self ) :

    return( self[0].toENDL( ) )

multiplicityModule.Component.toENDL = toENDL
