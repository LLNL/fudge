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

def shiftFloatABit( f, s, r_eps, a_eps, z_eps ) :
    "Only for internal use by shiftFloatDownABit and shiftFloatUpABit."

    if( a_eps is None ) : a_eps = r_eps
    if( z_eps is None ) : z_eps = r_eps
    if( r_eps < 0. ) : raise Exception( 'r_eps = %s < 0.' % r_eps )
    if( r_eps > 0.5 ) : raise Exception( 'r_eps = %s > 0.5' % r_eps )
    if( a_eps < 0. ) : raise Exception( 'a_eps = %s < 0.' % a_eps )
    if( a_eps > 0.5 ) : raise Exception( 'a_eps = %s > 0.5' % a_eps )
    if( z_eps < 0. ) : raise Exception( 'z_eps = %s < 0.' % z_eps )
    if( z_eps > 0.5 ) : raise Exception( 'z_eps = %s > 0.5' % z_eps )
    if( f < 0. ) : r_eps = -r_eps
    if( f == 0. ) : return( s * z_eps )
    return( ( 1 + s * r_eps ) * f + s * a_eps )

def shiftFloatDownABit( f, r_eps, a_eps = 0., z_eps = None ) :
    """
    Returns a float that is slightly less than f. The amount less depends on r_eps, a_eps and z_eps.
    If f is 0., returns -z_eps; otherwise returns ( 1 - r_eps ) * f - a_eps. r_eps, a_eps and z_eps
    all must be between 0. and 0.5 inclusive."""

    return( shiftFloatABit( f, -1, r_eps, a_eps, z_eps ) )

def shiftFloatUpABit( f, r_eps, a_eps = 0., z_eps = None ) :
    """
    Returns a float that is slightly more than f. The amount more depends on r_eps, a_eps and z_eps.
    If f is 0., returns z_eps; otherwise returns ( 1 + r_eps ) * f + a_eps. r_eps, a_eps and z_eps
    all must be between 0. and 0.5 inclusive."""

    return( shiftFloatABit( f, 1, r_eps, a_eps, z_eps ) )
