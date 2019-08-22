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

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    """
    mimics the math.isclose function introduced in Python3
    """
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def shiftFloatABit( f, s, r_eps, a_eps, z_eps ) :
    """Only for internal use by shiftFloatDownABit and shiftFloatUpABit."""

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
