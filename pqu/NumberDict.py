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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

"""
Dictionary storing numerical values.
"""

class NumberDict( dict ) :
    """
    An instance of this class acts like an array of numbers with generalized
    (non-integer) indices. A value of zero is assumed for undefined entries.
    NumberDict instances support addition, and subtraction with other NumberDict
    instances, and multiplication and division by scalars.

    This class is used by the PhysicalUnit class to track units and their power.
    For example, for the unit 'm**2 / kg**3', the dictionary items are [('kg', -3), ('m', 2)].
    """

    def __getitem__( self, item ) :

        try:
            return( dict.__getitem__( self, item ) )
        except KeyError:
            return( 0 )

    def __coerce__( self, other ) :

        if( isinstance( other, dict ) ) : other = NumberDict( other )
        other._removeUnneededDimensionlessUnit( )
        return( self, other )

    def __add__( self, other ) :

        sum_dict = NumberDict( )
        for key in self : sum_dict[key] = self[key]
        for key in other : sum_dict[key] = sum_dict[key] + other[key]
        sum_dict._removeUnneededDimensionlessUnit( )
        return( sum_dict )

    def __sub__( self, other ) :

        sum_dict = NumberDict( )
        for key in self : sum_dict[key] = self[key]
        for key in other : sum_dict[key] = sum_dict[key] - other[key]
        sum_dict._removeUnneededDimensionlessUnit( )
        return( sum_dict )

    def __mul__( self, other ) :

        new = NumberDict( )
        for key in self : new[key] = other * self[key]
        new._removeUnneededDimensionlessUnit( )
        return( new )

    __rmul__ = __mul__

    def __div__( self, other ) :

        new = NumberDict( )
        for key in self : new[key] = self[key] / other
        new._removeUnneededDimensionlessUnit( )
        return( new )

    def __deepcopy__( self ) :

        other = NumberDict( )
        for key in self : other[key] = self[key]
        return( other )

    __copy__ = __deepcopy__

    def _removeUnneededDimensionlessUnit( self ) :
        """For internal use only."""

        if( ( '' in self ) and ( len( self ) > 1 ) ) : del self['']
