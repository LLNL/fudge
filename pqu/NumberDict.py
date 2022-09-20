# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
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

    def __truediv__( self, other ) :

        new = NumberDict( )
        for key in self : new[key] = self[key] / other
        new._removeUnneededDimensionlessUnit( )
        return( new )

    __div__ = __truediv__

    def _removeUnneededDimensionlessUnit( self ) :
        """For internal use only."""

        if( ( '' in self ) and ( len( self ) > 1 ) ) : del self['']
