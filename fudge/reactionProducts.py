# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from PoPs import IDs as IDsPoPsModule

__doc__ = """
This module contains the ReactionProducts class which is used to store a unique list of products for a reaction in a dictionary-like object.
"""

class Category :

    particle = 0
    intermediate = 1
    process = 2
    reaction = 3

class ReactionProduct :

    def __init__( self, _category, _value ) :

        self.__category = _category
        self.__value = _value

    def __str__( self ) :

        return( 'c=%s : v=%s' % ( self.category, self.value ) )

    __repr__ = __str__

    def __eq__( self, other ) :

        if( not( isinstance( other, ReactionProduct ) ) ) : raise TypeError( 'Other must be a ReactionProduct instance: %s.' % other.__class__ )

        return( ( self.category == other.category ) and ( self.value == other.value ) )

    @property
    def category( self ) :

        return( self.__category )

    @property
    def value( self ) :

        return( self.__value )

    @value.setter
    def value( self, _value ) :

        self.__value = _value

class ReactionProducts( dict ) :
    """
    Class used to store a unique list of products for a reaction in a dictionary-like object. The key for each item in the dictionary
    is a product's particle id or an outputChannel's process.
    """

    def asSortedList( self, PoPs, excludePhotons = True ) :
        """
        """

        particles = []
        intermediates = []
        processes = []
        reaction = None
        for key in self :
            if( excludePhotons and ( key == IDsPoPsModule.photon ) ) : continue
            reactionProduct = self[key]
            if( reactionProduct.category == Category.particle ) :
                particles.append( [ reactionProduct.category, PoPs[key], key, reactionProduct.value ] )
            elif( reactionProduct.category == Category.intermediate ) :
                intermediates.append( [ reactionProduct.category, PoPs[key], key, reactionProduct.value ] )
            elif( reactionProduct.category == Category.process ) :
                processes.append( [ key, reactionProduct.value ] )
            else :
                reaction = reactionProduct.value

        particles = [ [ key, value ] for category, dummy, key, value in sorted( particles ) ]
        intermediates = [ [ key, value ] for category, dummy, key, value in sorted( intermediates ) ]

        return( particles, intermediates, processes, reaction )

    def cullPhotons( self ) :
        """
        Returns a copy of self but with photons removed from the list.
        """

        culled = self.copy( )
        if( IDsPoPsModule.photon in culled ) : culled.pop( IDsPoPsModule.photon )
        return( culled )
