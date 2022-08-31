# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import copy

from LUPY import enums as enumsModule
from PoPs import IDs as IDsPoPsModule
from PoPs import specialNuclearParticleID as specialNuclearParticleID_PoPsModule

__doc__ = """
This module contains the ReactionProducts class which is used to store a unique list of products for a reaction in a dictionary-like object.
"""

class Category(enumsModule.Enum):

    particle = enumsModule.auto()            # A particle that does not have an outputChannel.
    intermediate = enumsModule.auto()        # A particle with an outputChannel.
    process = enumsModule.auto()             # The process for the reaction.

class Mode(enumsModule.Enum):

    unique = enumsModule.auto()
    products = enumsModule.auto()
    productsWithPhotons = enumsModule.auto()
    familiar = enumsModule.auto()
    residual = enumsModule.auto()

class ReactionProduct :

    def __init__(self, category, multiplicity):

        self.__category = category
        self.__multiplicity = multiplicity

    def __str__( self ) :

        return( 'category=%s multiplicity=%s' % ( self.category, self.multiplicity) )

    __repr__ = __str__

    def __eq__( self, other ) :

        if not isinstance(other, ReactionProduct): raise TypeError('Other must be a ReactionProduct instance: %s.' % other.__class__)

        return self.category == other.category and self.multiplicity == other.multiplicity

    @property
    def category( self ) :

        return( self.__category )

    @property
    def multiplicity( self ) :

        return( self.__multiplicity)

    @multiplicity.setter
    def multiplicity(self, multiplicity):

        self.__multiplicity = multiplicity

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
                particles.append([reactionProduct.category, PoPs[key], key, reactionProduct.multiplicity])
            elif( reactionProduct.category == Category.intermediate ) :
                intermediates.append( [ reactionProduct.category, PoPs[key], key, reactionProduct.multiplicity] )
            elif( reactionProduct.category == Category.process ) :
                processes.append([key, reactionProduct.multiplicity])
            else :
                reaction = reactionProduct.multiplicity

        particles = [ [ key, multiplicity ] for category, dummy, key, multiplicity in sorted( particles ) ]
        intermediates = [ [ key, multiplicity ] for category, dummy, key, multiplicity in sorted( intermediates ) ]

        if len(intermediates) == 1:                 # The following is for some legacy GNDS 1.10 files which have reactions
            if intermediates[0][1] == 2:            # like "n + (O16 -> O16 + photon) [continuum]" which should have been
                intermediates[0][1] = 1             # "n + O16 + photon [continuum]".
                particles.append(intermediates[0])
                intermediates = []

        return( particles, intermediates, processes, reaction )

    def cullPhotons( self ) :
        """
        Returns a copy of self but with photons removed from the list.
        """

        culled = self.copy( )
        if( IDsPoPsModule.photon in culled ) : culled.pop( IDsPoPsModule.photon )
        return( culled )

    def specialReactionLabel(self, outputMode, pops, mode=specialNuclearParticleID_PoPsModule.Mode.familiar, recalled=False):
        """
        Returns the special label for the ReactionProducts that is unique.

        :param outputMode           A value from the Mode enum.
        :param mode:                A value from the special particle id Mode enum.
        :param recalled:            For internal use only. Notes that **specialReactionLabel** is called recursively.

        :return:                    The special label for *self*.
        """

        def update(label, a_sep, products, category):

            sep = ''
            leftToDo = []
            for product in products:
                reactionProduct = self[product]
                if outputMode != Mode.unique and category == Category.intermediate:
                    continue
                if reactionProduct.category != category:
                    leftToDo.append(product)
                    continue
                if len(label) > 0: sep = a_sep
                multiplicity = reactionProduct.multiplicity
                if multiplicity == 1 or category == Category.process:
                    multiplicity = ''
                elif multiplicity == 0:
                    multiplicity = '*'
                pid = specialNuclearParticleID_PoPsModule.specialNuclearParticleID(product, mode)
                if pid in metaStables: pid = metaStables[pid]
                label += '%s%s%s' % (sep, multiplicity, pid)

            return label, leftToDo

        metaStables = {}
        for alias in pops.aliases:
            if alias.isMetaStable: metaStables[alias.pid] = pops.final(alias.pid)

        if not recalled:
            for ReactionProduct in self:                        # Special check for fission.
                if 'Fission' in ReactionProduct: return 'f'

            if len(self) == 2:                                  # Special check for capture.
                ids = list(self.keys())
                if IDsPoPsModule.photon in ids:
                    ids.pop(ids.index(IDsPoPsModule.photon))
                    reactionProduct = self[ids[0]]
                    if reactionProduct.category == Category.particle and reactionProduct.multiplicity == 1:
                        if mode == specialNuclearParticleID_PoPsModule.Mode.familiar: return 'g'
                        return 'photon'
                    
        if outputMode == Mode.familiar or outputMode == Mode.residual:
            residualExluded = ReactionProducts()
            for pid in self:
                if self[pid].category != Category.particle:
                    continue
                residualExluded[pid] = copy.deepcopy(self[pid])
            products = specialNuclearParticleID_PoPsModule.sortLightParticle(residualExluded.keys())
            residual = products[-1]
            if outputMode == Mode.residual: return residual
            multiplicity = self[residual].multiplicity
            if multiplicity > 1:
                residualExluded[residual].multiplicity -= 1
            else:
                residualExluded.pop(residual)
            return residualExluded.specialReactionLabel(Mode.products, pops, mode=mode, recalled=True)

        products = list(self.keys())
        if outputMode != Mode.unique: products = specialNuclearParticleID_PoPsModule.sortLightParticle(products)
        if outputMode != Mode.productsWithPhotons:
            try:
                products.pop(products.index(IDsPoPsModule.photon))
            except:
                pass

        label, products = update('', '+', products, Category.particle)
        if outputMode == Mode.unique:
            label, products = update(label, '++', products, Category.process)
        label, products = update(label, '--', products, Category.intermediate)

        return label
