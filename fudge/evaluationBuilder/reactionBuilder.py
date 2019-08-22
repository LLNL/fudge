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

from . import crossSectionBuilder, productBuilder
from fudge import gnd
from pqu import PQU

class reactionBuilder:
    def __init__(self, crossSection=None, products=[]):
        if crossSection: self.addCrossSection( crossSection )
        else: self.crossSection = None
        if products:
            [self.addProduct(p) for p in products]
        else: self.products = []
        self.documentation = {}
        self.residualInfo = ()

    def addCrossSection( self, crossSection ):
        if isinstance( crossSection, crossSectionBuilder.crossSectionBuilder ):
            self.crossSection = crossSection.component
        elif isinstance( crossSection, gnd.reactionData.crossSection.component ):
            self.crossSection = crossSection
        else:
            raise Exception("addCrossSection() requires an evaluator.crossSectionBuilder instance")

    def addProduct( self, product ):
        if isinstance( product, productBuilder.productBuilder ):
            self.products.append( product.product )
        elif isinstance( product, gnd.product.product ):
            self.products.append( product )
        else:
            raise Exception("addProduct() requires a gnd.product or evaluator.productBuilder instance")

    def addResidual( self, excitationLevel = None, excitationEnergy = None ):
        self.residualInfo = (excitationLevel, excitationEnergy)

    def addDocumentation( self, label, documentation ):
        self.documentation = gnd.documentation.documentation( label, documentation )

    def finalize( self, reactionSuite, particleGenerator ):

        import datetime
        date = datetime.datetime.today().strftime("%Y-%m-%d")

        def calculateZA( ZA1, ZA2 ):    # compute the residual
            minus = (ZA2 < 0)
            ZA2 = abs( ZA2 )
            if( ( ZA1 % 1000 ) == 0 ): ZA2 = 1000 * ( ZA2 // 1000 )
            if minus: return ZA1 - ZA2
            return ZA1 + ZA2

        # add the residual particle (as computed from ZA of target, projectile and products)
        # If the reaction is 2-body, set the residual's distribution to 'recoil', otherwise 'none'

        try:
            twoBody = ( len( self.products ) == 1 and ( self.products[0].multiplicity.getConstant() == 1 ) )
        except: # in case multiplicity not constant
            twoBody = False

        ZAresidual = calculateZA( reactionSuite.target.getZ_A_SuffixAndZA()[-1],
                reactionSuite.projectile.getZ_A_SuffixAndZA()[-1] )
        for prod in self.products:
            if prod.getName() == 'gamma': continue
            ZAproduct = prod.particle.getZ_A_SuffixAndZA()[-1] * prod.multiplicity.getConstant()
            ZAresidual = calculateZA( ZAresidual, -ZAproduct )

        if self.residualInfo:
            from fudge.particles import nuclear
            residualName = nuclear.nucleusNameFromZA( ZAresidual )
            level, energy = self.residualInfo
            if level == 'continuum':
                residualName += '_c'
                level = 'c'
            residual = particleGenerator.newProduct( residualName, levelIndex = level,
                    levelEnergy = energy, multiplicity = 1 )
            if twoBody:
                residual.addDistribution( 'Recoil', self.products[0] )
            self.addProduct( residual )
            self.needResidual = False

        if twoBody:
            channel_ = gnd.channels.twoBodyOutputChannel( )
        else:
            channel_ = gnd.channels.NBodyOutputChannel( )
        self.reaction = gnd.reactions.reaction.reaction( channel_, label="-1", ENDF_MT=-1, attributes={'date':date} )
        self.reaction.setCrossSection( self.crossSection )

        # ensure each product has a unique label:
        from fudge.legacy.converting import endlToGND
        labler = endlToGND.tokenAndIndices()
        for prod in self.products:
            prod.label = labler.getLabel( prod.getToken() )

        for prod_ in self.products:
            self.reaction.outputChannel.addProduct( prod_ )

        if ZAresidual != 0 and not self.residualInfo:
            raise Exception("ZA doesn't balance for reaction with outputChannel %s! Was the residual added?" % self.reaction)

        return self.reaction

class summedReactionBuilder( reactionBuilder ):
    def __init__(self, name, crossSection=None):
        self.name = name
        if crossSection: self.addCrossSection( crossSection )
        else: self.crossSection = None
        self.documentation = {}

    def addProduct(self): raise Exception("No products allowed for summed channel")

    def finalize(self, reactionSuite, particleGenerator):
        pass

class fissionReactionBuilder( reactionBuilder ):
    def __init__(self, crossSection=None, products=[]):
        pass

