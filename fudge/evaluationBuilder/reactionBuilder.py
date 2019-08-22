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
            if( prod.name == 'gamma' ) : continue
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
        from fudge.legacy.converting import toGNDMisc
        labler = toGNDMisc.tokenAndIndices()
        for prod in self.products:
            prod.label = labler.getLabel( prod.name )

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

