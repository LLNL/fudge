# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Module for the 'reference' class for distributions."""

import xData.link as linkModule

from . import base as baseModule

__metaclass__ = type

class form( linkModule.link, baseModule.form ) :

    moniker = 'reference'
    subformAttributes = []

    def __init__( self, link = None, root = None, path = None, label = None, relative = False ) :

        linkModule.link.__init__( self, link = link, root = root, path = path, label = label, relative = relative )
        baseModule.form.__init__( self, label, None, [] )

    @property
    def referenceInstance( self ):

        if self.link is None: raise Exception("Unresolved link!")
        return self.link

    @property
    def productFrame( self ): return self.referenceInstance.productFrame

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        return( self.referenceInstance.calculateAverageProductData( style, indent = indent, **kwargs ) )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        pass

    def processMC_cdf( self, style, tempInfo, indent ) :

        # temporary solution:
        return( self.referenceInstance.processMC_cdf( style, tempInfo, indent ) )

        # better solution: add another link pointing to the processed version of what this points to:
        """
        newReference = form( label=style.label, relative=True )
        self.ancestor.add( newReference )

        tempInfo['brokenLinks'].append( [self, newReference] )
        """

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.referenceInstance.processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.referenceInstance.toPointwise_withLinearXYs( **kwargs ) )

class CoulombPlusNuclearElastic( form ) :

    moniker = 'CoulombPlusNuclearElastic'

class thermalNeutronScatteringLaw( form ) :

    moniker = 'thermalNeutronScatteringLaw'
