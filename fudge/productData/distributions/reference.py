# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Module for the 'reference' class for distributions."""

from xData import enums as xDataEnumsModule
from xData import link as linkModule

from . import base as baseModule

class Form( linkModule.Link, baseModule.Form ) :

    moniker = 'reference'
    subformAttributes = []

    def __init__( self, link = None, root = None, path = None, label = None, relative = False ) :

        linkModule.Link.__init__( self, link = link, root = root, path = path, label = label, relative = relative )
        baseModule.Form.__init__( self, label, xDataEnumsModule.Frame.none, [] )

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

    def energySpectrumAtEnergy(self, energyIn, frame, **kwargs):
        """Returns the energy spectrum in the lab frame for the specified incident energy."""

        return(self.referenceInstance.energySpectrumAtEnergy(energyIn, frame, **kwargs))

    def fixDomains(self, labels, energyMin, energyMax):
        """This method does nothing."""

        return 0

    def processMC_cdf( self, style, tempInfo, indent ) :

        # temporary solution:
        return( self.referenceInstance.processMC_cdf( style, tempInfo, indent ) )

        # better solution: add another link pointing to the processed version of what this points to:
        """
        newReference = Form( label=style.label, relative=True )
        self.ancestor.add( newReference )

        tempInfo['brokenLinks'].append( [self, newReference] )
        """

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.referenceInstance.processMultiGroup( style, tempInfo, indent ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        return( self.referenceInstance.toPointwise_withLinearXYs( **kwargs ) )

class CoulombPlusNuclearElastic(Form):

    moniker = 'CoulombPlusNuclearElastic'

class ThermalNeutronScatteringLaw( Form ) :

    moniker = 'thermalNeutronScatteringLaw'
