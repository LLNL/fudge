# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""This module contains the 'unspecified' form for distribution."""

from xData import enums as xDataEnumsModule
from . import base as baseModule

class Form( baseModule.Form ) :

    moniker = 'unspecified'
    subformAttributes = []

    def __init__(self, label, productFrame=xDataEnumsModule.Frame.lab):

        if productFrame is None:
            productFrame = xDataEnumsModule.Frame.lab      # This is a kludge and needs to be fixed. Ergo, None should never be allowed as a productFrame value.
        baseModule.Form.__init__( self, label, productFrame, [] )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )
        UC = cls( element.get( 'label' ), element.get( 'productFrame' ) )
        xPath.pop( )
        return( UC )

    @property
    def domainUnit( self ) :

        return( self.rootAncestor.domainUnit )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        return( [ None, None ] )

    def copy( self ):

        return Form( self.label, self.productFrame )

    def fixDomains(self, energyMin, energyMax, domainsToFix):
        """
        The method does nothing.
        """

        return 0

    def isSpecified( self ) :

        return( False )

    def processMC_cdf( self, style, tempInfo, indent ) :

        return( None )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( None )
