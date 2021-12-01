# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Coherent and incoherent photon scattering forms.
"""

from xData import link as linkModule

from . import base as baseModule

class baseForm( baseModule.form, linkModule.link ) :

    def __init__( self, link = None, root = None, path = None, relative = False, label = None ) :

        linkModule.link.__init__( self, link = link, root = root, path = path, relative = relative, label = label )
        baseModule.form.__init__( self, label, None, [] )

    @property
    def productFrame( self ) :

        return( self.link.productFrame )

    def check( self, info ) :

        return( [] )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        return( self.link.calculateAverageProductData( style, indent, **kwargs ) )

    def processMC_cdf( self, style, tempInfo, indent ) :

        return( self.link.processMC_cdf( style, tempInfo, indent ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.link.processMultiGroup( style, tempInfo, indent ) )

    def toXMLList( self, indent = '', **kwargs ) :

        return( linkModule.link.toXMLList( self, indent = indent, **kwargs ) )

class coherentPhotonScattering :

    class form( baseForm ) :

        moniker = 'coherentPhotonScattering'

class incoherentPhotonScattering :

    class form( baseForm ) :

        moniker = 'incoherentPhotonScattering'
