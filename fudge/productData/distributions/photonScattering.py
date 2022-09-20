# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Coherent and incoherent photon scattering forms.
"""

from xData import enums as xDataEnumsModule
from xData import link as linkModule

from . import base as baseModule

class BaseForm( baseModule.Form, linkModule.Link ) :

    def __init__( self, link = None, root = None, path = None, relative = False, label = None ) :

        linkModule.Link.__init__( self, link = link, root = root, path = path, relative = relative, label = label )
        baseModule.Form.__init__( self, label, xDataEnumsModule.Frame.none, [] )

    @property
    def productFrame( self ) :

        return( self.link.productFrame )

    def check( self, info ) :

        return( [] )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        return( self.link.calculateAverageProductData( style, indent, **kwargs ) )

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """This method does nothing."""

        return 0

    def processMC_cdf( self, style, tempInfo, indent ) :

        return( self.link.processMC_cdf( style, tempInfo, indent ) )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( self.link.processMultiGroup( style, tempInfo, indent ) )

    def toXML_strList( self, indent = '', **kwargs ) :

        return( linkModule.Link.toXML_strList( self, indent = indent, **kwargs ) )

class CoherentPhotonScattering :

    class Form(BaseForm):

        moniker = 'coherentPhotonScattering'

class IncoherentPhotonScattering :

    class Form(BaseForm):

        moniker = 'incoherentPhotonScattering'
