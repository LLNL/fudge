# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""This module contains the 'unspecified' form for distribution."""

from xData import standards as standardsModule

from . import base as baseModule

__metaclass__ = type

class form( baseModule.form ) :

    moniker = 'unspecified'
    subformAttributes = []
    ancestryMembers = ( '', )

    def __init__( self, label, productFrame = standardsModule.frames.labToken ) :

        baseModule.form.__init__( self, label, productFrame, [] )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )
        UC = form( element.get( 'label' ), element.get( 'productFrame' ) )
        xPath.pop( )
        return( UC )

    @property
    def domainUnit( self ) :

        return( self.getRootAncestor( ).domainUnit )

    def calculateAverageProductData( self, style, indent = '', **kwargs ) :

        return( [ None, None ] )

    def copy( self ):

        return form( self.label, self.productFrame )

    def isSpecified( self ) :

        return( False )

    def processMC_cdf( self, style, tempInfo, indent ) :

        return( None )

    def processMultiGroup( self, style, tempInfo, indent ) :

        return( None )
