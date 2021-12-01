# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import ancestry as ancestryModule

from .. import suite as suiteModule

class decayConstant( ancestryModule.ancestry ) :

    moniker = 'decayConstant'

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

        self.__rate = rateModule.suite( )
        self.__rate.setAncestor( self )

        self.__product = productModule.product( )
        self.__product.setAncestor( self )

    @property
    def rate( self ) :

        return( self.__rate )

    @property
    def product( self ) :

        return( self.__product )

    def convertUnits( self, unitMap ) :

        self.__rate.convertUnits( unitMap )
        self.__product__rate.convertUnits( unitMap )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.rate.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.product.toXMLList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

class suite( suiteModule.sortedSuite ) :

    moniker = 'decayConstants'

    def __init__( self ) :

        suiteModule.sortedSuite.__init__( self, allowedClasses = ( decayConstant, ) )
