# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

# FIXME CMM this appears to be unused.  Likely replaced by the fissionFramentData.rate module
# Remove?

from LUPY import ancestry as ancestryModule

from .. import suite as suiteModule

class DecayConstant(ancestryModule.AncestryIO):

    moniker = 'decayConstant'

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__(self)

        self.__rate = rateModule.Suite( )
        self.__rate.setAncestor( self )

        self.__product = productModule.product( )       # FIXME; productModule does not exist.
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

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.rate.toXML_strList( indent = indent2, **kwargs )
        XMLStringList += self.product.toXML_strList( indent = indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

class Suite( suiteModule.SortedSuite ) :

    moniker = 'decayConstants'

    def __init__( self ) :

        suiteModule.SortedSuite.__init__( self, allowedClasses = ( DecayConstant, ) )
