# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import ancestry as ancestryModule
from xData import standards as standardsModule
from xData import physicalQuantity as physicalQuantityModule

from .. import IDs as IDsPoPsModule
from .. import suite as suiteModule

class product( ancestryModule.ancestry ) :

    moniker = 'product'

    def __init__( self, pid, productFrame ) :

        self.__pid = pid
        self.__productFrame = productFrame

        self.__multiplicity = multiplicity.suite( )
        self.__multiplicity.setAncestor( self )

        self.__spectrum = spectrumModule.suite( )
        self.__spectrum.setAncestor( self )

    @property
    def pid( self ) :

        return( self.__pid )

    @property
    def productFrame( self ) :

        return( self.__productFrame )

    @property
    def multiplicity( self ) :

        return( self.__multiplicity )

    @property
    def spectrum( self ) :

        return( self.__spectrum )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s pid="%s" productFrame="%s">' % ( indent, self.moniker, self.__pid, self.__productFrame ) ]
        XMLStringList += self.__multiplicity.toXMLList( indent2, **kwargs )
        XMLStringList += self.__spectrum.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )
