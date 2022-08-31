# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from .. import IDs as IDsPoPsModule

class Product(ancestryModule.AncestryIO):

    moniker = 'product'

    def __init__( self, pid, productFrame ) :

        ancestryModule.AncestryIO.__init__(self)

        self.__pid = pid
        self.__productFrame = productFrame

        self.__multiplicity = multiplicity.Suite( )
        self.__multiplicity.setAncestor( self )

        self.__spectrum = spectrumModule.Suite( )
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

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s pid="%s" productFrame="%s">' % ( indent, self.moniker, self.__pid, self.__productFrame ) ]
        XMLStringList += self.__multiplicity.toXML_strList( indent2, **kwargs )
        XMLStringList += self.__spectrum.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )
