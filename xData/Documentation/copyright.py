# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the class for the GNDS documentation child node copyright.
"""

from .. import text as textModule

class Copyright( textModule.Text ) :
    """A class representing a GNDS documentation/copyright node."""

    moniker = 'copyright'

    def __init__( self, href = '' ) :

        textModule.Text.__init__( self )

        self.href = href

    @property
    def href( self ) :
        """Returns self's href instance."""

        return( self.__href )

    @href.setter
    def href( self, value ) :

        if( not( isinstance( value, str ) ) ) : raise TypeError( 'href must be a str instance.' )
        self.__href = value

    def XML_extraAttributes( self, **kwargs ) :

        if( self.href == '' ) : return ''

        return ' href="%s"' % self.href

    def parseNode(self, node, xPath, linkData, **kwargs):
        """
        Parses a copyright node.
        """

        textModule.Text.parseNode(self, node, xPath, linkData, **kwargs)
        xPath.append( node.tag )

        self.__href = node.get( 'href', '' )

        xPath.pop( )
