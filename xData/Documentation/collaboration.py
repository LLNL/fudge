# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from .. import suite as suiteModule

"""
This module contains the class for the GNDS documentation child node collaboration.
"""

from .. import ancestry as ancestryModule

class Collaboration( ancestryModule.Ancestry2 ) :
    """A class representing a GNDS documentation/collaborations/collaboration node."""

    moniker = 'collaboration'

    def __init__( self, name, href = '' ) :

        self.__name = name
        self.__href = href

    @property
    def name( self ) :
        """Returns self's label instance."""

        return( self.__name )

    @property
    def href( self ) :
        """Returns self's href instance."""

        return( self.__href )

    def toXMLList( self, indent = '', **kwargs ) :

        href = ''
        if( len( self.__href ) > 0 ) : href = ' href="%s"' % self.__href

        return( [ '%s<%s name="%s"%s/>' % ( indent, self.moniker, self.__name, href ) ] )

    @staticmethod
    def parseConstructBareNodeInstance( node, xPath, linkData, **kwargs ) :

        name = node.get( 'name' )
        href = node.get( 'href' )

        return( Collaboration( name, href ) )

class Collaborations( suiteModule.Suite ) :

    moniker = 'collaborations'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ Collaboration ] )
