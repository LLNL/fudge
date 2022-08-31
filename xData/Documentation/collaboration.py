# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from .. import suite as suiteModule

"""
This module contains the class for the GNDS documentation child node collaboration.
"""

from LUPY import ancestry as ancestryModule

class Collaboration(ancestryModule.AncestryIO):
    """A class representing a GNDS documentation/collaborations/collaboration node."""

    moniker = 'collaboration'
    keyName = 'label'

    def __init__(self, name, href=''):

        ancestryModule.AncestryIO.__init__(self)

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

    def toXML_strList( self, indent = '', **kwargs ) :

        href = ''
        if( len( self.__href ) > 0 ) : href = ' href="%s"' % self.__href

        return( [ '%s<%s name="%s"%s/>' % ( indent, self.moniker, self.__name, href ) ] )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        name = node.get( 'name' )
        href = node.get( 'href' )

        return cls( name, href)

class Collaborations( suiteModule.Suite ) :

    moniker = 'collaborations'
    suiteName = 'label'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ Collaboration ] )
