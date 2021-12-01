# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the GNDS documentation child nodes keywords and keywork classes.
"""

from .. import ancestry as ancestryModule
from .. import suite as suiteModule
from .. import text as textModule

class Keyword( textModule.Text ) :
    """A class representing a GNDS documentation/abstract node."""

    moniker = 'keyword'

    def __init__( self, label, type, text ) :

        textModule.Text.__init__( self, text, label = label )

        self.__type = type

    @property
    def type( self ) :

        return( self.__type )

    def XML_extraAttributes( self, **kwargs ) :

        if( self.__type == '' ) : return ''

        return ' type="%s"' % self.__type

    @staticmethod
    def parseConstructBareNodeInstance( node, xPath, linkData, **kwargs ) :

        label = node.get( 'label' )
        type = node.get( 'type' )

        return( Keyword( label, type, None ) )

class Keywords( suiteModule.Suite ) :

    moniker = 'keywords'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ Keyword ] )
