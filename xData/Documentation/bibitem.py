# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the GNDS documentation child node endfCompatible class.
"""

from .. import text as textModule
from .. import suite as suiteModule

class Bibitem( textModule.Text ) :
    """A class representing a GNDS authors/bibitems/bibitem node."""

    moniker = 'bibitem'

    def __init__( self, xref, text = None ) :

        textModule.Text.__init__( self, text=text)

        self.__xref = textModule.raiseIfNotString( xref, 'xref' )

    @property
    def xref( self ) :

        return ( self.__xref )

    def XML_extraAttributes( self, **kwargs ) :

        if( self.__xref == '' ) : return ''

        return ' xref="%s"' % self.__xref

    @staticmethod
    def parseConstructBareNodeInstance( node, xPath, linkData, **kwargs ) :

        _xref= node.get( 'xref' )

        return( Bibitem( _xref ) )

class Bibliography( suiteModule.Suite ) :
    """A class representing a GNDS authors/bibitems node."""

    moniker = 'bibliography'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ Bibitem ] )
