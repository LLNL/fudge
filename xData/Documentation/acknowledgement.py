# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the GNDS documentation child node endfCompatible class.
"""

from .. import suite as suiteModule
from .. import text as textModule
from . import abstractClasses as abstractClassesModule

class Acknowledgement( textModule.Text ) :
    """A class representing a GNDS authors/acknowledgement node."""

    moniker = 'acknowledgement'

    def __init__( self, _label ) :

        textModule.Text.__init__( self )

        self.__label = textModule.raiseIfNotString( _label, 'label' )

    @property
    def label( self ) :

        return( self.__label )

    def XML_extraAttributes( self, **kwargs ) :

        if( self.__label == '' ) : return ''

        return ' label="%s"' % self.__label

    @staticmethod
    def parseConstructBareNodeInstance( node, xPath, linkData, **kwargs ) :

        label = node.get( 'label', '' )

        return( Acknowledgement( label ) )

class Acknowledgements( suiteModule.Suite ) :

    moniker = 'acknowledgements'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ Acknowledgement ] )
