# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Multi-group distribution classes.
"""
import abc

from xData import gridded as griddedModule

from . import base as baseModule

noExtrapolationToken = 'noExtrapolation'

__metaclass__ = type

class subform( baseModule.subform ) :
    """
    Abstract base class for multiGroup subforms.
    """
    __metaclass__ = abc.ABCMeta

    def __init__( self ) :

        baseModule.subform.__init__( self )

class gridded3d( subform, griddedModule.gridded3d ) :

    def __init__( self, **kwargs ) :

        subform.__init__( self )
        griddedModule.gridded3d.__init__( self, **kwargs )

    @staticmethod
    def dataToString( values, self, indent = '', **kwargs ) :

        valueFormatter = kwargs['valueFormatter']
        sep = kwargs.get( 'sep', values.sep )

        XMLStrList = [ ]
        start = 0
        for length in self.lengths :
            strList = [ valueFormatter( self.data[i1] ) for i1 in range( start, start + length ) ]
            XMLStrList.append( indent + sep.join( strList ) )
            start += length
        return( XMLStrList )

class form( baseModule.form ) :

    moniker = 'multiGroup3d'
    subformAttributes = ( 'multiGroupSubform', )

    def __init__( self, label, productFrame, multiGroupSubform ) :

        if( not( isinstance( multiGroupSubform, subform ) ) ) :
            raise TypeError( 'instance is not a multiGroup subform: moniker = %s' % multiGroupSubform.moniker )
        baseModule.form.__init__( self, label, productFrame, ( multiGroupSubform, ) )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        xPath.append( element.tag )

        subformElement = element[0]
        subformClass = {
                gridded3d.moniker : gridded3d,
            }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( "unknown gridded3d subform: %s" % subformElement.tag )
        SubForm = subformClass.parseXMLNode( subformElement, xPath, linkData )

        xPath.pop( )
        return form( element.get( 'label' ), element.get( 'productFrame' ), SubForm )
