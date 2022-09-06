# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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

class Subform( baseModule.Subform, abc.ABC ) :
    """
    Abstract base class for multiGroup subforms.
    """

    def __init__( self ) :

        baseModule.Subform.__init__( self )

class Gridded3d( Subform, griddedModule.Gridded3d ) :

    def __init__(self, axes, array, **kwargs):

        Subform.__init__(self)
        griddedModule.Gridded3d.__init__(self, axes, array, **kwargs)

    @staticmethod
    def dataToString( values, self, indent = '', **kwargs ) :

        valueFormatter = kwargs['valueFormatter']

        XMLStrList = [ ]
        start = 0
        for length in self.lengths :
            strList = [ valueFormatter( self.data[i1] ) for i1 in range( start, start + length ) ]
            XMLStrList.append( indent + ' '.join( strList ) )
            start += length
        return( XMLStrList )

class Form( baseModule.Form ) :

    moniker = 'multiGroup3d'
    subformAttributes = ( 'multiGroupSubform', )

    def __init__( self, label, productFrame, multiGroupSubform ) :

        if( not( isinstance( multiGroupSubform, Subform ) ) ) :
            raise TypeError( 'instance is not a multiGroup subform: moniker = %s' % multiGroupSubform.moniker )
        baseModule.Form.__init__( self, label, productFrame, ( multiGroupSubform, ) )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )

        subformElement = element[0]
        subformClass = {
                Gridded3d.moniker : Gridded3d,
            }.get( subformElement.tag )
        if( subformClass is None ) : raise Exception( "unknown Gridded3d subform: %s" % subformElement.tag )
        SubForm = subformClass.parseNodeUsingClass(subformElement, xPath, linkData, **kwargs)

        instance = cls( element.get( 'label' ), element.get( 'productFrame' ), SubForm )

        xPath.pop( )

        return  instance
