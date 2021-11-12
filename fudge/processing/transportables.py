# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__metaclass__ = type

from xData import ancestry as ancestryModule
from xData import axes as axesModule

from fudge import suites as suitesModule

from . import group as groupModule

class conserve :

    number = 'number'
    energyOut = 'energyOut'

class transportable( ancestryModule.ancestry ) :
    """
    This class stores the product conserve and the group for one particle.
    """

    moniker = 'transportable'

    def __init__( self, particle, _conserve, group ) :

        if( not( isinstance( particle, str ) ) ) : raise TypeError( 'particle must only be a string instance.' )
        if( _conserve not in [ conserve.number, conserve.energyOut ] ) :
            raise TypeError( 'Invalid product conserve' )
        if( not( isinstance( group, groupModule.group ) ) ) : raise TypeError( 'Group boundaries must only be a grid instance.' )

        self.__particle = particle
        self.__conserve = _conserve
        self.__group = group

    @property
    def label( self ) :

        return( self.particle )

    @property
    def particle( self ) :

        return( self.__particle )

    @property
    def conserve( self ) :

        return( self.__conserve )

    @property
    def group( self ) :

        return( self.__group )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s label="%s" conserve="%s">' % ( indent, self.moniker, self.label, self.conserve ) ]
        xmlStringList += self.group.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        label = element.get('label')
        xPath.append( '%s[@label="%s"]' % (element.tag,label) )
        group = groupModule.group.parseXMLNode( element.find(
            groupModule.group.moniker ), xPath, linkData )
        conserve = element.get('conserve')
        xPath.pop()
        return transportable( label, conserve, group )


class transportables( suitesModule.suite ) :
    """
    This class stores a transportable instanse for each particle type.
    """

    moniker = 'transportables'

    def __init__( self ) :

        suitesModule.suite.__init__( self, ( transportable, ), allow_href = True )
