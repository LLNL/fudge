# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from fudge import enums as enumsModule
from fudge import suites as suitesModule

from . import group as groupModule

class Transportable(ancestryModule.AncestryIO):
    """
    This class stores the product conserve and the group for one particle.
    """

    moniker = 'transportable'

    def __init__(self, particle, conserve, group):

        ancestryModule.AncestryIO.__init__(self)

        if( not( isinstance( particle, str ) ) ) : raise TypeError( 'particle must only be a string instance.' )
        self.__particle = particle

        self.__conserve = enumsModule.Conserve.checkEnumOrString(conserve)

        if( not( isinstance( group, groupModule.Group ) ) ) : raise TypeError( 'Group boundaries must only be a grid instance.' )
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

    def convertUnits(self, unitMap):
        '''
        Converts unit per *unitMap*.
        '''

        self.__group.convertUnits(unitMap)

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s label="%s" conserve="%s">' % ( indent, self.moniker, self.label, self.conserve ) ]
        xmlStringList += self.group.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append('%s[@label="%s"]' % (node.tag, node.get('label')))

        label = node.get('label')
        group = groupModule.Group.parseNodeUsingClass(node.find(groupModule.Group.moniker), xPath, linkData, **kwargs)
        conserve = node.get('conserve')

        xPath.pop()

        return cls( label, conserve, group )

class Transportables( suitesModule.Suite ) :
    """
    This class stores a Transportable instanse for each particle type.
    """

    moniker = 'transportables'

    def __init__( self ) :

        suitesModule.Suite.__init__( self, ( Transportable, ), allow_href = True )
