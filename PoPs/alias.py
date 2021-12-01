# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the alias classes.
"""

import abc

from . import misc as miscModule
from . import suite as suiteModule

from .groups import misc as chemicalElementMiscModule

class alias( miscModule.classWithIDKey ) :

    __metaclass__ = abc.ABCMeta

    def __init__( self, id, pid ) :

        miscModule.classWithIDKey.__init__( self, id )

        if( not( isinstance( pid, str ) ) ) : raise TypeError( 'pid not str' )
        self.__pid = pid

    @property
    def pid( self ) :

        return( self.__pid )

    def check( self, info ):
        from . import warning as warningModule
        warnings = []
        if self.pid not in info['PoPs']:
            warnings.append( warningModule.AliasToNonExistentParticle(self.id, self.pid, self))
        return warnings

    def copy( self ) :

        return( self.__class__( self.id, self.pid ) )

    def isAlias( self ) :

        return( True )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        return( [ '%s<%s id="%s" pid="%s"/>' % ( indent, self.moniker, self.id, self.pid ) ] )

class particle( alias ) :

    moniker = 'particle'

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        self = cls( element.get('id'), element.get('pid') )

        xPath.pop( )
        return( self )

class metaStable( alias ) :

    moniker = 'metaStable'

    def __init__( self, id, pid, metaStableIndex ) :

        alias.__init__( self, id, pid )

        if( not( isinstance( metaStableIndex, int ) ) ) :
            TypeError( 'metaStableIndex must be an int: %s' % miscModule.toLimitedString( metaStableIndex ) )
        self.__metaStableIndex = metaStableIndex

    @property
    def metaStableIndex( self ) :

        return( self.__metaStableIndex )

    def toXMLList( self, indent = '', **kwargs ) :

        return( [ '%s<%s id="%s" pid="%s" metaStableIndex="%s"/>' % ( indent, self.moniker, self.id, self.pid, self.metaStableIndex ) ] )

    def copy( self ) :

        return( self.__class__( self.id, self.pid, self.metaStableIndex ) )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        self = cls( element.get('id'), element.get('pid'), int( element.get('metaStableIndex') ) )

        xPath.pop( )
        return( self )

    @staticmethod
    def metaStableNameFromNuclearLevelNameAndMetaStableIndex( nuclideName, metaStableIndex ) :

        if( not( isinstance( metaStableIndex, int ) ) ) :
            TypeError( 'metaStableIndex must be an int: %s' % miscModule.toLimitedString( metaStableIndex ) )
        if( metaStableIndex < 1 ) : raise ValueError( 'metaStableIndex must be greater than 0 got "%s".' % metaStableIndex )

        baseID, chemicalElementID, A, levelID, isNucleus, anti, qualifier = chemicalElementMiscModule.chemicalElementALevelIDsAndAnti( nuclideName )

        if( isNucleus ) : chemicalElementID = chemicalElementID[0].lower( ) + chemicalElementID[1:]
        return( "%s_m%d" % ( chemicalElementMiscModule.isotopeSymbolFromChemicalElementIDAndA( chemicalElementID, A ), metaStableIndex ) )

    @staticmethod
    def nuclideNameAndMetaStableIndexFromName( name ) :
        """
        This function returns the nuclide name and meta-stable index from its argument name. If name does not appear to be a 
        meta-stable, ( name, 0 ) are returned. This function splits on the string '_m'. If the number of sub-strings returned 
        is not 2, the name is considered not to be a meta-stable and ( name, 0 ) are returned. For example, name = 'O16' will 
        return ( 'O16', 0 ), name = 'Am242_m1' will return ( 'Am242_m1', 1 ) and name = 'Am242_m1_m2' will return ( 'Am242_m1_m2', 0 ).
        """

        if( name.count( '_m' ) != 1 ) : return( name, 0 )
        nuclideName, metaStableIndex = name.split( '_m' )
        try :
            metaStableIndex = int( metaStableIndex )
        except :
            metaStableIndex = 0

        return( nuclideName, metaStableIndex )

class suite( suiteModule.suite ) :

    moniker = 'aliases'

    def __init__( self ) :

        suiteModule.suite.__init__( self, ( alias, ) )

    def has_pid( self, ParticleID ) :
        """Returns True if one of the aliases has pid equal to ParticleID and False otherwise."""

        for alias in self :
            if( alias.pid == ParticleID ) : return( True )

        return( False )

    def parseXMLNode( self, element, xPath, linkData ) :

        for child in element :
            if( child.tag == metaStable.moniker ) :
                self.add( metaStable.parseXMLNodeAsClass( child, xPath, linkData ) )
            else :
                self.add( particle.parseXMLNodeAsClass( child, xPath, linkData ) )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        self = cls( )
        self.parseXMLNode( element, xPath, linkData )

        xPath.pop( )
        return( self )
