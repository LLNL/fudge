# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the Affiliation, Affiliations, Note and AuthorAbstract classes.
"""

from LUPY import ancestry as ancestryModule

from .. import suite as suiteModule
from .. import text as textModule

class Affiliation( textModule.Text ) :
    """A class representing a GNDS author/affiliations/affiliation node."""

    moniker = 'affiliation'
    keyName = 'label'

    def __init__( self, name = '', href = '' ) :

        textModule.Text.__init__( self )

        self.name = name
        self.href = href

    @property
    def name( self ) :
        """Returns self's name instance."""

        return( self.__name )

    @name.setter
    def name( self, value ) :

        self.__name = textModule.raiseIfNotString( value, 'name' )

    @property
    def href( self ) :
        """Returns self's href instance."""

        return( self.__href )

    @href.setter
    def href( self, value ) :

        self.__href = textModule.raiseIfNotString( value, 'href' )

    def XML_extraAttributes( self, **kwargs ) :

        attributes = ''
        if( len( self.__name ) > 0 ) : attributes += ' name="%s"' % self.__name
        if( len( self.__href ) > 0 ) : attributes += ' href="%s"' % self.__href

        return attributes

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        name = node.get( 'name' )
        href = node.get( 'href' )

        return cls(name, href)

class Affiliations( suiteModule.Suite ) :

    moniker = 'affiliations'
    suiteName = 'label'

    def __init__( self ) :

        suiteModule.Suite.__init__( self, [ Affiliation ] )

class Note( textModule.Text ) :
    """A class representing a GNDS authors/author/note node."""

    moniker = 'note'

class AuthorAbstract(ancestryModule.AncestryIO):

    ancestryMembers = ( 'affiliations', 'note' )

    def __init__( self, name, orcid, email ) :

        ancestryModule.AncestryIO.__init__(self)

        self.__name = textModule.raiseIfNotString( name, 'name' )
        self.__orcid = textModule.raiseIfNotString( orcid, 'orcid' )
        self.__email = textModule.raiseIfNotString( email, 'email' )

        self.__affiliations = Affiliations( )
        self.__affiliations.setAncestor( self )

        self.__note = Note( )
        self.__note.setAncestor( self )

    @property
    def name( self ) :
        """."""

        return( self.__name )

    @property
    def orcid(self):

        return self.__orcid

    @property
    def email( self ) :
        """."""

        return( self.__email )

    @property
    def affiliations( self ) :
        """."""

        return( self.__affiliations )

    @property
    def note( self ) :
        """."""

        return( self.__note )

    def XML_extraAttributes( self, **kwargs ) :

        attributes = ' name="%s"' % self.__name
        if len(self.__orcid) > 0: attributes += ' orcid="%s"' % self.__orcid
        if len(self.__email) > 0: attributes += ' email="%s"' % self.__email

        return attributes

    def toXML_strList(self, indent = '', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XMLList = [ '%s<%s%s>' % ( indent, self.moniker, self.XML_extraAttributes(**kwargs) ) ]
        XMLList += self.__affiliations.toXML_strList( indent = indent2, **kwargs )
        XMLList += self.__note.toXML_strList( indent = indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker

        return( XMLList )
