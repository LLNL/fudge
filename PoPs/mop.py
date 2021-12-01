# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__doc__ = """
This module contains classes to handle a MoP (Map of PoPs) instance.
"""

import os
import re

from xData import ancestry as ancestryModule
from . import database as databaseModule

class FormatVersion :

    version_0_1 = '0.1'

    allowed = ( version_0_1, )

class Base( ancestryModule.ancestry ) :
    """Base class used by all other mop related classes. Mainly inherits class ancestryModule.ancestry and defines toXML."""

    def __init( self ) :
        """
        Constructor for Base class.
        """

        ancestryModule.ancestry.__init__( self )

    def toXML( self, indent = '', **kwargs ) :
        """
        Returns an XML string of self.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param **kwargs:        A keyword list passed to self's toXML_list method.
        """

        return( '\n'.join( self.toXML_list( indent, **kwargs ) ) )

    def saveToOpenedFile( self, fOut, **kwargs ) :
        """
        Saves the contents of self to an already opened python file object (e.g., one returned by the open function).

        :param fOut:            A python file object.
        :param **kwargs:        A keyword list passed to self's toXML method.
        """

        xmlString = self.toXML( **kwargs )
        fOut.write( xmlString )
        fOut.write( '\n' )

    def saveToFile( self, fileName, **kwargs ) :
        """
        Saves the contents of self to the specified path (i.e., fileName).

        :param fileName:        The path to save the contents of self to.
        :param **kwargs:        A keyword list passed to self's saveToOpenedFile method.
        """

        dirname = os.path.dirname( fileName )
        if( ( len( dirname ) > 0 ) and not( os.path.exists( dirname ) ) ) : os.makedirs( dirname )
        with open( fileName, "w" ) as fout :
            fout.write( '<?xml version="1.0" encoding="UTF-8"?>\n' )
            self.saveToOpenedFile( fout, **kwargs )

class Mop( Base ) :

    moniker = 'mop'

    def __init__( self, name, path = './', format = FormatVersion.version_0_1 ) :

        """
        Constructor for a Mop instance.

        :param name:            The name of the library.
        :param path:            The path to the mop file if read in from disk.
        :param format:          The format version for the mop.
        """

        Base.__init__( self )

        if( not( isinstance( name, str ) ) ) : raise TypeError( 'Library must be a string.' )
        self.__name = name 

        if( not( isinstance( path, str ) ) ) : raise TypeError( 'Library must be a string.' )
        self.__path = path

        if( format not in FormatVersion.allowed ) : raise ValueError( 'Invalid format "%s."' % format )
        self.__format = format

        self.__entries = []
        self.__pops = databaseModule.database( 'Internal', '0.1' )

    def __len__( self ) :
        """Returns the number of entries in self."""

        return( len( self__entries ) )

    def __getitem__( self, index ) :
        """
        Returns the (index-1)^th item of self.

        :param index:       The index of the iter to return.
        """

        return( self.__entries[index] )

    @property
    def format( self ) :
        """Returns the format self is stored as."""

        return( self.__format )

    @property
    def name( self ) :
        """Returns self's name."""

        return( self.__name )

    @property
    def path( self ) :
        """Returns the path to self. If the data in self was read from disk, this is the path to the file."""

        return( self.__path )

    @property
    def fileName( self ) :
        """Same the path method."""

        return( self.__path )

    def append( self, entry ) :
        """Appends *entry* to the end of self's list of entries."""

        if( not( isinstance( entry, EntryBase ) ) ) : raise TypeError( 'Invalid instance to add.' )

        self.__entries.append( entry )
        entry.setAncestor( self )

    def fetch( self, re_pid, readPoPs ) :
        """
        Returns a list of items. Each item contains information about a PoPs database where one or more matches to *re_pid* exists. 
        Each item contains the list [ ids, all_ids, fileName, pops ]) where:

            -) ids:         The list of PoPs ids in a PoPs database that matches *re_pid*,
            -) all_ids:     The list of all PoPs ids in the PoPs database,
            -) fileName:    The file path to the PoPs database,
            -) pops:        An instance of the PoPs database if *readPoPs* is *True* and *None* otherwise.

        The *ids* are a list of str ids (e.g., [ 'Pu234', 'Pu234e_1' ]).  *re_pid* can be a str instance or a *re.compile* instance.
        Note, even if *readPoPs* is False, the returned pops will not be None if it had previously been read in.
        
        :param re_pid:          A str or a *re.compile* pattern representing the ids to match.
        :param readPoPs:        If *True*, the pops file is read and returned as the fourth element of each returned item.
        """

        matches = []
        for entry in self.__entries : matches += entry.fetch( re_pid, readPoPs )

        return( matches )

    def toXML_list( self, indent = "", **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of self.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param **kwargs:        A keyword list.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        format = kwargs.get( 'format', FormatVersion.version_0_1 )
        if( format not in FormatVersion.allowed ) : raise TypeError( 'Invalid format = "%s".' % format )

        XML_list = [ '%s<%s name="%s" format="%s">' % ( indent, self.moniker, self.name, format ) ]
        for entry in self.__entries : XML_list += entry.toXML_list( indent2, **kwargs )
        XML_list[-1] +=  '</%s>' % self.moniker

        return( XML_list )

    @staticmethod
    def readXML( path ) :
        """
        Reads in a XML mop file.

        :param path:        Path to a mop file.

        :return:            Mop instance containing all entries from the mop file.
        """

        from xml.etree import cElementTree
        from LUPY.xmlNode import xmlNode

        node = cElementTree.parse( path ).getroot( )
        node = xmlNode( node, xmlNode.etree )

        return( Mop.parseXML_node( node, path ) )

    @staticmethod
    def parseXML_node( node, path ) :
        """
        Creates a Mop instance from an XML mop element.

        :param node:        XML element to parse.

        :return:            Mop instance.
        """

        if( node.tag != Mop.moniker ) : raise TypeError( 'Invalid node name. Expected "%s" got "%s".' % ( Mop.moniker, node.tag ) )

        format = node.get( 'format', None )
        if( format is None ) : raise ValueError( 'Mop node does not have a "format" attribute.' )

        name = node.get( 'name', None )
        if( name is None ) : raise ValueError( 'Mop node does not have "name" attribute.' )

        mop = Mop( name, path = path )

        settings = { 'format' : format }
        for child in node :
            if( child.tag == Import.moniker ) :
                mop.append( Import.parseXML_node( settings, child ) )
            elif( child.tag == PoPs.moniker ) :
                mop.append( PoPs.parseXML_node( settings, child ) )
            else :
                raise ValueError( 'Invalid child tag "%s" for mop file.' % child.tag )

        return( mop )

class EntryBase( Base ) :
    """
    Base class for mop entries.
    """

    def __init__( self, path ) :
        """
        Base constructor for mop entries.

        :param path:            Path to the file for the entry.
        """

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( path, str ) ) ) : raise TypeError( 'Path must be a string.' )
        self.__path = path

    @property
    def fileName( self ) :
        """Returns the absolute path to the file pointed to by the 'path' attribute."""

        filename = self.path
        if( filename != os.sep ) : filename = os.path.join( os.path.dirname( os.path.realpath( self.ancestor.fileName ) ), filename )
        return( filename )

    @property
    def path( self ) :
        """Returns the path of the mop file. This may be absolute or relative."""

        return( self.__path )

class PoPs( EntryBase ) :
    """
    Class for storing a GNDS mop's PoPs instance.
    """

    moniker = 'pops'

    def __init__( self, path, ids ) :
        """
        Constructor for a mop's PoPs instance.

        :param path:            Path to the GNDS PoPs.database file.
        :param ids:             The list of particles ids in the PoPs.database file.
        """

        EntryBase.__init__( self, path )
        self.__pops = None
        self.__ids = ids

    @property
    def ids( self ) :
        """Returns the list of ids."""

        return( self.__ids )

    @property
    def pops( self ) :
        """Returns a reference to the PoPs.database instance. If the PoPs file has not be read in (see method *readIfNeeded*), None is returned."""

        return( self.__pops )

    def fetch( self, re_pid, readPoPs ) :
        """
        Returns the list [ ids, all_ids, fileName, pops ] if matches to *re_pid* are found in *self* and an empty list otherwise. 

            -) ids:         The list of PoPs ids in the PoPs database that matches *re_pid*,
            -) all_ids:     The list of all PoPs ids in the PoPs database,
            -) fileName:    The path to the PoPs database,
            -) pops:        An instance of the PoPs database if *readPoPs* is *True* and *None* otherwise.

        The *ids* are a list of str ids (e.g., [ 'Pu234', 'Pu234e_1' ]). *re_pid* can be a str instance or a *re.compile* instance.
        
        :param re_pid:          A str or *re.compile* pattern representing ids to match.
        :param readPoPs:        If *True*, the pops file is read and returned as the fourth element.
        """

        if( isinstance( re_pid, str ) ) : re_pid = re.compile( re_pid )
        matches = [ pid for pid in self.ids if re_pid.fullmatch( pid ) ]

        if( len( matches ) > 0 ) :
            matches = [ [ matches, self.ids, self.fileName, self.readIfNeeded( readPoPs ) ] ]

        return( matches )

    def readIfNeeded( self, readPoPs = True ) :
        """
        Reads in the PoPs.database referenced by the *path* member if it has not been read in by a prior call and *readPoPs* 
        is True (default). Also, returns the value of *self.__pops*. The returned value will only be *None* if *readPoPs* 
        is *False* and all prior calls to this method had "*readPoPs* = False".

        :param readPoPs:        If *True*, the pops file is read.
        :return:                The value if *self.__pops*.
        """

        if( ( self.__pops is None ) and readPoPs ) :
            self.__pops = databaseModule.database.readFile( self.fileName )
            self.__pops.setAncestor( self )

        return( self.__pops )

    def toXML_list( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of self.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param **kwargs:        A keyword list.

        :return:                List of str instances representing the XML lines of self.
        """

        return( [ '%s<%s path="%s">' % ( indent, self.moniker, self.path ) + ' '.join( self.__ids ) + '</%s>' % self.moniker ] )

    @staticmethod
    def parseXML_node( settings, node ) :
        """
        Creates an Import instance from an XML import element.

        :param settings:        Information from parent map instances needed to parse the XML element.
        :param node:            XML element to parse.

        :return:                Import instance.
        """

        if( node.tag != PoPs.moniker ) : raise TypeError( 'Invalid node name. Expected "%s" got "%s".' % ( PoPs.moniker, node.tag ) )

        path = node.get( 'path', None )
        ids = node.text.split( )

        pops = PoPs( path, ids )

        return( pops )

class Import( EntryBase ) :
    """
    Class for storing a GNDS mop's Import instance.
    """

    moniker = 'import'

    def __init__( self, path ) :
        """
        Constructor for a mop's import instance.

        :param path:            Path to the GNDS mop file.
        """

        EntryBase.__init__( self, path )
        self.__mop = None

    @property
    def mop( self ) :
        """Returns a reference to self.__mop."""

        return( self.__mop )

    def fetch( self, re_pid, readPoPs ) :
        """
        See documentation for Mop.fetch.
        
        :param re_pid:          A str or *re* pattern representing ids to match.
        :param readPoPs:        If *True*, the pops file is read and returned as the fourth element of each returned item.
        """

        self.readIfNeeded( )
        return( self.__mop.fetch( re_pid, readPoPs ) )

    def readIfNeeded( self ) :
        """
        Reads in the PoPs.database referenced by the *path* member if it has not been read in.
        """

        if( self.__mop is None ) :
            self.__mop = Mop.readXML( self.fileName )
            self.__mop.setAncestor( self )

    def toXML_list( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of self.

        :param indent:          The amount of indentation for each line. Child nodes and text may be indented more.
        :param **kwargs:        A keyword list.

        :return:                List of str instances representing the XML lines of self.
        """

        return( [ '%s<%s path="%s"/>' % ( indent, self.moniker, self.path ) ] )

    @staticmethod
    def parseXML_node( settings, node ) :
        """
        Creates an Import instance from an XML import element.

        :param settings:        Information from parent map instances needed to parse the XML element.
        :param node:            XML element to parse.

        :return:                Import instance.
        """

        if( node.tag != Import.moniker ) : raise TypeError( 'Invalid node name. Expected "%s" got "%s".' % ( Import.moniker, node.tag ) )

        path = node.get( 'path', None )

        _import = Import( path )

        return( _import )
