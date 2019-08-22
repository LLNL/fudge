# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

__metaclass__ = type

import string
import abc
import inspect

from xData import ancestry as ancestryModule

class suite( ancestryModule.ancestry ) :
    """
    Abstract base class for holding a list of objects that can be indexed by an integer or a string. 
    For integer indexing, the indexing is by order in which an object was put into the array or, 
    if a sortedSuite instance, by its sort order in the suite.
    For string indexing, this class acts like the python dict class. This class
    is like the python OrderedDict class. However, this class does not inherit the OrderedDict class
    as it does act differently. For example, only objects in the list of allowedClasses can be
    added to an suite instance. The allowedClasses is an argument to the constructor.

    The constructor for this class has 3 (not counting self) arguments. They are::

        + ------------------+-------------------------------------------------------+
        | allowedClasses    | A list of classes. Only objects matching one of the   |
        |                   | classes can be add.                                   |
        + ------------------+-------------------------------------------------------+
        | key               | The name of the member in an object that is used as   |
        |                   | the key. Default is 'label'.                          |
        + ------------------+-------------------------------------------------------+
        | replace           | If False, the add method will only allow an object to |
        |                   | be insert if no object object in the list has the     |
        |                   | same key. If True, and an object of the same key      |
        |                   | exists, the added object replaces the old object at   |
        |                   | position in the list. Default is True.                |
        + ------------------+-------------------------------------------------------+
    """

    __metaclass__ = abc.ABCMeta

    def __init__( self, allowedClasses, key = 'label', replace = True ) :

        ancestryModule.ancestry.__init__( self )

        __allowedClasses = []
        for i1, cls in enumerate( allowedClasses ) :
            if( not( inspect.isclass( cls ) ) ) : raise TypeError( 'Item at index %d is not a class' % i1 )
            __allowedClasses.append( cls )
        self.__allowedClasses = tuple( __allowedClasses )   # Make it a tuple so that it cannot be changed.

        if( not( isinstance( replace, bool ) ) ) : raise TypeError( 'replace must be a bool instance.' )
        self.__replace = replace

        if( not( isinstance( key, str ) ) ) : raise TypeError( 'key must be a string instance.' )
        self.__key = key

        self.__items = []

    def __contains__( self, key ) :

        for item in self.__items :
            if( item.key == key ) : return( True )
        return( False )

    def __getitem__( self, key ) :
        if not self.__items: raise KeyError('suite is empty')
        if( isinstance( key, int ) ) : return( self.__items[key] )
        if( not( isinstance( key, str ) ) ) : raise TypeError( 'key must be a string' )

        for item in self.__items :
            if( item.key == key ) : return( item )
        raise KeyError( 'item with key "%s" not found in suite "%s"' % ( key, self.moniker ) )

    def __iter__( self ) :

        n1 = len( self )
        for i1 in range( n1 ) : yield self.__items[i1]

    def __len__( self ) :

        return( len( self.__items ) )

    @property
    def allowedClasses( self ) :

        return( self.__allowedClasses )     # No fear, user cannot alter as __allowedClasses is a tuple.

    @property
    def key( self ) :

        return( self.__key )

    @property
    def replace( self ) :

        return( self.__replace )

    def add( self, item ) :

        classAllowed = False
        for cls in self.__allowedClasses :
            if( isinstance( item, cls ) ) :
                classAllowed = True
                break
        if( not( classAllowed ) ) : raise TypeError( 'Invalid class "%s" for suite "%s"' % ( item.__class__, self.moniker ) )

        index, replace = self.addIndex( item )
        item.setAncestor( self, attribute = self.__key )
        if( replace ) : del self.__items[index]
        self.__items.insert( index, item )

    def addIndex( self, item ) :
        """
        Returns the index at which item will be added to self.
        """

        index = len( self )
        for i1, _item in enumerate( self.__items ) :
            if( _item.key == item.key ) :
                if( self.replace ) : return( i1, True )
                raise KeyError( 'item with key = "%s" already present in non-replaceable suite' % getattr( item.key ) )
        return( index, False )

    def check( self, info ):

        from . import warning as warningModule
        warnings = []
        for child in self:
            childWarnings = child.check( info )
            if childWarnings:
                warnings.append( warningModule.context( '%s: %s' % (child.moniker, str(child)), childWarnings ) )
        return warnings

    def convertUnits( self, unitMap ) :

        for item in self.__items : item.convertUnits( unitMap )

    def remove( self, key ) :

        for i1, item in enumerate( self.__items ) :
            if( item.key == key ) :
                del self.__items[i1]
                return( True )
        return( False )

    def toXML( self, indent = "", **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        if( len( self ) == 0 ) : return( [] )
        xmlString = [ '%s<%s>' % ( indent, self.moniker ) ]
        for item in self : xmlString += item.toXMLList( indent2, **kwargs )
        xmlString[-1] += '</%s>' % self.moniker
        return( xmlString )

    def uniqueLabel( self, item ) :
        """
        If item's key is the same as another item's key in self, construct a new unique key
        based on item's key appended with '__' and one or more lower case letters (i.e., 'a' to 'z').
        """

        def incrementSuffix( suffix ) :

            if( len( suffix ) == 0 ) : return( 'a' )
            index = string.ascii_lowercase.find( suffix[-1] ) + 1
            if( index != 26 ) : return( suffix[:-1] + string.ascii_lowercase[index] )
            return( incrementSuffix( suffix[:-1] ) + 'a' )

        if( item.key in self ) :
            key__ = item.key + '__'
            n1 = len( key__ )
            l1 = 0
            suffixes = []
            for _item in self :          # Find list of longest keys that start with key__.
                subkey = _item.key[:n1]
                if( subkey == key__ ) :
                    suffix = _item.key[n1:]
                    if( not( suffix.islower( ) ) ) : continue       # Ignore non-standard keys.
                    l2 = len( suffix )
                    if( l2 < l1 ) : continue
                    if( l2 > l1 ) :
                        l1 = l2
                        suffixes = []
                    suffixes.append( suffix )
            if( len( suffixes ) == 0 ) :
                suffix = 'a'
            else :
                suffix = incrementSuffix( sorted( suffixes )[-1] )
            item.key = key__ + suffix

        return( item )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        for child in element :
            parseClass = None
            for _class in self.__allowedClasses :
                if( child.tag == _class.moniker ) :
                    parseClass = _class
                    break
            if( parseClass is None ) :
                raise TypeError( 'Invalid element "%s" encountered in suite "%s"' % ( child.tag, self.moniker ) )

            self.add( parseClass.parseXMLNodeAsClass( child, xPath, linkData ) )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        kwargs = element.attrib.copy( )
        if( 'Z' in kwargs ) : kwargs['Z'] = int( kwargs['Z'] )
        return( cls( **kwargs ).parseXMLNode( element, xPath, linkData ) )

    @classmethod
    def parseXMLStringAsClass( cls, string ) :

        from xml.etree import cElementTree

        element = cElementTree.fromstring( string )
        kwargs = element.attrib.copy( )
        return( cls( **kwargs ).parseXMLNode( element, [], [] ) )

class sortedSuite( suite ) :

    def addIndex( self, item ) :

        for i1, _item in enumerate( self ) :
            cmp = item.sortCompare( _item )
            if( cmp == 0 ) :
                if( self.replace ) : return( i1, True )
                raise KeyError( 'item with key = "%s" already present in non-replaceable suite' % getattr( item, self.__key ) )
            elif( cmp < 0 ) :
                return( i1, False )
        return( len( self ), False )
