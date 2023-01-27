# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>


import string
import abc
import inspect

from LUPY import ancestry as ancestryModule

class Suite(ancestryModule.AncestryIO_bare, abc.ABC):
    """
    Abstract base class for holding a list of objects that can be indexed by an integer or a string.
    For integer indexing, the indexing is by order in which an object was put into the array or,
    if a sortedSuite instance, by the index assigned when the object was added to the suite.
    For string indexing, this class acts like the python dictionary.

    The suite also has a tuple of 'allowedClasses' controlling what type(s) of objects it can contain.
    The allowedClasses must be passed in as an argument to the constructor.

    """

    def __init__( self, allowedClasses, replace = True ) :
        """
        :param allowedClasses: tuple of classes. Only instances of these classes can be added to the suite.
        :param replace: boolean (default=True). If False, the add method can only add objects with keys
            not already present in the suite. If True, the original object will be replaced by the new object
            having the same key.
        """

        ancestryModule.AncestryIO.__init__(self)

        __allowedClasses = []
        for i1, cls in enumerate( allowedClasses ) :
            if( not( inspect.isclass( cls ) ) ) : raise TypeError( 'Item at index %d is not a class' % i1 )
            __allowedClasses.append( cls )
        self.__allowedClasses = tuple( __allowedClasses )   # Make it a tuple so that it cannot be changed.

        if( not( isinstance( replace, bool ) ) ) : raise TypeError( 'replace must be a bool instance.' )
        self.__replace = replace

        self.__items = []

    def __contains__( self, key ) :

        for item in self.__items :
            if( item.key == key ) : return( True )
        return( False )

    def __getitem__( self, key ) :

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
        """Return the tuple of classes that are allowed in this suite."""

        return( self.__allowedClasses )     # No fear, user cannot alter as __allowedClasses is a tuple.

    @property
    def replace( self ) :

        return( self.__replace )

    def add( self, item ) :
        """
        Attempt to add item to self.
        Raises TypeError if item is not an instance of any class in self.allowedClasses
        Raises KeyError if self.replace == False and the item shares a key with another item already in the suite.

        :param item: instance to add to self
        """

        classAllowed = False
        for cls in self.__allowedClasses :
            if( isinstance( item, cls ) ) :
                classAllowed = True
                break
        if( not( classAllowed ) ) : raise TypeError( 'Invalid class "%s" for suite "%s"' % ( item.__class__, self.moniker ) )

        index, replace = self.addIndex( item )
        item.setAncestor(self)
        if( replace ) : del self.__items[index]
        self.__items.insert( index, item )

    def addIndex( self, item ) :
        """
        Returns the index at which item will be added to self.
        Raises KeyError if self.replace == False and the item shares a key with another item already in the suite.
        """

        index = len( self )
        for i1, _item in enumerate( self.__items ) :
            if( _item.key == item.key ) :
                if( self.replace ) : return( i1, True )
                raise KeyError( 'item with key = "%s" already present in non-replaceable suite' % item.key )
        return( index, False )

    def check( self, info ):
        """See documentation in database.check"""

        from . import warning as warningModule

        warnings = []
        for child in self:
            childWarnings = child.check( info )
            if childWarnings:
                warnings.append( warningModule.Context( '%s: %s' % (child.moniker, str(child)), childWarnings ) )
        return warnings

    def convertUnits( self, unitMap ) :
        """See documentation in database.convertUnits"""

        for item in self.__items : item.convertUnits( unitMap )

    def keys(self):
        '''Returns the list of all keys in *self*.'''

        return [item.keyValue for item in self.__items]

    def remove( self, key ) :
        """Remove object with specified key from the suite."""

        for i1, item in enumerate( self.__items ) :
            if( item.key == key ) :
                del self.__items[i1]
                return( True )
        return( False )

    def replicate( self, other ) :
        """
        Add copies of all items in other to self.
        :param other: suite or other iterable
        """

        for item in other : self.add( item.copy( ) )

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlString = ['%s<%s>' % (indent, self.moniker)]
        for item in self :
            xmlString += item.toXML_strList(indent2, **kwargs)

        if len(xmlString) == 1:
            if kwargs.get('showEmptySuites', False):
                xmlString = [xmlString[0][:-1] + '/>']
            else:
                xmlString = []
        else:
            xmlString[-1] += '</%s>' % self.moniker

        return xmlString

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

    def parseNode(self, node, xPath, linkData, **kwargs):

        xPath.append( node.tag )

        for child in node:
            parseClass = None
            for _class in self.__allowedClasses:
                if child.tag == _class.moniker:
                    parseClass = _class
                    break
            if parseClass is None:
                raise TypeError( 'Invalid node "%s" encountered in suite "%s"' % ( child.tag, self.moniker ) )

            self.add(parseClass.parseNodeUsingClass(child, xPath, linkData, **kwargs))

        xPath.pop( )
        return( self )

class SortedSuite( Suite ) :

    def addIndex( self, item ) :
        """
        Find the location to insert item into a sorted suite.
        :param item: instance to add. Must define a 'sortCompare' method (used to determine insertion point).
        :return: tuple(insertion index (int), whether to replace existing item at that index (boolean))
        """

        for i1, _item in enumerate( self ) :
            cmp = item.sortCompare( _item )
            if( cmp == 0 ) :
                if( self.replace ) : return( i1, True )
                raise KeyError( 'item with key = "%s" already present in non-replaceable suite' % item.key )
            elif( cmp < 0 ) :
                return( i1, False )
        return( len( self ), False )
