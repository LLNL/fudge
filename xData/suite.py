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

class Suite(ancestryModule.AncestryIO_bare):
    """
    Abstract base class for holding a list of objects that can be indexed by an integer or a string.
    For integer indexing, the indexing is by order in which an object was put into the array or,
    if a sortedSuite instance, by its sort order in the suite.  For string indexing, this class acts 
    like the python dict class. 

    This class is like the python OrderedDict class. However, this class does not inherit the OrderedDict 
    class as it does act differently. For example, only objects in the list of **allowedClasses** can be
    added to an suite instance. The **allowedClasses** is a required argument to the constructor.

    The constructor for this class has 2 (not counting self) arguments. They are::

        + ------------------+-----------------------------------------------------------+
        | allowedClasses    | A list of classes. Only objects matching one of the       |
        |                   | classes can be added to the suite.                        |
        + ------------------+-----------------------------------------------------------+
    """

    __metaclass__ = abc.ABCMeta

    def __init__( self, allowedClasses ) :

        ancestryModule.AncestryIO.__init__( self )

        __allowedClasses = []
        for i1, cls in enumerate( allowedClasses ) :
            if( not( inspect.isclass( cls ) ) ) : raise TypeError( 'Item at index %d is not a class.' % i1 )
            __allowedClasses.append( cls )
        self.__allowedClasses = tuple( __allowedClasses )   # Make it a tuple so that it cannot be changed.

        self.__items = []

    def __contains__( self, key ) :
        """
        Returns True if the suite has an item whose **keyValue** is *key* and False otherwise.

        :param str key: The key.
        :rtype:     bool
        """

        if( not( isinstance( key, str ) ) ) : raise TypeError( 'key must be a str instance.' )

        for item in self.__items :
            if( item.keyValue == key ) : return( True )
        return( False )

    def __delitem__( self, key ) :
        """
        Removes an element from self's list. If *key* is an **int**, the item at that index is removed.
        If *key* is a **str**, the item with that **keyValue** is removed.

        :param key: The key of the item to remove.
        :type key: int or str.
        """

        if( isinstance( key, int ) ) :
            del self.__items[key]
            return

        if( not( isinstance( key, str ) ) ) : raise TypeError( 'key must be an int or str instance.' )
        status = self.remove( key )
        if( not( status ) ) : KeyError( 'key not in suite' )

    def __getitem__( self, key ) :
        """
        Returns an item from the suite. If *key* is an int, the item at that index is returned.
        If *key* is a str, the item with that keyValue is returned.

        :param key: The key of the item to return.
        :type key: int or str.
        """

        if( isinstance( key, int ) ) : return( self.__items[key] )

        if( not( isinstance( key, str ) ) ) : raise TypeError( 'key must be a str or an int.' )
        for item in self.__items :
            if( item.keyValue == key ) : return( item )
        raise KeyError( 'key "%s" not found in suite "%s"' % ( key, self.moniker ) )

    def __iter__( self ) :
        """Iterates over self using the index order of items (i.e. int indexing)."""

        n1 = len( self )
        for i1 in range( n1 ) : yield self.__items[i1]

    def __len__( self ) :
        """
        Returns the number of items in self.

        :rtype:     int
        """

        return( len( self.__items ) )

    @property
    def allowedClasses( self ) :
        """Returns the tuple of allowed classes."""

        return( self.__allowedClasses )     # Ok to let user touch this since it cannot be alter as __allowedClasses is a tuple.

    def add( self, item ) :
        """
        Adds **item** to the suite. **item** must be an allowed class, otherwise a **TypeError** exception is raised.

        If another item in the suite has the same keyValue a **KeyError** exception is raised.

        :param item:    The item to add to the suite.
        """

        classAllowed = False
        for cls in self.__allowedClasses :
            if( isinstance( item, cls ) ) :
                classAllowed = True
                break
        if( not( classAllowed ) ) : raise TypeError( 'Invalid class "%s" for suite "%s"' % ( item.__class__, self.moniker ) )

        if( not( isinstance( item.keyValue, str ) ) ) : TypeError( "item.keyValue must be a string." )

        index, exists = self.addIndex( item )
        if( exists ) : raise KeyError( 'item with keyValue = "%s" already present in non-replaceable suite' % item.keyValue )
        item.setAncestor( self )
        self.__items.insert( index, item )

    def addIndex( self, item ) :
        """
        Returns the index at which *item* will be added to *self*.

        :param item:    The item to add to the suite.

        :return:    the tuple *( index, exists )*. *index* is the index in *self* at which *item* will be added. 
                    *exists* is **True** if an *item* with the same **keyValue** already exist and **False** otherwise.
        :rtype:         ( int, bool )
        """

        index = len( self )
        for i1, _item in enumerate( self.__items ) :
            if( _item.keyValue == item.keyValue ) : return( i1, True )

        return( index, False )

    def amendForPatch( self, fromLabel, toLabel ) :

        for item in self.__items : item.amendForPatch( fromLabel, toLabel )

    def checkAncestry( self, verbose = 0, level = 0 ) :

        for item in self : item.checkAncestry( verbose = verbose, level = level )

    def clear( self ):
        """
        Remove all members of the suite.
        """

        self.__items = []

    def convertUnits( self, unitMap ) :
        """
        Calls each items convertUnits function.

        :param unitMap:     A dictionary of keys of "from" units and values of "to" units.
        """

        for index, item in enumerate( self.__items ) : item.convertUnits( unitMap )

    def diff( self, other, diffResults ) :

        keys1 = [ item.keyValue for item in self ]
        keys2 = [ item.keyValue for item in other ]
        keys = [ key for key in keys1 ]
        for key in keyValues2 :
            if( key not in keyValues ) : keyValues.append( key )

        for key in keys :
            if( key not in keys1 ) :
                diffResults.append( '%s missing - 1' % other.moniker, '', '', other[key].toXLink( ) )
            elif( key not in keys2 ) :
                diffResults.append( '%s missing - 2' % self.moniker, '', self[key].toXLink( ), '' )
            else :
                if( type( self ) != type( other ) ) :
                    print( '    Cannot diff type "%s" with type "%s".' % ( type( self ), type( other ) ) )
                else :
                    if( hasattr( self[key], 'diff' ) ) :
                        self[key].diff( other[key], diffResults )
                    else :
                        print( '    Method "diff" missing for object %s' % type( self[key] ) )

    def pop( self, key, *default ) :
        """
        Removes *item* with **keyValue** *key*, and returns it. If *key* is not found, *default* is returned if given, 
        otherwise **KeyError** is raised.

        :param key:         KeyValue of item to Remove.
        :param default:     Item to return if no item with *key* found.

        :return:            Returns the *item* removed or *default* if one is given.
        """

        if( len( default ) > 2 ) : raise Exception( 'Only one default value is allowed: got %s' % len( default ) )

        for i1, item in enumerate( self.__items ) :
            if( item.keyValue == key ) : return( self.__items.pop( i1 ) )
        if( len( default ) == 1 ) : return( default[0] )

        raise KeyError( key )

    def remove( self, key ) :
        """
        Removes item by with keyValue *key*. Returns True if key was present, otherwise returns False.

        :param key: str

        :return: bool
        """

        for i1, item in enumerate( self.__items ) :
            if( item.keyValue == key ) :
                del self.__items[i1]
                return( True )

        return( False )

    def replace( self, item ) :
        """
        Replace an existing item with *item*. If no item in the suite has the same keyValue as *item*, a KeyError is raised.
        If *item* is not an allowed class, a TypeError is raised.

        :param item:
        :return:
        """

        if( not( isinstance( item.keyValue, str ) ) ) : raise ValueError( '''Item's keyValue must be a string.''' )

        found = False
        for cls in self.__allowedClasses :
            if( isinstance( item, cls ) ) :
                found = True
                break
        if( not( found ) ) : raise TypeError( 'Invalid class "%s" for suite "%s"' % ( item.__class__, self.moniker ) )

        for i1, item2 in enumerate( self.__items ) :
            if( item2.keyValue == item.keyValue ) :
                del self.__items[i1]
                self.__items.insert( i1, item )
                item.setAncestor( self, attribute = item.keyName )
                return

        raise KeyError( 'item with %s = "%s" not present in suite' % ( item.keyName, item.keyValue ) )

    def replicate( self, other ) :
        """
        Clear self and then makes a copy of all items in *other* and adds the copy to *self*.
        """

        self.clear( )
        for item in other : self.add( item.copy( ) )

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns an list of XML strings representing *self* and its items.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        showEmptySuites = kwargs.get( 'showEmptySuites', False )

        if( ( len( self ) == 0 ) and not( showEmptySuites ) ) : return( [] )
        XML_stringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        for item in self : XML_stringList += item.toXML_strList( indent = indent2, **kwargs )
        XML_stringList[-1] += '</%s>' % self.moniker

        return( XML_stringList )

    def uniqueKey( self, item ) :
        """
        If item's keyValue is the same as another item's keyValue in self, construct a new unique key
        based on item's keyValue appended with '__' and one or more lower case letters (i.e., 'a' to 'z').
        """

        def incrementSuffix( suffix ) :

            if( len( suffix ) == 0 ) : return( 'a' )
            index = string.ascii_lowercase.find( suffix[-1] ) + 1
            if( index != 26 ) : return( suffix[:-1] + string.ascii_lowercase[index] )
            return( incrementSuffix( suffix[:-1] ) + 'a' )

        if( item.keyValue in self ) :
            key__ = item.keyValue + '__'
            n1 = len( key__ )
            l1 = 0
            suffixes = []
            for _item in self :          # Find list of longest keys that start with key__.
                if( _item.keyValue[:n1] == key__ ) :
                    suffix = _item.keyValue[n1:]
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
            item.keyValue = key__ + suffix

        return( item )

    def parseNode(self, node, xPath, linkData, **kwargs):

        xPath.append( node.tag )

        for child in node :
            parseClass = None
            for _class in self.__allowedClasses :
                if( child.tag == _class.moniker ) :
                    parseClass = _class
                    break
            if( parseClass is None ) :
                raise TypeError( 'Invalid node "%s" encountered in suite "%s".' % ( child.tag, self.moniker ) )

            instance = parseClass.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            self.add( instance )

        xPath.pop( )
        return( self )

class SortedSuite( Suite ) :
    """
    Like the **suite** class except the index at which an *item* is added is determined by the
    *item*'s **sortCompare** method.
    """

    def addIndex( self, item ) :
        """
        Uses *item*'s **sortCompare** method to determine the index in *self* at which *item* will be added.
        Overrides the **addIndex** method in the **suite** class.

        :return:    the tuple *( index, exists )*. *index* is the index in *self* at which *item* will be added. 
                    *exists* is **True** if an *item* with the same **keyValue** already exist and **False** otherwise.
        :rtype:         ( int, bool )
        """

        for i1, _item in enumerate( self ) :
            cmp = item.sortCompare( _item )
            if( cmp == 0 ) :
                return( i1, True )
            elif( cmp < 0 ) :
                return( i1, False )

        return( len( self ), False )

class Suite2( Suite ) :
    """
    Like the **suite** class except the method **toXML** are added.
    """

    pass

class SortedSuite2( SortedSuite ) :
    """
    Like the **sortSuite** class except the method **toXML** are added.
    """

    pass
