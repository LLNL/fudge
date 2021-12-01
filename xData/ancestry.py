# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys
import abc

from xml.etree import cElementTree

from . import formatVersion as formatVersionModule

__metaclass__ = type

class ancestry :
    """
    This class is designed to be a base class for a class (instance) that is a member in another
    class (instance). The function of this class is to aid in tracking a class's ancestors. That is,
    if an instance is part of a hierarchy, this class provides methods, for example, that list 
    the position of the instance within the hierarchy or give its position relative to another member 
    in the hierarchy using xlinks.  For example, if rs is a class with moniker 'reactionSuite' containing 
    reactionA and reactionB, each with moniker 'reaction', then
        
        >>>reactionA.ancestor -> rs
        >>>reactionB.toXLink() -> "/reactionSuite/reaction[@label='1']"

    This class defines three members:

        moniker                     Name of the class.
        ancestor                    Instance which self is a child of.
        attribute                   Additional qualifier for xlink string.
        keyName                     The name of the key for instance if it goes in a suite.
        ancestryMembers             Tuple of names of members added by sub-classes.
        legacyMemberNameMapping     A map whose keys are legacy member names and whose associated values are the current memeber names.

    For an instance, the xlink is '/the/list/of/ancestors/self' if attribute is None or
    '/the/list/of/ancestors/self[@attribute="value"]' if attribute is not None. For example,
    for a hierarchy consisting of class A with moniker 'nameA' and member mB of class B, class B 
    with moniker 'nameB' and member mC of class C and class C with moniker 'nameC', then the xlink 
    for an instance of a C class in the hierarchy is '/nameA/nameB/nameC'. If the mC instance
    set attribute to be the member 'greeting' that, for this example as value 'Hi', then the
    xlink for the C class is '/nameA/nameB/nameC[@greeting="Hi"]'.
    """

    ancestryMembers = tuple( )
    legacyMemberNameMapping = {}
    keyName = None
    formatVersion = formatVersionModule.default
    monikerByFormat = {}

    def __init__( self ) :

        self.__ancestor = None
        self.__attribute = None

    def __str__( self ) :

        return( self.toXLink( ) )

    @property
    @abc.abstractmethod
    def moniker( self ) :

        pass

    @property
    def ancestor( self ) :
        """Returns self's ancestor."""

        return( self.__ancestor )

    @property
    def attribute( self ) :

        return( self.__attribute )

    @property
    def keyValue( self ) :
        """Returns self's keyValue."""

        if( self.keyName is None ) :
            if( hasattr( self, 'label' ) ) : return( self.label )           # For legacy, this is deprecated.
            return( None )
        return( getattr( self, self.keyName ) )

    @property
    def rootAncestor( self ) :
        """Traverse up the ancestry tree to the root ancestor and return it. The root ancestor is the instance whose ancestor is None."""

        ancestor = self
        while( ancestor.__ancestor is not None ) : ancestor = ancestor.__ancestor
        return( ancestor )

    def checkAncestry( self, verbose = 0, level = 0 ) :

        def check( child ) :

            if( child is None ) : return
            if( self.isChild( child ) ) :
                child.checkAncestry( verbose = verbose, level = level )
            else :
                print('WARNING from checkAncestry: member "%s" not a child of %s' % (member, self.toXLink()))
                print('    Its ancestry is: %s' % child.toXLink())

        if( self.ancestryMembers == ( '', ) ) : return

        prefix = ( level + 1 ) * '    '
        if( len( self.ancestryMembers ) == 0 ) :
            if( verbose != 0 ) : print("%s---- no items in ancestryMembers for %s" % (prefix, self.toXLink()))
        for member in self.ancestryMembers :
            if( member == '' ) : continue
            if( verbose > 0 ) : print("%s%s" % (prefix, member))
            doLoop = False
            if( member[0] == '[' ) :
                member = member[1:]
                doLoop = True
            if( hasattr( self, member ) ) :
                m1 = getattr( self, member )
                if( doLoop ) :
                    for child in m1 : check( child )
                else :
                    check( m1 )
            else :
                print('WARNING from checkAncestry: %s does not have member "%s"' % (self.toXLink(), member))

    def copy( self ):
        """
        Use deepcopy, but don't copy self.ancestor
        :return: copy of self
        """
        import copy
        memodict = {id(self.ancestor):None} # TODO: add 'unresolvedLinks' to memodict
        return copy.deepcopy(self, memo=memodict)

    def findAttributeInAncestry( self, attributeName ) :

        if( hasattr( self, attributeName ) ) : return( getattr( self, attributeName ) )
        if( self.__ancestor is None ) : raise Exception( 'Could not find attribute name = %s in ancestry' % attributeName )
        return( self.__ancestor.findAttributeInAncestry( attributeName ) )

    def findClassInAncestry( self, class_ ) :

        if( isinstance( self, class_ ) ) : return( self )
        if( self.__ancestor is None ) : raise Exception( 'Could not find class name = %s in ancestry' % class_.__name__ )
        return( self.__ancestor.findClassInAncestry( class_ ) )

    def findEntity( self, entityName, attribute = None, value = None ) :
        """
        Default findEntity method. In general, sub-classes should over-ride this method.
        This method uses the following algorithm to find entity. Firstly, if 'attribute' is None, then self is assumed
        to have an attribute named entityName which is taken to be the desired entity. Otherwise, self is iterated
        over until an item with an attribute named attribute with value value is found. In either case, if
        an entity is found, its moniker value must be entityName. If no entity is found, raise AttributeError.
        """

        if entityName in self.legacyMemberNameMapping: entityName = self.legacyMemberNameMapping[entityName]

        if( entityName in ( '.', self.moniker ) ) :
            return self
        elif( entityName == '..' ) :
            return self.__ancestor
        entity = None
        if( attribute is None ) :
            entity = getattr( self, entityName )
        else :
            try :                       # try needed in case self cannot be iterated.
                for entity_ in iter(self) :
                    if( str( getattr( entity_, attribute, None ) ) == value ) :
                        entity = entity_
                        break
            except TypeError:
                pass
        if( entity is None or entityName != getattr( entity, 'moniker') ):
            raise AttributeError( "Can't find entity %s in %s" % (entityName,self) )
        return entity

    def findInstancesOfClassInChildren( self, cls, level = 9999 ) :
        """
        Finds all instances of class *cls* in self's children, grand-children, etc.
        """

        foundInstances = []
        level -= 1
        if( level < 0 ) : return( foundInstances )
        for ancestryMember in self.ancestryMembers :
            if( ancestryMember == '' ) : continue
            if( ancestryMember[0] == '[' ) : ancestryMember = ancestryMember[1:]
            instance = getattr( self, ancestryMember )
            if( instance is None ) : continue
            if( isinstance( instance, cls ) ) : foundInstances.append( instance )
            foundInstances += instance.findInstancesOfClassInChildren( cls, level = level )

        return( foundInstances )

    def getRootAncestor( self ) :
        """This function is deprecated. See and use rootAncestor instead."""

        return( self.rootAncestor )

    def isChild( self, child ) :

        if( isinstance( child, ancestry ) ) : return( child.__ancestor == self )
        return( False )

    def isParent( self, parent ) :

        return( self.__ancestor == parent )

    def setAncestor( self, ancestor, attribute = None ) :
        """Sets self's ancestor to ancestor."""

        self.__ancestor = ancestor 
        self.__attribute = attribute

    def toRelativeXLink( self, other = None, formatVersion = None ) :
        """
        Returns a string that is a relative xlink to another element (using XML xpath syntax).
        Both elements must reside in the same hierarchy.  For a description of xpath, see 
        http://en.wikipedia.org/wiki/XPath_1.0#Syntax_and_semantics.
        """

        if( other is None ) :
            if( self.__ancestor is None ) : return( '' )
            return( self.__ancestor.toRelativeXLink( ) + '../' )
        else :
            if( self.getRootAncestor( ) is not other.getRootAncestor( ) ) : 
                raise Exception( 'Root ancestors not the same ("%s" != "%s")' % ( self.toXLink( ), other.toXLink( ) ) )
            thisPath = self.toXLink( formatVersion = formatVersion ).split( '/' )
            othersPath = other.toXLink( formatVersion = formatVersion ).split( '/' )
            for i1, tag in enumerate( thisPath ) :
                if( i1 >= len( othersPath ) ) : break
                if( tag != othersPath[i1] ) : break
            relativePath = ''
            for i2 in range( len( thisPath ) - i1) : relativePath += '../'
            relativePath += '/'.join( othersPath[i1:] )
            return( relativePath )

    @abc.abstractmethod
    def toXMLList( self, **kwargs ) :
        pass

    def toXLink( self, attributeName = None, attributeValue = None, formatVersion = None ) :
        """
        Returns a string that is an xlink to self (using XML xpath syntax).  The resulting 
        xlink starts at the root element. For a description of xpath, see 
        http://en.wikipedia.org/wiki/XPath_1.0#Syntax_and_semantics.
        """

        s1, attribute = '', ''
        if( self.__ancestor is not None ) : s1 = self.__ancestor.toXLink( formatVersion = formatVersion )
        if( ( attributeName is None ) and ( self.__attribute is not None ) ) : 
            attributeName, attributeValue = self.__attribute, getattr( self, self.__attribute )
        if( attributeName is not None ) :
            if( attributeValue is None ) : raise Exception( 'attributeValue is None while attributeName is not' )
            attribute = "[@%s='%s']" % ( attributeName, attributeValue )
        elif( attributeValue is not None ) :
            raise Exception( 'attribute name is None but value is not: value = %s' % attributeValue )
        moniker = self.monikerByFormat.get(formatVersion, self.moniker)

        return( s1 + '/%s%s' % ( moniker, attribute ) )

    @abc.abstractmethod
    def parseXMLNode( self, element, xPath, linkData ) :
        pass

    def followXPath( self, xPath ):
        """
        :param xPath: string xPath, e.g. "/reactionSuite/reaction[@label='2']"
        :return: class instance pointed to by xPath

        Uses ancestry.findEntity to find each element
        """

        def follow2( xPathList, node ):
            """For internal use. Recursive helper function: descend the path to find the correct element."""

            if len(xPathList)==0: return node

            xPathNext = xPathList[0]

            try:
                if "[@" in xPathNext:
                    r1,r2 = xPathNext.split("[@",1)
                    r2,r3 = r2.split("=",1)
                    r3 = r3.rstrip(']')[1:-1]
                    nodeNext = node.findEntity( r1, r2, r3 )
                else:
                    nodeNext = node.findEntity( xPathNext )
            except:
                raise XPathNotFound( )
            return follow2(xPathList[1:], nodeNext)

        # FIXME refactor to use xml.etree.ElementPath.xpath_tokenizer?
        xPathList = xPath.split('/')

        while not xPathList[0]: # trim empty sections from the beginning
            xPathList = xPathList[1:]
        if '{' in xPath:        # careful, qualifiers may contain '/'
            xpl2 = []
            for val in xPathList:
                if '}' in val and '{' not in val and '{' in xpl2[-1] and '}' not in xpl2[-1]:
                    xpl2[-1] += '/' + val
                else:
                    xpl2.append( val )
            xPathList = xpl2
        try:
            return follow2( xPathList, self )
        except XPathNotFound :
            raise XPathNotFound( "Cannot locate path '%s'" % xPath )

class Ancestry2( ancestry ) :

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def parseAncestryMembers( self, node, xPath, linkData, **kwargs ) :

        for child in node :
            for memberName in self.ancestryMembers :
                if( memberName == '' ) : continue
                if( memberName[:1] == '[' ) : memberName = memberName[1:]
                member = getattr( self, memberName )
                if( child.tag == member.moniker ) :
                    member.parseNode( child, xPath, linkData, **kwargs )
                    break

    def parseNode( self, node, xPath, linkData, **kwargs ) :

        xPath.append( node.tag )

        self.parseAncestryMembers( node, xPath, linkData, **kwargs )

        xPath.pop( )

    @classmethod
    def parseXMLString( cls, string, **kwargs ) :
        """
        Parses a XML string using class cls.
        """

        from LUPY import xmlNode as xmlNodeMode        # Wrapper around the xml parser.

        types = ( str, )
        if( not( isinstance( string, types ) ) ) : raise TypeError( 'Invalid string.' )

        node = cElementTree.fromstring( string )
        node = xmlNodeMode.xmlNode( node, xmlNodeMode.xmlNode.etree )

        return( cls.parseNodeUsingClass( node, [], {}, **kwargs ) )

    @classmethod
    def readXMLFile( cls, fileName, **kwargs ) :
        """
        Reads and parses an XML file using class cls.
        """

        from LUPY import xmlNode as xmlNodeMode        # Wrapper around the xml parser.

        if( not( isinstance( fileName, str ) ) ) : raise TypeError( 'Invalid file name.' )

        node = cElementTree.parse( fileName ).getroot( )
        node = xmlNodeMode.xmlNode( node, xmlNodeMode.xmlNode.etree )

        return( cls.parseNodeUsingClass( node, [], {}, **kwargs ) )

    @classmethod
    def parseNodeUsingClass( cls, node, xPath, linkData, **kwargs ) :
        """
        """

        xPath.append( node.tag )
        instance = cls.parseConstructBareNodeInstance( node, xPath, linkData, **kwargs )
        xPath.pop( )

        instance.parseNode( node, xPath, linkData, **kwargs )

        return( instance )

    @staticmethod
    def parseConstructBareNodeInstance( node, xPath, linkData, **kwargs ) :
        """
        For internal use.  This class should only be called from parseNodeUsingClass.  This class method must be overriden by 
        sub-class. The overrided method should only parse the parts of node needed to construct cls and then return the instance.
        """

        raise Exception( 'This static method of "%s" must be overloaded by sub-class and is not for moniker "%s".' % ( Ancestry2, node.tag ) )

class XPathNotFound( Exception ):

    pass

if( __name__ == '__main__' ) :

    class person( ancestry ) :

        def __init__( self, name ) :

            ancestry.__init__( self )
            self.name = name
            self.children = []

        def __getitem__( self, index ) :

            return( self.children[index] )

        def addChild( self, child ) :

            self.children.append( child )
            self.children[-1].setAncestor( self, attribute = 'name' )

    class parent( person ) :

        moniker = 'parent'

        def __init__( self, name ) :

            person.__init__( self, name )

    class child( person ) :

        moniker = 'child'

        def __init__( self, name ) :

            person.__init__( self, name )

    class grandson( person ) :

        moniker = 'grandson'

        def __init__( self, name ) :

            person.__init__( self, name )

        def toXLink( self ) :

            return( ancestry.toXLink( self, 'hi', 'bye' ) )

    class granddaughter( person ) :

        moniker = 'granddaughter'

        def __init__( self, name ) :

            person.__init__( self, name )

    p = parent( 'Fred' )
    c = child( 'Mary' )
    p.addChild( c )
    gc1 = grandson( 'Joe' )
    c.addChild( gc1 )
    print(str(p))
    print(str( c ))
    print(str( gc1 ))
    gc2 = granddaughter( 'Tami' )
    c.addChild( gc2 )
    print(gc2)

    print(p.findEntity( 'child', 'name', 'Mary' ))
    print(c.findEntity( 'child', 'name', 'Tom' ))
    print(c.findEntity( 'granddaughter', 'name', 'Tami' ))
    print(c.findEntity( 'grandson', 'name', 'Joe' ))

    print()
    print('Checking ancestry:')
    p.checkAncestry( )
