# <<BEGIN-copyright>>
# <<END-copyright>>

__metaclass__ = type

class ancestry :
    """This class keeps track of its 'parent' class. Thus, in a nested hierarchy (such as GND),
    by inheriting from ancestry child classes know where they live in the hierarchy,
    can print an 'xLink' indicating their position in the hierarchy, and also have access to their parent.
    For example, if rs is a reactionSuite containing reactionA and reactionB,
        
        >>>reactionA.getParent() -> rs
        >>>reactionB.toXLink() -> "/reactionSuite/reaction[@label='1']"

    Nearly all classes in GND inherit from this class.
    """

    def __init__( self, moniker, parent, attribute = None ) :

        if( not( hasattr( self, 'name' ) ) ) : self.name = moniker
        if( not( hasattr( self, 'moniker' ) ) ) : self.moniker = moniker
        self.parent = parent
        self.attribute = attribute

    def __str__( self ) :

        return( self.toXLink( ) )

    def findAttributeInAncestry( self, attributeName ) :

        if( hasattr( self, attributeName ) ) : return( getattr( self, attributeName ) )
        if( self.parent is None ) : raise Exception( 'Could not find attribute name = %s in ancestry' % attributeName )
        return( self.parent.findAttributeInAncestry( attributeName ) )

    def findClassInAncestry( self, class_ ) :

        if( isinstance( self, class_ ) ) : return( self )
        if( self.parent is None ) : raise Exception( 'Could not find class name = %s in ancestry' % class_.__name__ )
        return( self.parent.findClassInAncestry( class_ ) )

    def getParent( self ) :
        """Returns self's parent."""

        return( self.parent )

    def getRootParent( self ) :
        """Traverse up the ancestry tree to the parent that does not have a parent, and returns it."""

        parent = self
        while( parent.parent is not None ) : parent = parent.parent
        return( parent )

    def setParent( self, parent ) :
        """Sets self's parent to parent."""

        self.parent = parent

    def toRelativeXLink( self, other = None ) :
        """Returns a string that is a relative xlink to another element (using XML xpath syntax).
        Both elements must reside in the same hierarchy.
        For a description of xpath, see http://en.wikipedia.org/wiki/XPath_1.0#Syntax_and_semantics."""

        if( other is None ) :
            if( self.parent is None ) : return( '' )
            return( self.parent.toRelativeXLink( ) + '../' )
        else :
            if( self.getRootParent( ) is not other.getRootParent( ) ) : raise Exception( 'Root ancesties not the same ("%s" != "%s")' % ( self.toXLink( ), other.toXLink( ) ) )
            thisPath = self.toXLink( ).split( '/' )
            othersPath = other.toXLink( ).split( '/' )
            for i, tag in enumerate( thisPath ) :
                if( i >= len( othersPath ) ) : break
                if( tag != othersPath[i] ) : break
            relativePath = ''
            for j in xrange( len( thisPath ) - i ) : relativePath += '../'
            relativePath += '/'.join( othersPath[i:] )
            thisPath = '/'.join( thisPath[i:] )
            return( relativePath )

    def toXLink( self, attributeName = None, attributeValue = None ) :
        """Returns a string that is an xlink to self (using XML xpath syntax).
        The resulting xlink starts at the root element.
        For a description of xpath, see http://en.wikipedia.org/wiki/XPath_1.0#Syntax_and_semantics."""

        s, attribute = '', ''
        if( self.parent is not None ) : s = self.parent.toXLink( )
        if( ( attributeName is None ) and ( self.attribute is not None ) ) : attributeName, attributeValue = self.attribute, getattr( self, self.attribute )
        if( attributeName is not None ) :
            attribute = "[@%s='%s']" % ( attributeName, attributeValue )
        elif( attributeValue is not None ) :
            raise Exception( 'attribute name is None but value is not: value = %s' % attributeValue )

        return( s + '/%s%s' % ( self.moniker, attribute ) )

if( __name__ == '__main__' ) :

    class parent( ancestry ) :

        def __init__( self ) :

            ancestry.__init__( self, 'parent', None )

    class child( ancestry ) :

        def __init__( self, parent ) :

            ancestry.__init__( self, 'child', parent )

    class grandson( ancestry ) :

        def __init__( self, parent ) :

            ancestry.__init__( self, 'grandson', parent )

        def toXLink( self ) :

            return( ancestry.toXLink( self, 'hi', 'bye' ) )

    class granddaughter( ancestry ) :

        def __init__( self, parent ) :

            ancestry.__init__( self, 'child', parent, attribute = 'hello' )
            self.hello = 'ciao'

    p = parent( )
    c = child( p )
    gc1 = grandson( c )
    print str( p )
    print str( c )
    print str( gc1 )
    gc2 = granddaughter( c )
    print gc2
