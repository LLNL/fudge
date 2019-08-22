# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
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
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

__metaclass__ = type

class ancestry :
    """
    This class is designed to be a base class for a class (instance) that is a member in another
    class (instance). The function of this class is to aid in tracking a class's ancestors. That is,
    if an instance is part of a hierarchy, this class provides methods, for example, that list 
    the position of the instance within the hierarchy or give its position relative to another member 
    in the hierarchy using xlinks.  For example, if rs is a class with moniker 'reactionSuite' containing 
    reactionA and reactionB, each with moniker 'reaction', then
        
        >>>reactionA.getAncestor() -> rs
        >>>reactionB.toXLink() -> "/reactionSuite/reaction[@label='1']"

    This class defines three members:

        moniker     Name of the class.
        ancestor    Instance which self is a child of.
        attribute   Additional qualifier for xlink string.

    For an instance, the xlink is '/the/list/of/ancestors/self' if attribute is None or
    '/the/list/of/ancestors/self[@attribute="value"]' if attribute is not None. For example,
    for a hierarchy consisting of class A with moniker 'nameA' and member mB of class B, class B 
    with moniker 'nameB' and member mC of class C and class C with moniker 'nameC', then the xlink 
    for an instance of a C class in the hierarchy is '/nameA/nameB/nameC'. If the mC instance
    set attribute to be the member 'greeting' that, for this example as value 'Hi', then the
    xlink for the C class is '/nameA/nameB/nameC[@greeting="Hi"]'.
    """

    def __init__( self ) :

        if( not( hasattr( self, 'moniker' ) ) ) : raise Exception( 'class "%s" does not have "moniker" static member' % self.__class__ )
        self.ancestor = None
        self.attribute = None

    def __str__( self ) :

        return( self.toXLink( ) )

    def findAttributeInAncestry( self, attributeName ) :

        if( hasattr( self, attributeName ) ) : return( getattr( self, attributeName ) )
        if( self.ancestor is None ) : raise Exception( 'Could not find attribute name = %s in ancestry' % attributeName )
        return( self.ancestor.findAttributeInAncestry( attributeName ) )

    def findClassInAncestry( self, class_ ) :

        if( isinstance( self, class_ ) ) : return( self )
        if( self.ancestor is None ) : raise Exception( 'Could not find class name = %s in ancestry' % class_.__name__ )
        return( self.ancestor.findClassInAncestry( class_ ) )

    def findEntity( self, entityName, attribute = None, value = None ) :
        """
        Default findEntity method. In general, this method should be over written by sub-class. This method 
        uses the follow algorithm to find entity. Firstly, if 'attribute' is None, then self is assumed to 
        have a attribute named entityName which is taken to be the desired entity. Otherwise, self is iterated 
        over until an item with an attribute named attribute with value value is found. In either case, if 
        an entity is found, its moniker value must be entityName. If no entity is found, raise AttributeError.
        """

        entity = None
        if( entityName in ( '.', self.moniker ) ) :
            entity = self
        elif( entityName == '..' ) :
            entity = self.ancestor
        elif( attribute is None ) :
            entity = getattr( self, entityName )
        else :
            try :                       # try needed in case self cannot be iterated.
                for entity_ in self :
                    if( getattr( entity_, attribute, None ) == value ) :
                        entity = entity_
                        break
            except :
                pass
        if( entity is None or entityName != getattr( entity, 'moniker') ):
            raise AttributeError( "Can't find entity %s in %s" % (entityName,self) )
        return entity

    def getAncestor( self ) :
        """Returns self's ancestor."""

        return( self.ancestor )

    def getRootAncestor( self ) :
        """Traverse up the ancestry tree to the root ancestor and return it. The root ancestor is the instance whose ancestor is None."""

        ancestor = self
        while( ancestor.ancestor is not None ) : ancestor = ancestor.ancestor
        return( ancestor )

    def setAncestor( self, ancestor, attribute = None ) :
        """Sets self's ancestor to ancestor."""

        self.ancestor = ancestor 
        self.attribute = attribute

    def toRelativeXLink( self, other = None ) :
        """
        Returns a string that is a relative xlink to another element (using XML xpath syntax).
        Both elements must reside in the same hierarchy.  For a description of xpath, see 
        http://en.wikipedia.org/wiki/XPath_1.0#Syntax_and_semantics.
        """

        if( other is None ) :
            if( self.ancestor is None ) : return( '' )
            return( self.ancestor.toRelativeXLink( ) + '../' )
        else :
            if( self.getRootAncestor( ) is not other.getRootAncestor( ) ) : 
                raise Exception( 'Root ancestors not the same ("%s" != "%s")' % ( self.toXLink( ), other.toXLink( ) ) )
            thisPath = self.toXLink( ).split( '/' )
            othersPath = other.toXLink( ).split( '/' )
            for i1, tag in enumerate( thisPath ) :
                if( i1 >= len( othersPath ) ) : break
                if( tag != othersPath[i1] ) : break
            relativePath = ''
            for i2 in xrange( len( thisPath ) - i1) : relativePath += '../'
            relativePath += '/'.join( othersPath[i1:] )
            return( relativePath )

    def toXLink( self, attributeName = None, attributeValue = None ) :
        """
        Returns a string that is an xlink to self (using XML xpath syntax).  The resulting 
        xlink starts at the root element. For a description of xpath, see 
        http://en.wikipedia.org/wiki/XPath_1.0#Syntax_and_semantics.
        """

        s1, attribute = '', ''
        if( self.ancestor is not None ) : s1 = self.ancestor.toXLink( )
        if( ( attributeName is None ) and ( self.attribute is not None ) ) : 
            attributeName, attributeValue = self.attribute, getattr( self, self.attribute )
        if( attributeName is not None ) :
            if( attributeValue is None ) : raise Exception( 'attributeValue is None while attributeName is not' )
            attribute = "[@%s='%s']" % ( attributeName, attributeValue )
        elif( attributeValue is not None ) :
            raise Exception( 'attribute name is None but value is not: value = %s' % attributeValue )

        return( s1 + '/%s%s' % ( self.moniker, attribute ) )

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
    print str( p )
    print str( c )
    print str( gc1 )
    gc2 = granddaughter( 'Tami' )
    c.addChild( gc2 )
    print gc2

    print p.findEntity( 'child', 'name', 'Mary' )
    print c.findEntity( 'child', 'name', 'Tom' )
    print c.findEntity( 'granddaughter', 'name', 'Tami' )
    print c.findEntity( 'grandson', 'name', 'Joe' )
