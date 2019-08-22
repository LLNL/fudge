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

import abc

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

    __metaclass__ = abc.ABCMeta
    ancestryMembers = tuple( )

    def __init__( self ) :

        self.ancestor = None
        self.attribute = None

    def __str__( self ) :

        return( self.toXLink( ) )

    @abc.abstractproperty
    def moniker( self ) :

        pass

    def checkAncestry( self, verbose = 0, level = 0 ) :

        def check( child ) :

            if( child is None ) : return
            if( self.isChild( child ) ) :
                child.checkAncestry( verbose = verbose, level = level )
            else :
                print 'WARNING from checkAncestry: member "%s" not a child of %s' % ( member, self.toXLink( ) )
                print '    Its ancestry is: %s' % child.toXLink( )

        if( self.ancestryMembers == ( '', ) ) : return

        prefix = ( level + 1 ) * '    '
        if( len( self.ancestryMembers ) == 0 ) :
            if( verbose != 0 ) : print "%s---- no items in ancestryMembers for %s" % ( prefix, self.toXLink( ) )
        for member in self.ancestryMembers :
            if( member == '' ) : continue
            if( verbose > 0 ) : print "%s%s" % ( prefix, member )
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
                print 'WARNING from checkAncestry: %s does not have member "%s"' % ( self.toXLink( ), member )

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
        Default findEntity method. In general, sub-classes should over-ride this method.
        This method uses the following algorithm to find entity. Firstly, if 'attribute' is None, then self is assumed
        to have an attribute named entityName which is taken to be the desired entity. Otherwise, self is iterated
        over until an item with an attribute named attribute with value value is found. In either case, if
        an entity is found, its moniker value must be entityName. If no entity is found, raise AttributeError.
        """

        if( entityName in ( '.', self.moniker ) ) :
            return self
        elif( entityName == '..' ) :
            return self.ancestor
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

    def getAncestor( self ) :
        """Returns self's ancestor."""

        return( self.ancestor )

    def getRootAncestor( self ) :
        """Traverse up the ancestry tree to the root ancestor and return it. The root ancestor is the instance whose ancestor is None."""

        ancestor = self
        while( ancestor.ancestor is not None ) : ancestor = ancestor.ancestor
        return( ancestor )

    def isChild( self, child ) :

        if( isinstance( child, ancestry ) ) : return( child.ancestor == self )
        return( False )

    def isParent( self, parent ) :

        return( self.ancestor == parent )

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

    def followXPath( self, xPath ):
        """
        :param xPath: string xPath, e.g. "/reactionSuite/reaction[@label='2']"
        :return: class instance pointed to by xPath

        Uses ancestry.findEntity to find each element
        """
        def follow2( xPathList, node ):
            # recursive helper function: descend the path to find the correct element
            if len(xPathList)==0: return node

            xPathNext = xPathList[0]
            try:
                if "[@" in xPathNext:
                    r1,r2 = xPathNext.split("[@")
                    r2,r3 = r2.split("=")
                    r3 = r3.rstrip(']')[1:-1]
                    nodeNext = node.findEntity( r1, r2, r3 )
                else:
                    nodeNext = node.findEntity( xPathNext )
            except:
                raise XPathNotFound()
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
        except XPathNotFound:
            raise XPathNotFound( "Cannot locate path '%s'" % xPath )


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

    print
    print 'Checking ancestry:'
    p.checkAncestry( )
