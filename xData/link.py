# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__metaclass__ = type

from . import base as baseModule
from . import ancestry as ancestryModule

class link( baseModule.xDataCoreMembers ) :
    """
    This class contains the path to another element, that could be stored in a different file.
    The 'follow' member function is used to get a pointer to the linked-to element,
    which will then be stored as self.link.
    
    Member data ::
    
        * label : the xml tag name (type: str)
        * link  : actual pointer to other data (type: str)
        * root  : file name where other data is stored, only needed if it's not the current file (type: str)
        * path  : self.path is an xPath expression that uniquely identifies the other data (type: str)

    sample links::
    <element href='/reactionSuite/reaction[@label="102"]'>
        or
    <element href='/covariances.xml#covId'>
    """

    moniker = 'link'
    ancestryMembers = ( '', )

    def __init__( self, link = None, root = None, path = None, label = None, relative = False ) :
        """
        Creates a link to another instance.

        :param link: pointer to another instance. That instance should inherit from ancestry
        :param root: external file identifier, required if the link points outside the current file
        :param path: xPath expression that uniquely identifies the other instance
        :param label: unique label, generally associated with a style (e.g. 'eval')
        :param relative: boolean, whether to use relative link when writing out xPath (default = False)
        """

        baseModule.xDataCoreMembers.__init__( self, self.moniker, index = None, label = label )
        self.link = link
        self.root = root
        self.path = path
        self.__relative = relative

            # careful with nesting single/double quotes.
        if self.path is not None: self.path = self.path.replace('"',"'")

    def __str__( self ) :

        return( "%s" % self.path )

    def __deepcopy__( self, memodict={} ):

        if( self.path == None ) : self.updateXPath( )
        link_ = self.link
        if id(link_) in memodict:
            link_ = memodict[id(link_)]
        copy_ = self.__class__( link = link_, root = self.root, path = self.path, label = self.label,
            relative = self.__relative )
        return( copy_ )

    copy = __deepcopy__

    def build_href( self, **kwargs ) :
        """Builds the href (using "formatVersion" if present in **kwargs) and returns the href. The href will include the root if the root is not None."""

        if( self.link is None ) :           # Should only happen when linking is to another protare's data.
            XPath = self.path
            if( self.root is not None ) : XPath = '%s#%s' % ( self.root, XPath )
        else :
            formatVersion = kwargs.get( 'formatVersion' )
            if self.__relative and hasattr(self, 'toRelativeXLink'):
                XPath = self.toRelativeXLink( self.link, formatVersion = formatVersion )
            else:
                XPath = self.link.toXLink( formatVersion = formatVersion )
                if( self.root is not None ) : XPath = '%s#%s' % ( self.root, XPath )

        return( XPath )

    def updateXPath( self ):
        """Ensure the xPath agrees with the linked-to data."""

        if self.__relative and hasattr(self, 'toRelativeXLink'):
            self.path = self.toRelativeXLink( self.link )
        elif hasattr(self.link, 'toXLink'):
            self.path = self.link.toXLink()

    def follow( self, startNode ):
        """
        :param startNode: instance corresponding to the beginning of self.path (must inherit from ancestry)
        :return: class instance pointed to by self.path
        """
        return startNode.followXPath( self.path )

    def toXML( self, indent = '' , **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        attributesStr = self.attributesToXMLAttributeStr( )
        return( [ '%s<%s%s href="%s"/>' % ( indent, self.moniker, attributesStr, self.build_href( **kwargs ) ) ] )

    @classmethod
    def parseXMLNode( cls, linkElement, xPath, linkData ):
        """
        Parse the xml-represented link back to python. The resulting link points to None,
        and must be resolved by the calling function.
        Link attributes can be converted to desired type by passing a 'typeConversion' dictionary in linkData
        """

        xPath.append( linkElement.tag )
        path = linkElement.get('href')
        if '#' in path: root, path = path.split('#')
        else: root = None
        relative = not path.startswith('/')
        # derived classes may add new (non-xlink) attributes:
        kwargs = dict( [v for v in list( linkElement.items( ) ) if not v[0].startswith( 'href' )] )
        for key in kwargs:
            if key in linkData.get('typeConversion',{}):
                kwargs[key] = linkData['typeConversion'][key](kwargs[key])
        result = cls( link=None, root=root, path=path, relative=relative, **kwargs )
        if( 'unresolvedLinks' in linkData ) : linkData['unresolvedLinks'].append( result )
        xPath.pop()

        return result

class link2( ancestryModule.ancestry ) :
    """
    This class contains an href to another node as a xPath. The href could reference a node in a different file.
    The 'follow' member function is used to get a pointer to the linked-to node,
    which will then be stored as self.instance.
    
    Member data ::
    
        * moniker   : the moniker of the instance pointed to by instance (type: str)
        * instance  : python reference to the linked instance (type: reference)
        * root      : file name where other node exists, only needed if it's not the current file (type: str)
        * href      : an xPath expression that uniquely identifies the other instance (type: str)
        * relative  : Whether to use relative xlink when writing out xPath (type: bool)

    sample links::
    <element href='/reactionSuite/reaction[@label="102"]'/>
        or
    <element href='/covariances.xml#covId'/>
    """

    moniker = 'link'
    ancestryMembers = ( '', )

    def __init__( self, moniker, instance = None, root = None, href = None, relative = False, keyName = None, keyValue = None, label = None ) :

        self.__moniker = moniker
        ancestryModule.ancestry.__init__( self )

        self.__instance = instance      # Pointer to linked instance.
        self.__root = root              # File name where other data is stored, only needed if it's not the current file
        self.__href = href              # xPath expression to linked instance.
        self.__relative = relative      # Whether to use relative xlink when writing out xPath
        self.__keyName = None
        self.__keyValue = None
        self.updateKey( keyName = keyName, keyValue = keyValue )
        self.__label = label

# BRB is this what we want or something else.
        if( self.__href is not None ) : self.__href = self.__href.replace( '"', "'" )        # Careful with nesting single/double quotes.

    def __str__( self ) :

        return( "%s" % self.__href )

    @property
    def moniker( self ) :

        return( self.__moniker )

    @property
    def instance( self ) :

        return( self.__instance )

    @instance.setter
    def instance( self, value ) :

        self.__instance = value

    @property
    def link( self ) :

        return( self.__instance )

    @link.setter
    def link( self, value ) :

        self.instance = value

    @property
    def keyName( self ) :

        return( self.__keyName )

    @property
    def keyValue( self ) :

        return( self.__keyValue )

    @property
    def label( self ) :

        return( self.__label )

    @property
    def root( self ) :

        return( self.__root )

    @property
    def href( self ) :

        return( self.__href )

    @property
    def path( self ) :

        return( self.__href )

    @property
    def relative( self ) :

        return( self.__relative )

    def copy( self, *args, **kwargs ) :

        try :
            if(self.__href is None) : self.updateXPath()
        except :
            pass
        return( link2( self.__moniker, instance = self.__instance, root = self.__root, href = self.__href, keyName = self.__keyName, 
                keyValue = self.__keyValue, label = self.__label ) )

    def build_href( self, **kwargs ) :
        """Builds the href (using "formatVersion" if present in **kwargs) and returns the href. The href will include the root if the root is not None."""

        if( self.link is None ) :           # Should only happen when linking is to another protare's data.
            XPath = self.path
            if( self.root is not None ) : XPath = '%s#%s' % ( self.root, XPath )
        else:
            formatVersion = kwargs.get( 'formatVersion' )
            if self.__relative and hasattr(self, 'toRelativeXLink'):
                XPath = self.toRelativeXLink( self.instance, formatVersion = formatVersion )
            else:
                XPath = self.instance.toXLink( formatVersion = formatVersion )
                if( self.root is not None ) : XPath = '%s#%s' % ( self.root, XPath )

        return( XPath )

    def updateKey( self, keyName = None, keyValue = None ) :

        if( self.__keyName is not None ) :
            if( self.__keyName != keyName ) : raise Exception( '''self's Keyname "%s" != "%s"''' % ( self.__keyName, keyName ) )

        if( self.__instance is not None ) :
            if( hasattr( self.__instance, 'keyName' ) ) :
                if( keyName is None ) :
                    keyName = self.instance.keyName
                else :
                    if( self.instance.keyName != keyName ) : raise Exception( '''Keyname "%s" != "%s"''' % ( self.instance.keyName, keyName ) )

        self.__keyName = keyName
        self.__keyValue = keyValue

    def updateXPath( self ) :
        """
        Ensure the xPath agrees with the linked instance.
        """

        if( self.__relative ) :
            self.__href = self.toRelativeXLink( self.__instance )
        else :
            if( self.__instance is not None ) :         # FIXME
                self.__href = self.__instance.toXLink( )

    def follow( self, startNode ) :
        """
        :param startNode: instance corresponding to the beginning of self.href (must inherit from ancestry)
        :return: instance pointed to by self.href
        """

        return( startNode.followXPath( self.__href ) )

    def toXML( self, indent = '' , **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :
        """
        Pointers show up in the attributes list on an xml element (i.e., <element href="...xpath" anotherAttribute="foo" .../>).
        """

        attributes = ''
        if( ( self.__keyName is not None ) and ( self.__keyValue is not None ) ) : attributes += ' %s="%s"' % ( self.__keyName, self.__keyValue )
        if( self.__label is not None ) : attributes += ' label="%s"' % self.__label

        return( [ '%s<%s%s href="%s"/>' % ( indent, self.moniker, attributes, self.build_href( **kwargs ) ) ] )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :
        """
        Parse the xml-represented link back to python. The resulting link points to None, and must be resolved by the calling function.
        Link attributes can be converted to desired type by passing a 'typeConversion' dictionary in linkData
        """

        xPath.append( element.tag )

        root = None
        href = None
        keyName = None
        keyValue = None
        label = None

        for key in list( element.keys( ) ) :
            if( key == 'href' ) :
                href = element.get( key )
            elif( key == 'label' ) :
                label = element.get( key )
            else :
                if( keyName is not None ) : raise Exception( 'Too many attributes for link: %s' % key )
                keyName = key
                keyValue = element.get( key )

        if( '#' in href ) : root, href = href.split( '#' )

        _link = cls( element.tag, instance = None, root = root, href = href, relative = href[0] == '/', keyName = keyName, keyValue = keyValue,
            label = label )

        if( 'unresolvedLinks' in linkData ) : linkData['unresolvedLinks'].append( _link )

        xPath.pop( )
        return( _link )
