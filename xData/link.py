# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from . import base as baseModule

class Link(baseModule.XDataCoreMembers):
    """
    An instance of this class represents a GNDS xPath and a link to to another instance (typically an instance of this class
    represents a form that links to another form).  The :py:func:`follow` member function is used to get a pointer to the 
    linked-to instance, which will then be stored as self's link member.

    The following table list the primary members of this class:

    +-----------+-----------------------------------------------------------+
    | Member    | Description                                               |
    +===========+===========================================================+
    | label     | XML tag name.                                             |
    +-----------+-----------------------------------------------------------+
    | link      | Actual pointer to other instance (form).                  |
    +-----------+-----------------------------------------------------------+
    | root      | File name where the other instance is stored. Only used   |
    |           | if the linked-to instance was read / will be stored in a  |
    |           | different file.                                           |
    +-----------+-----------------------------------------------------------+
    | path      | An xPath expression that identifies the other data.       |
    +-----------+-----------------------------------------------------------+
    | relative  | If True, the xPath is relative to *self*, otherwise       |
    |           | it is an absolute xPath.                                  |
    +-----------+-----------------------------------------------------------+
    
    example links::

        <form href='/reactionSuite/reaction[@label="102"]'>
        or
        <form href='/covariances.xml#covId'>
    """

    moniker = 'link'

    def __init__(self, link=None, root=None, path=None, label=None, relative=False):
        """
        Creates a link to another instance.

        :param link:        This is a reference to the linked-to instance. That instance must inherit from ancestry.
        :param root:        External file identifier, required if the link points outside the current file.
        :param path:        This is the xPath expression that uniquely identifies the lined-to instance.
        :param label:       The label to *self*.
        :param relative:    If True, use relative xPath when writing to a file (default = False).
        """

        baseModule.XDataCoreMembers.__init__(self, index=None, label=label)

        self.link = link
        self.root = root
        self.path = path
        self.__relative = relative

            # careful with nesting single/double quotes.
        if self.path is not None: self.path = self.path.replace('"',"'")

    def __str__( self ) :
        """
        This method returns the xPath of the linked-to instance.

        :returns:       A python str.
        """

        return( "%s" % self.path )

    def __deepcopy__(self, memo = {}):
        """
        This method returns a copy of *self*.
        """

        if self.path is None: self.updateXPath()
        theCopy = self.__class__( link = self.__link, root = self.root, path = self.path, label = self.label, relative = self.__relative )

        return theCopy

    @property
    def link(self):
        """
        This method returns a reference to the lined-to instance. If the reference is None then :py:func:`updateLink` is 
        called before the reference is returned.  Also see :py:func:`linkWithoutUpdating`.

        :returns:           The linked-to instance.
        """

        if self.__link is None:
            self.updateLink()

        return self.__link

    @link.setter
    def link(self, instance):
        """
        Thie method set's the linked-to instance to *instance*.

        :param instance:    The new linked-to instance.
        """

        self.__link = instance

    @property
    def linkWithoutUpdating(self):
        """
        This method returns a reference to the lined-to instance. If reference is None then None is returned.
        Also see :py:func:`link`.

        :returns:           The linked-to instance.
        """

        return self.__link

    @property
    def relative( self ) :
        """
        This method returns the value of the relative member.

        :returns:           A boolean.
        """

        return( self.__relative )

    def build_href( self, **kwargs ) :
        """
        This member builds the xPath and returns it (using "formatVersion" if present in *kwargs*). 
        The xPath will include the root if the root is not None.

        :param kwargs:      A dictionary of extra arguments that controls how *self* creates the href.

        :returns:           A python str.
        """

        if( self.__link is None ) :           # Should only happen when linking is to another protare's data.
            XPath = self.path
            if( self.root is not None ) : XPath = '%s#%s' % ( self.root, XPath )
        else :
            formatVersion = kwargs.get( 'formatVersion' )
            if self.__relative and hasattr(self, 'toRelativeXLink'):
                XPath = self.toRelativeXLink( self.__link, formatVersion = formatVersion )
            else:
                XPath = self.__link.toXLink( formatVersion = formatVersion )
                if( self.root is not None ) : XPath = '%s#%s' % ( self.root, XPath )

        return( XPath )

    def updateLink(self):
        """
        This method uses the xPath to reset *self*'s link.
        """

        if self.path[0] == '.':
            self.link = self.follow(self)
        elif not self.root:
            self.link = self.follow(self.rootAncestor)
        else:
            # link points to an external resource. Check in externalFiles section:
            from fudge import suites as suitesModule
            externalFiles = self.findAttributeInAncestry(suitesModule.ExternalFiles.moniker)
            if externalFiles is None:
                print("WARNING: unable to locate external resource '%s#%s': no externalFiles section found!" %
                      (self.root, self.path))
                return

            externalFile = externalFiles[self.root[1:]] # skip '$' character from start of self.root
            if externalFile.instance is None:
                print("WARNING: link '%s#%s' points to an external file that has not been loaded yet." %
                      (self.root, self.path))
                return
            self.link = self.follow(externalFile.instance)

    def updateXPath( self ):
        """This method sets *self*'s path member to the linked-to instance."""

        if self.__relative and hasattr(self, 'toRelativeXLink'):
            self.path = self.toRelativeXLink( self.link )
        elif hasattr(self.link, 'toXLink'):
            self.path = self.link.toXLink()

    def follow( self, startNode ):
        """
        This method returns the instance reference by *self*'s xPath.

        :param startNode: instance corresponding to the beginning of self.href (must inherit from ancestry)

        :return: instance pointed to by self.href
        """

        return startNode.followXPath( self.path )

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        attributeStr = self.attributesToXMLAttributeStr()
        return ['%s<%s%s href="%s"/>' % (indent, self.moniker, attributeStr, self.build_href(**kwargs))]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.
        The resulting link points to None, and must be resolved by the calling function.
        Link attributes can be converted to desired type by passing a 'typeConversion' dictionary in linkData.

        :param cls:         Class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :returns:           An instance of *cls* representing *node*.
        """

        xPath.append( node.tag )
        path = node.get('href')
        if '#' in path: root, path = path.split('#')
        else: root = None
        relative = not path.startswith('/')
        # derived classes may add new (non-xlink) attributes:
        kwargs = dict( [v for v in list( node.items( ) ) if not v[0].startswith( 'href' )] )
        for key in kwargs:
            if key in linkData.get('typeConversion',{}):
                kwargs[key] = linkData['typeConversion'][key](kwargs[key])
        result = cls(link=None, root=root, path=path, relative=relative, **kwargs)

        if( 'unresolvedLinks' in linkData ) : linkData['unresolvedLinks'].append( result )

        xPath.pop()

        return result

class Link2(ancestryModule.AncestryIO):
    """
    This class is like :py:class:`link` but adds members *moniker*, *keyName, and *keyValue*.
    """

    moniker = 'link'

    def __init__( self, moniker, instance = None, root = None, href = None, relative = False, keyName = None, keyValue = None, label = None ) :

        self.__moniker = moniker
        ancestryModule.AncestryIO.__init__(self)

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
        """
        This method returns the xPath of the linked-to instance.

        :returns:       A python str.
        """

        return( "%s" % self.__href )

    def __deepcopy__(self, memo = {}):
        """
        This method returns a copy of *self*.
        """

        if self.path is None:
            self.updateXPath()
        theCopy = self.__class__(moniker=self.moniker, instance=self.instance, root=self.root,
                                 href=self.href, relative=self.relative, keyName=self.keyName,
                                 label=self.label )

        return theCopy

    @property
    def moniker( self ) :
        """
        This method returns *self*'s moniker.

        :returns:       A python str.
        """

        return( self.__moniker )

    @property
    def instance(self) :
        """
        This method returns a reference to the lined-to instance. If the reference is None then :py:func:`updateLink` is 
        called before the reference is returned.  Also see :py:func:`linkWithoutUpdating`.

        :returns:           The linked-to instance.
        """

        if self.__instance is None:
            self.updateLink()

        return self.__instance

    @instance.setter
    def instance(self, value):

        self.__instance = value

    @property
    def link(self):
        """
        This method returns a reference to the lined-to instance. If the reference is None then :py:func:`updateLink` is 
        called before the reference is returned.  Also see :py:func:`linkWithoutUpdating`.

        :returns:           The linked-to instance.
        """

        return self.instance

    @link.setter
    def link(self, value):
        """
        Thie method set's the linked-to instance to *value*.

        :param value:       The new linked-to instance.
        """

        self.instance = value

    @property
    def keyName( self ) :
        """
        This method returns *self*'s keyName.

        :returns:           A python str.
        """

        return( self.__keyName )

    @property
    def keyValue( self ) :
        """
        This method returns *self*'s keyValue.

        :returns:           A python str.
        """

        return( self.__keyValue )

    @property
    def index(self):
        """
        This method returns *self*'s index.

        :returns:           A python str.
        """

        if self.__keyName == 'index': return int(self.keyValue)

        raise Exception('%s does not have an index member' % self.moniker)

    @index.setter
    def index(self, value):
        """
        Thie method set's *self*'s index to *value*.

        :param value:       The new index for *self*.
        """

        if self.__keyName == 'index':
            self.__keyValue = int(value)
            return

        raise Exception('%s does not have an index member' % self.moniker)

    @property
    def label( self ) :
        """
        This method returns *self*'s label.

        :returns:           A python str.
        """

        return( self.__label )

    @property
    def root( self ) :
        """
        This method returns *self*'s root.

        :returns:           A python str.
        """

        return( self.__root )

    @property
    def href( self ) :
        """
        This method returns *self*'s href.

        :returns:           A python str.
        """

        return( self.__href )

    @property
    def path( self ) :
        """
        This method returns *self*'s path.

        :returns:           A python str.
        """

        return( self.__href )

    @property
    def relative( self ) :
        """
        This method returns the value of the relative member.

        :returns:           A boolean.
        """

        return( self.__relative )

    def build_href( self, **kwargs ) :
        """
        This member builds the xPath and returns it (using "formatVersion" if present in *kwargs*). 
        The xPath will include the root if the root is not None.

        :param kwargs:      A dictionary of extra arguments that controls how *self* creates the href.

        :returns:           A python str.
        """

        if self.__instance is None:                 # Should only happen when linking to another protare's data.
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
        """
        Thie method updates *self*'s *keyName* and *keyValue* members. This is mainly for use by the constructor :py:func:`__init__`.

        :param keyName:     The new keyName.
        :param keyValue:    The new keyValue.
        """

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

    def updateLink(self):
        """
        This method uses the xPath to reset *self*'s link.
        """

        if self.path[0] == '.':
            self.link = self.follow(self)
        else:
            self.link = self.follow(self.rootAncestor)

    def updateXPath( self ) :
        """This method sets *self*'s path member to the linked-to instance."""

        if( self.__relative ) :
            self.__href = self.toRelativeXLink( self.__instance )
        else :
            if( self.__instance is not None ) :         # FIXME
                self.__href = self.__instance.toXLink( )

    def follow( self, startNode ) :
        """
        This method returns the instance reference by *self*'s xPath.

        :param startNode: instance corresponding to the beginning of self.href (must inherit from ancestry)

        :return: instance pointed to by self.href
        """

        return( startNode.followXPath( self.__href ) )

    def toXML_strList( self, indent = '', **kwargs ):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        attributes = ''
        if self.__keyName is not None and self.__keyValue is not None: attributes += ' %s="%s"' % ( self.__keyName, self.__keyValue )
        if self.__label is not None: attributes += ' label="%s"' % self.__label

        return [ '%s<%s%s href="%s"/>' % ( indent, self.moniker, attributes, self.build_href( **kwargs ) ) ]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.
        The resulting link points to None, and must be resolved by the calling function.
        Link attributes can be converted to desired type by passing a 'typeConversion' dictionary in linkData.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :returns:           An instance of *cls* representing *node*.
        """

        xPath.append( node.tag )

        root = None
        href = None
        keyName = None
        keyValue = None
        label = None

        for key in list( node.keys( ) ) :
            if( key == 'href' ) :
                href = node.get( key )
            elif( key == 'label' ) :
                label = node.get( key )
            else :
                if( keyName is not None ) : raise Exception( 'Too many attributes for link: %s' % key )
                keyName = key
                keyValue = node.get( key )

        if '#' in href: root, href = href.split('#')

        _link = cls( node.tag, instance = None, root = root, href = href, relative = href[0] != '/', keyName = keyName, keyValue = keyValue, label = label )

        if( 'unresolvedLinks' in linkData ) : linkData['unresolvedLinks'].append( _link )

        xPath.pop( )
        return( _link )
