# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import sys

from LUPY import enums as enumsModule
from LUPY import ancestry as ancestryModule

def isString( string ) :
    """
    Executes a TypeError raise if string is not a valid python2 str or unicode instance, or if is not a valid python3 str instance.
    Returns the string argument unchanged."""

    if( sys.version_info.major == 2 ) :
        if( not( isinstance( string, ( str, unicode ) ) ) ) : raise TypeError( 'Text must be a str or unicode.' )
    else :
        if( not( isinstance( string, str ) ) ) : raise TypeError( 'Text must be a str' )

    return( string )

class Encoding(enumsModule.Enum):

    ascii = enumsModule.auto()
    utf8 = enumsModule.auto()

class Markup(enumsModule.Enum):

    none = enumsModule.auto()
    xml = enumsModule.auto()
    html = enumsModule.auto()
    latex = enumsModule.auto()

class Text(ancestryModule.AncestryIO):

    moniker = 'text'
    # FIX-ME Should add allowed encodings as an argument with default of all. Like encodings = Encoding.allowed. Also should check in body setter.

    def __init__( self, text = None, encoding = Encoding.ascii, markup = Markup.none, label = None ) :

        ancestryModule.AncestryIO.__init__(self)

        self.encoding = Encoding.checkEnumOrString(encoding)
        self.markup = Markup.checkEnumOrString(markup)

        if( label is not None ) :
            if( not( isinstance( label, str ) ) ) : raise TypeError( 'Invalid label instance.' )
        self.__label = label

        self.body = text

    def __len__( self ) :

        return( len( self.__body ) )

    def __getitem__( self, index ) :

        return( self.__body[index] )

    @property
    def label( self ) :

        return( self.__label )

    @property
    def encoding( self ) :

        return( self.__encoding )

    @encoding.setter
    def encoding( self, value ) :

        self.__encoding = Encoding.checkEnumOrString(value)

    @property
    def markup( self ) :

        return( self.__markup )

    @markup.setter
    def markup( self, value ) :

        self.__markup = Markup.checkEnumOrString(value)

    @property
    def body( self ) :

        return( self.__body )

    @body.setter
    def body(self, text):
        """Set the text to *text*, over riding the current text."""

        self.__filled = False
        self.__body = ''

        text2 = text
        if text is None: text2 = ''

        isString(text2)

        if self.encoding == Encoding.ascii:                                 # If encoding is ascii verify text is.
            try:
                text2.encode('ascii')
            except UnicodeEncodeError as ex:
                print("    ERROR: encountered non-ascii character '%s' at index %d of text" %
                      (ex.object[ex.start], ex.start))
                print("    Bad character in context:\n", ex.object[ex.start-20:ex.end+20], "\n")
                raise ex

        if text is not None: self.__filled = True
        self.__body = text2

    @property
    def filled( self ) :

        return( self.__filled )

    def copy( self ) :

            # No need to copy text has it is immutable.
        return( Text( self.body, self.encoding, self.markup ) )

    __copy__ = copy

    def numberOfLines( self ) :

        _numberOfLines = self.__body.count( '\n' )
        if( len( self.__body ) > 0 ) :
            if( self.__body[-1] != '\n' ) : _numberOfLines += 1

        return( _numberOfLines )

    def XML_extraAttributes( self, **kwargs ) :

        return( '' )

    def bodyToXML_CDATA( self, **kwargs ) :

        return( '<![CDATA[%s]]>' % self.__body )

    def toXML_strList(self, indent = '', **kwargs):

        kwargs['addExtraAttributes'] = self.__filled
        extraAttributes = self.XML_extraAttributes(**kwargs)

        if not self.__filled:
            if extraAttributes == '' and not kwargs.get('showEmptyText', False):
                return []
            return ['%s<%s%s/>' % (indent, self.moniker, extraAttributes)]

        attributeStr = ''
        if self.__label is not None: attributeStr += ' label="%s"' % self.__label
        if self.__encoding != Encoding.ascii: attributeStr += ' encoding="%s"' % self.__encoding
        if self.__markup != Markup.none: attributeStr += ' markup="%s"' % self.__markup
        attributeStr += extraAttributes

        return [ '%s<%s%s>%s</%s>' % ( indent, self.moniker, attributeStr, self.bodyToXML_CDATA(**kwargs), self.moniker ) ]

    def parseNode(self, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        self.encoding = node.get('encoding', Encoding.ascii)
        self.markup = node.get('markup', Markup.none)

        text = node.find('text')
        if text is None:                                # Format >= '2.0.LLNL_4' parsing.
            if node.text is None:
                self.body = None
            elif len(node.text) == 0:
                self.body = None
            else:
                self.body = node.text
        else:                                           # Pre format '2.0.LLNL_4' parsing where the actual text node is a child node named 'text'.
            self.parseNode(text, xPath, linkData, **kwargs)

        xPath.pop()

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        instance = cls()
        instance.parseNode(node, xPath, linkData, **kwargs)

        return instance

def raiseIfNotString( string, name ) :

    if( not( isinstance( string, str ) ) ) : TypeError( 'Invalid %s instance.' % name )
    return( string )
