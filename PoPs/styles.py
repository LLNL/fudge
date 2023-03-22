# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import abc
import datetime

from LUPY import misc as LUPY_miscModule
from LUPY import ancestry as ancestryModule
from xData.Documentation import documentation as documentationModule

class Styles(ancestryModule.AncestryIO_base):
    """
    Stores the list of PoPs styles that appear inside a database.
    """
    # FIXME merge with fudge.styles class?

    moniker = 'styles'

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__(self)
        self.__styles = []

    def __contains__( self, label ) :

        for item in self :
            if( item.label == label ) : return( True )
        return( False )

    def __len__( self ) :

        return( len( self.__styles ) )

    def __iter__( self ) :

        n1 = len( self )
        for i1 in range( n1 ) : yield self.__styles[i1]

    def __getitem__( self, label ) :

        if( isinstance( label, int ) ) : return( self.__styles[label] )
        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        for _style in self.__styles :
            if( _style.label == label ) : return( _style )
        raise IndexError( 'No style labelled == "%s"' % label )

    def add( self, _style ) :
        """
        Append a style to the list of styles.

        :param _style: style instance
        """

        if( not( isinstance( _style, Style ) ) ) : raise TypeError( 'invalid style instance' )

        for __style in self :
            if( __style.label == _style.label ) : raise ValueError( 'style labeled "%s" already exists' % _style.label )

        self.__styles.append( _style )
        _style.setAncestor(self)

    def convertUnits( self, unitMap ) :
        """See database.convertUnits"""

        for _style in self : _style.convertUnits( unitMap )

    def remove( self, _style ):
        """
        Remove the specified style. Raises KeyError if specified style not found.
        :param _style: may be a string or a style instance
        """

        index = None
        for idx,tmp in enumerate(self):
            if( ( tmp is _style ) or ( tmp.label == _style ) ) : index = idx
        if index is None:
            raise KeyError("style '%s' not found in styles" % _style)
        self.__styles.pop(index)

    def getEvaluatedStyle(self):
        _evaluateds = [_style for _style in self if isinstance(_style, Evaluated)]
        if len(_evaluateds) == 0:
            raise Exception("No evaluated style found.")
        if len(_evaluateds) > 1:
            raise Exception("Multiple (%d) evaluated styles found." % len(_evaluateds))
        return _evaluateds[0]

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        if( len( self ) == 0 ) : return( [] )

        xmlStringList = ['%s<%s>' % (indent, self.moniker)]
        for _style in self : xmlStringList += _style.toXML_strList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    def parseNode(self, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        classDict = {}
        for _style in ( Evaluated, ) : classDict[_style.moniker] = _style
        for child in node:
            _class = classDict.get( child.tag, None )
            if( _class is None ) :
                raise TypeError( 'encountered unknown style "%s"' % child.tag )
            self.add(_class.parseNodeUsingClass(child, xPath, linkData, **kwargs))


        xPath.pop()

class Style(ancestryModule.AncestryIO, abc.ABC):
    """ Abstract base class for all 'style' classes in PoPs """

    keyName = 'label'

    def __init__( self, label, derivedFrom, date = None ) :
        """
        :param label: string identifying this style
        :param derivedFrom: may be a string (equal to another style's label), or None if this is an original evaluation
        :param date: optional date when style was generated. Defaults to today's date
        """

        ancestryModule.AncestryIO.__init__( self )

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a str instance.' )
        self.__label = label

        if( date is None ) : date = str( datetime.date.today( ) )
        self.__date = date

        self.derivedFrom = derivedFrom

        self.__documentation = documentationModule.Documentation()
        self.__documentation.setAncestor(self)

    @property
    def date( self ) :
        """Returns date when style was generated"""

        return( self.__date )

    @property
    def derivedFrom( self ) :
        """Returns string label of the style that self derives from"""

        return( self.__derivedFrom )

    @derivedFrom.setter
    def derivedFrom(self, derivedFrom):
        '''Sets *self*'s derivedFrom to *derivedFrom*.'''

        self.__derivedFrom = LUPY_miscModule.isString(derivedFrom, 'derivedFrom must be a string.')

    @property
    def derivedFromStyle( self ) :
        """Returns style instance that self derives from"""

        for _style in self.ancestor :
            if( _style.label == self.__derivedFrom ) : return( _style )
        return( None )

    @property
    def documentation(self):
        return self.__documentation

    @property
    def label( self ) :

        return( self.__label )

    def findFormMatchingDerivedStyle( self, component, styleFilter = None ) :
        """
        This method searches backwards through the chain of derived styles to find a matching form in
        the specified component. It starts with the style that self derives from, then the style that derives from, etc.
        The first form matching one of those styles is returned (unless styleFilter is used to screen out some styles).

        :param component: a PoPs component, e.g. 'mass', 'spin' or 'energy'.
        :param styleFilter: optional function that can be used to filter out some matches. styleFilter should
            take a style as input and return False if that style should be omitted from matches
        :return: matching form, or None if no match is found
        """

        def alwaysTrue( a_style ) : return( True )

        if( styleFilter is None ) : styleFilter = alwaysTrue

        parent = self.ancestor
        derivedFrom = self.derivedFrom
        while( derivedFrom != '' ) :
            for form in component :
                if( ( derivedFrom == form.label ) and styleFilter( parent[derivedFrom] ) ) : return( form )
            derivedFrom = parent[derivedFrom].derivedFrom
        return( None )
        
    def findDerivedFromStyle( self, cls ) :
        """
        This method searches backwards through the chain of derived styles to find an instance of the specified
        class (which should inherit from the base style class).

        :param cls: desired class
        :return: matching class instance, or None if no match found
        """

        _style = self
        while( _style is not None ) :
            _style = _style.derivedFromStyle
            if( ( _style is not None ) and ( isinstance( _style, cls ) ) ) : return( _style )
        return( None )

    def sibling( self, label ) :
        """Returns the sibling of self with label."""

        return( self.ancestor[label] )

    def convertUnits( self, unitMap ) :

        pass

    def XMLCommonAttributes( self ) :
        """Helper function for writing to XML"""

        XMLCommon = 'label="%s"' % self.label
        if( self.derivedFrom != "" ) : XMLCommon += ' derivedFrom="%s"' % self.derivedFrom
        if( self.date is not None ) : XMLCommon += ' date="%s"' % self.date
        return( XMLCommon )

    @staticmethod
    def parseBaseNodeCommons(node, xPath, **kwargs):
        """
        Parse *label*, *derivedFrom* and *date* that are common to all style nodes. Also, set xPath.
        """

        label = node.get('label')
        xPath.append('%s[@label="%s"]' % (node.tag, label))

        derivedFrom = node.get('derivedFrom', '')
        date = node.get('date', None)

        return label, derivedFrom, date

class Evaluated(Style):

    moniker = 'evaluated'

    def __init__( self, label, derivedFrom, library = '', version = '', date = None ) :
        """
        :param label: a unique string identifying this style
        :param derivedFrom: may be a string (equal to another style's label). For an original evaluation use an empty string.
        :param library: optional string identifying the library, e.g. 'ENSDF'
        :param version: optional string identifying the library version
        :param date: optional date when style was generated. Defaults to today's date
        """

        Style.__init__( self, label, derivedFrom, date = date )

        if( not( isinstance( library, str ) ) ) : raise TypeError( 'library must be a string' )
        self.__library = library

        if( not( isinstance( version, str ) ) ) : raise TypeError( 'version must be a string' )
        self.__version = version

    @property
    def library( self ) :

        return( self.__library )

    @library.setter
    def library( self, value ) :

        self.__library = value

    @property
    def version( self ) :

        return( self.__version )

    @version.setter
    def version( self, value ) :

        self.__version = value

    def copy( self ) :

        return( Evaluated( self.label, self.derivedFrom, library = self.library, version = self.version, date = self.date ) )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlStringList = [ '%s<%s %s library="%s" version="%s">' %
                ( indent, self.moniker, self.XMLCommonAttributes( ), self.library, self.version ) ]

        if kwargs.get('formatVersion') != '0.1':  # PoPs version corresponding to GNDS 1.10
            xmlStringList += self.documentation.toXML_strList( indent2, **kwargs )

        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        label, derivedFrom, date = Style.parseBaseNodeCommons(node, xPath, **kwargs)        # Also adds to xPath.

        library = node.get('library', '')
        version = node.get('version', '')

        _evaluated = cls( label, derivedFrom, library, version, date = date )

        _documentation = node.find(documentationModule.Documentation.moniker)
        if _documentation is not None:
            _evaluated.documentation.parseNode(_documentation, xPath, linkData, **kwargs)

        xPath.pop()                 # Need to pop as parseBaseNodeCommons added to it.

        return _evaluated
