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

from fudge.core.utilities import brb

import xData.ancestry as ancestryModule
from fudge.gnd import suites as suitesModule
from pqu.PQU import PQU as PQUModule

__metaclass__ = type

class styles( ancestryModule.ancestry ) :
    """
    Stores the list of nuclear data styles that appear inside a file.
    The list generally includes one 'evaluated' style, plus 0 or more derived styles (heated, grouped, etc.)
    """

    moniker = 'styles'

    def __init__( self ) :

        self.__styles = []
        ancestryModule.ancestry.__init__( self )

    def __len__( self ) :

        return( len( self.__styles ) )

    def __iter__( self ) :

        n1 = len( self )
        for i1 in range( n1 ) : yield self.__styles[i1]

    def __getitem__( self, label ) :

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a string' )
        for style in self.__styles :
            if( style.label == label ) : return( style )
        raise IndexError( 'No style labelled == %s' % label )

    def add( self, _style ) :
        """
        Append a style to the list of styles.
        """

        if( not( isinstance( _style, style ) ) ) : raise TypeError( 'invalid style instance' )
        if( isinstance( _style, evaluated ) ) :
            for __style in self :
                if( isinstance( __style, evaluated ) ) : raise Exception( 'multiple %s styles not allowed' % __style.moniker )
        for __style in self :
            if( __style.label == _style.label ) : raise Exception( 'style labeled "%s" already exists' )
        self.__styles.append( _style )
        _style.setAncestor( self, 'label' )

    def remove( self, _style ):
        """
        Remove the specified style. Accepts either a string or a style instance
        """
        index = None
        for idx,tmp in enumerate(self):
            if( ( tmp is _style ) or ( tmp.label == _style ) ) : index = idx
        if index is None:
            raise KeyError("style '%s' not found in styles" % _style)
        self.__styles.pop(index)

    def getStylesOfClass( self, cls ) :
        """
        Returns a list of all styles of class cls.
        """

        styleList = []
        for style in self :
            if( isinstance( style, cls ) ) : styleList.append( style )
        return( styleList )

    def getStyleOfClass( self, cls ) :
        """
        Returns a list of all styles of class cls.
        """

        styleList = self.getStylesOfClass( cls )
        if( len( styleList ) == 0 ) : return( None )
        if( len( styleList )  > 1 ) : raise Exception( '''multiple (%d) style's of type "%s" found''' % \
                ( len( styleList ), cls.moniker ) )
        return( styleList[0] )

    def getTempStyleNameOfClass( self, cls ):

        styleNames = [ style.label for style in self ]
        index = 0
        while( True ) :
            tmpLabel = 'tmp%d_%s' % ( index, cls.moniker )
            if( tmpLabel not in styleNames ) : return( tmpLabel )
            index += 1

    def getEvaluatedStyle( self ) :

        return( self.getStyleOfClass( evaluated ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = ['%s<%s>' % (indent, self.moniker)]
        for style in self : xmlStringList += style.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    def parseXMLNode( self, stylesElement, xPath, linkData ) :

        xPath.append( stylesElement.tag )
        classDict = {
            evaluated.moniker : evaluated,
            crossSectionReconstructed.moniker : crossSectionReconstructed,
            heated.moniker : heated,
            averageProductData.moniker : averageProductData }
        for styleElement in stylesElement :
            class_ = classDict.get( styleElement.tag, None )
            if class_ is None :
                raise TypeError( 'encountered unknown style "%s"' % styleElement.tag )

            self.add( class_.parseXMLNode( styleElement, xPath, linkData ) )
        xPath.pop()

class style( ancestryModule.ancestry ) :

    # FIXME should be abstract base class

    def __init__( self, label, date ) :

        ancestryModule.ancestry.__init__( self )
        self.__label = label
        self.__date = date
        self.__derivedStyles = derivedStyles( )

    @property
    def date( self ) :

        return( self.__date )

    @property
    def derivedStyles( self ) :

        return( self.__derivedStyles )

    @property
    def label( self ) :

        return( self.__label )

    def findFormMatchingDerivedStyles( self, components ) :
        """
        This methods checks sequentially through each style in self's list of styles.
        For each style, this method loops over each component in components. If a component
        with that style is found that component is returned. Otherwise, None is returned.
        """

        for style in self.__derivedStyles :
            for form in components :
                if( style.label == form.label ) : return( form )
        return( None )

    @staticmethod
    def requiredArgs() :

        return( 'label', 'date' )

    @classmethod
    def parseXMLNode( cls, styleElement, xPath, linkData ) :

        xPath.append( styleElement.tag )
        attrs = dict.fromkeys( cls.requiredArgs(), None )
        for attribute in attrs : attrs[attribute] = styleElement.get( attribute )
        for name in attrs :
            if( attrs[name] is None ) : raise Exception( 'required attribute "%s" is missing' % name )
        xPath.pop()
        return cls( **attrs )

class evaluated( style ) :

    moniker = 'evaluated'

    def __init__( self, label, date, temperature, library, version ) :

        style.__init__( self, label, date )
        if( not( isinstance( library, str ) ) ) : raise TypeError( 'library must be a string' )
        self.__library = library
        if( not( isinstance( version, str ) ) ) : raise TypeError( 'version must be a string' )
        self.__version = version
        self.__temperature = PQUModule( temperature )

    @property
    def library( self ) :

        return( self.__library )

    @property
    def temperature( self ) :

        return( self.__temperature )

    @property
    def version( self ) :

        return( self.__version )

    def toXMLList( self, indent = '', **kwargs ) :

        xmlStringList = [ '%s<%s label="%s" library="%s" version="%s" date="%s" temperature="%s">' % 
                ( indent, self.moniker, self.label, self.library, self.version, self.date, self.temperature ) ]
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def requiredArgs() :

        return ( 'label', 'date', 'temperature', 'library', 'version' )

class crossSectionReconstructed( style ) :

    moniker = 'crossSectionReconstructed'

    def __init__( self, label, date ) :

        style.__init__( self, label, date )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s label="%s" date="%s">' % ( indent, self.moniker, self.label, self.date ) ]
        xmlStringList += self.derivedStyles.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

class heated( style ) :

    moniker = 'heated'

    def __init__( self, label, date, temperature ) :

        style.__init__( self, label, date, temperature )
        self.__temperature = PQUModule( temperature )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s label="%s" date="%s" temperature="%s">' % ( indent, self.moniker, self.label, self.date,
                self.temperature ) ]
        xmlStringList += self.derivedStyles.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @property
    def temperature( self ) :

        return( self.__temperature )

class averageProductData( style ) :

    moniker = 'averageProductData'

    def __init__( self, label, date ) :

        style.__init__( self, label, date )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s label="%s" date="%s">' % ( indent, self.moniker, self.label, self.date ) ]
        xmlStringList += self.derivedStyles.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

class derivedStyles( suitesModule.suite ) :

    moniker = 'derivedStyles'

    def __init__( self ) :

        suitesModule.suite.__init__( self, [ style ] )

    def toXMLList( self, indent = '', **kwargs ) :

        if( len( self ) == 0 ) : return( [] )
        xmlStringList = '%s<%s>' % ( indent, self.moniker )
        xmlStringList += ' '.join( [ item.label for item in self ] )
        xmlStringList += '</%s>' % self.moniker
        return( [ xmlStringList ] )

def findEvaluated( forms ) :

    return( findStyle( evaluated, forms ) )

def findAllOfStyle( cls, forms ) :

    formsFound = []
    for form in forms :
        if( isinstance( form, cls ) ) : formsFound.append( form )
    return( formsFound )

def findStyle( cls, forms ) :

    formsFound = findAllOfStyle( cls, forms )
    if( len( formsFound )  > 1 ) : raise Exception( 'multiple styles found' )
    if( len( formsFound ) == 0 ) : raise Exception( 'no style found' )
    return( formsFound[0] )
