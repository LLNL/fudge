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
import datetime

import xData.ancestry as ancestryModule

__metaclass__ = type

class styles( ancestryModule.ancestry ) :
    """
    Stores the list of PoPs styles that appear inside a file.
    """

    moniker = 'styles'
    ancestryMembers = ( 'styles', )

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )
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

        @param _style: style instance
        """

        if( not( isinstance( _style, style ) ) ) : raise TypeError( 'invalid style instance' )

        for __style in self :
            if( __style.label == _style.label ) : raise ValueError( 'style labeled "%s" already exists' % _style.label )

        self.__styles.append( _style )
        _style.setAncestor( self, 'label' )

    def convertUnits( self, unitMap ) :

        for _style in self : _style.convertUnits( unitMap )

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

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        if( len( self ) == 0 ) : return( [] )

        xmlStringList = ['%s<%s>' % (indent, self.moniker)]
        for _style in self : xmlStringList += _style.toXMLList( indent2, **kwargs )
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    def parseXMLNode( self, stylesElement, xPath, linkData ) :

        xPath.append( stylesElement.tag )

        classDict = {}
        for _style in ( evaluated, ) : classDict[_style.moniker] = _style
        for styleElement in stylesElement :
            _class = classDict.get( styleElement.tag, None )
            if( _class is None ) :
                raise TypeError( 'encountered unknown style "%s"' % styleElement.tag )
            self.add( _class.parseXMLNode( styleElement, xPath, linkData ) )

        xPath.pop()

class style( ancestryModule.ancestry ) :

    __metaclass__ = abc.ABCMeta

    def __init__( self, label, derivedFrom, date = None ) :

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( label, str ) ) ) : raise TypeError( 'label must be a str instance.' )
        self.__label = label

        if( date is None ) : date = str( datetime.date.today( ) )
        self.__date = date

        if( not( isinstance( derivedFrom, str ) ) ) : raise TypeError( 'label must be a str instance.' )
        self.__derivedFrom = derivedFrom

    @property
    def date( self ) :

        return( self.__date )

    @property
    def derivedFrom( self ) :

        return( self.__derivedFrom )

    @property
    def derivedFromStyle( self ) :

        for _style in self.ancestor :
            if( _style.label == self.__derivedFrom ) : return( _style )
        return( None )

    @property
    def label( self ) :

        return( self.__label )

    def findFormMatchingDerivedStyle( self, component, styleFilter = None ) :
        """
        This method searches the link of derivedFroms, starting with self's derivedFrom,
        to find a form in component matching one of the derivedFroms. If a form is found
        matching one of the derivedFroms, that form is returned. If no match is found, None is returned.
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

        _style = self
        while( _style is not None ) :
            _style = _style.derivedFromStyle
            if( ( _style is not None ) and ( isinstance( _style, cls ) ) ) : return( _style )
        return( None )

    def sibling( self, label ) :
        """Returns the sibling of self's with label."""

        return( self.ancestor[label] )

    def convertUnits( self, unitMap ) :

        pass

    def XMLCommonAttributes( self ) :

        XMLCommon = 'label="%s"' % self.label
        if( self.derivedFrom != "" ) : XMLCommon += ' derivedFrom="%s"' % self.derivedFrom
        if( self.date is not None ) : XMLCommon += ' date="%s"' % self.date
        return( XMLCommon )

    @staticmethod
    def parseXMLNodeBase( element, xPath ) :

        label = element.get( 'label' )
        xPath.append( '%s[@label="%s"]' % ( element.tag, label ) )

        derivedFrom = element.get( 'derivedFrom', "" )
        date = element.get( 'date', None )

        return( label, derivedFrom, date )

class evaluated( style ) :

    moniker = 'evaluated'

    def __init__( self, label, derivedFrom, library = '', version = '', date = None ) :

        style.__init__( self, label, derivedFrom, date = date )

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

        return( evaluated( self.label, self.derivedFrom, library = self.library, version = self.version, date = self.date ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        xmlStringList = [ '%s<%s %s library="%s" version="%s">' %
                ( indent, self.moniker, self.XMLCommonAttributes( ), self.library, self.version ) ]
        xmlStringList[-1] += '</%s>' % self.moniker
        return( xmlStringList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        label, derivedFrom, date = style.parseXMLNodeBase( element, xPath )

        library = element.get( 'library', '' )
        version = element.get( 'version', '' )

        _evaluated = evaluated( label, derivedFrom, library, version, date = date )

        xPath.pop( )
        return( _evaluated )
