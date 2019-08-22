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

import base as baseModule

__metaclass__ = type

from pqu import PQU
from standards import types

class values( baseModule.xDataCoreMembers )  :

    moniker = 'values'

    def __init__( self, values_, valueType = types.float64Token, sep = ' ', start = 0, size = None,
                index = None, label = None ) :

        baseModule.xDataCoreMembers.__init__( self, self.moniker, index = index, label = label )

        if( valueType == types.float64Token ) :
            checker = float
        elif( valueType == types.float32Token ) :
            checker = float
        elif( valueType == types.integer32Token ) :
            checker = int
        else :
            raise TypeError( 'Invalid valueType = "%s"' % valueType )

        self.__values = [ checker( value ) for value in values_ ]
        self.__valueType = valueType 

        if( not( isinstance( sep, str ) ) or ( len( sep ) != 1 ) ) : raise TypeError( 'sep must be a string of length 1: sep = "%s"' % sep )
        self.__sep = sep

        self.__start = int( start )
        length = len( values_ )
        if( size is not None ) :
            size = int( size )
            if( size < ( self.__start + length ) ) : raise ValueError( 'size = %d < start + length = %d' % 
                    ( size, self.__start + length ) )
            self.__size = size
        else :
            self.__size = self.__start + length

    def __len__( self ) :

        return( len( self.__values ) )

    def __getitem__( self, index ) :

        return( self.__values[index] )

    def copy( self, untrimZeros = False, label = None ) :

        if( label is None ) : label = self.label
        if( ( untrimZeros ) and ( ( self.start != 0 ) or ( self.end != self.size ) ) ) :
            return( values( self.start * [ 0 ] + self.__values + ( self.size - self.end ) * [ 0 ], self.__valueType, sep = self.__sep, label = label ) )
        else :
            return( values( self.__values, self.__valueType, sep = self.sep, start = self.start, size = self.size, label = label ) )

    __copy__ = copy
    __deepcopy__ = __copy__

    @property
    def end( self ) :

        return( self.__start + len( self ) )

    @property
    def sep( self ) :

        return( self.__sep )

    @property
    def size( self ) :

        return( self.__size )

    @property
    def start( self ) :

        return( self.__start )

    @property
    def valueType( self ) :

        return( self.__valueType )

    def toString( self ) :

        strList = [ "%s" % value for value in self.__values ]
        return( self.sep.join( strList ) )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        valuesPerLine = kwargs.get( 'valuesPerLine', 100 )
        valueFormatter = kwargs.get( 'valueFormatter', PQU.toShortestString )
        sep = kwargs.get( 'sep', None )
        outline = kwargs.get( 'outline', False )
        if( len( self ) < 14 ) : outline = False

        if( sep is None ) : sep = self.sep

        strList, line = [], None
        attributeStr = ''

        length = len( self )
        if( self.start != 0 ) : attributeStr += ' start="%d"' % self.start
        if( self.end != self.size ) : attributeStr += ' size="%d"' % self.size

        if( not( isinstance( sep, str ) ) ) : raise TypeError( 'sep must be a of type string' )
        if( len( sep ) != 1 ) : raise TypeError( 'sep string length must be 1, its length is %d' % len( sep ) )
        if( sep != ' ' ) :
            attributeStr += ' sep = "%s"' % sep
            sep += ' '
        if( self.__valueType != types.float64Token ) : attributeStr += ' valueType="%s"' % self.__valueType
        attributeStr += baseModule.xDataCoreMembers.attributesToXMLAttributeStr( self )
        XMLList = [ '%s<%s length="%d"%s>' % ( indent, self.moniker, len( self.__values ), attributeStr ) ]
        if( outline ) :                     # Logic above guarantees more than 14 elements in self.
            line = []
            for i1 in range( 6 ) : line.append( valueFormatter( self[i1] ) )
            line.append( ' ... ' )
            for i1 in range( -6, 0 ) : line.append( valueFormatter( self[i1] ) )
            strList.append( sep.join( line ) )
        else :
            for index, value in enumerate( self.__values ) :
                if( ( index % valuesPerLine ) == 0 ) :
                    if( line is not None ) : strList.append( sep.join( line ) )
                    line = []
                line.append( valueFormatter( value ) )
            if( ( line is not None ) and ( len( line ) > 0 ) ) : strList.append( sep.join( line ) )
        if( len( strList ) > 0 ) : XMLList[-1] += strList.pop( 0 )

        for line in strList : XMLList.append( "%s%s" % ( indent2, line ) )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @staticmethod
    def parseXMLNode( element, xPath = [] ) :

        from numericalFunctions import listOfDoubles_C

        attrs = { 'sep' : ' ', 'valueType' : types.float64Token, 'start' : 0, 'size' : None, 'label' : None }
        attributes = { 'length' : int, 'sep' : str, 'valueType' : str, 'start' : int, 'size' : int, 'label' : str }
        for key, item in element.items( ) :
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key]( item )

        if( element.text is None ) : element.text = ''
# BRB: Need to check extraCharaters
        values1, extraCharaters = listOfDoubles_C.createFromString( element.text, sep = attrs['sep'] )
        length = attrs.pop( 'length', len( values1 ) )
        if( length != len( values1 ) ) : raise Exception( 'length = %d != len( values1 ) = %d' % ( length, len( values1 ) ) )

        values1 = [ value for value in values1 ]

        return( values( values1, **attrs ) )

    @staticmethod
    def parseXMLString( cls, XMLString ) :

        from xml.etree import cElementTree

        return( values.parseXMLNode( cElementTree.fromstring( XMLString ) ) )

    @staticmethod
    def valuesWithTrimmedZeros( values_, valueType = types.float64Token, sep = ' ' ) :

        length = len( values_ )
        start, end = 0, length - 1
        while( end >= start ) :
            if( values_[end] != 0 ) : break
            end -= 1
        end += 1
        if( end != 0 ) :
            for value in values_ :
                if( value != 0 ) : break
                start += 1
        if( ( start != 0 ) or ( end != length ) ) : values_ = values_[start:end]
        return( values( values_, valueType = valueType, sep = sep, start = start, size = length ) )
