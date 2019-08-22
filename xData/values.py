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

__metaclass__ = type

from . import base as baseModule
from . import standards as standardsModule

from numericalFunctions import pointwiseXY_C as pointwiseXY_CModule
floatToShortestString = pointwiseXY_CModule.floatToShortestString

class values( baseModule.xDataCoreMembers )  :

    moniker = 'values'

    def __init__( self, _values, valueType = standardsModule.types.float64Token, sep = ' ', start = 0, size = None,
                index = None, label = None ) :

        baseModule.xDataCoreMembers.__init__( self, self.moniker, index = index, label = label )

        if( valueType not in ( standardsModule.types.float64Token, standardsModule.types.float32Token,
                standardsModule.types.integer32Token ) ) :
            raise TypeError( 'Invalid valueType = "%s"' % valueType )
        self.__valueType = valueType 

        if( not( isinstance( sep, str ) ) or ( len( sep ) != 1 ) ) : raise TypeError( 'sep must be a string of length 1: sep = "%s"' % sep )
        self.__sep = sep

        self.__start = int( start )
        length = len( _values )
        if( size is not None ) :
            size = int( size )
            if( size < ( self.__start + length ) ) : raise ValueError( 'size = %d < start + length = %d' % 
                    ( size, self.__start + length ) )
            self.__size = size
        else :
            self.__size = self.__start + length

        self.__values = _values

    def __len__( self ) :

        return( len( self.__values ) )

    def __getitem__( self, index ) :

        return( self.__values[index] )

    def copy( self, untrimZeros = False ) :

        import copy
        start = self.start
        size = self.size
        _values = copy.copy( self.__values )
        if( ( untrimZeros ) and ( ( self.start != 0 ) or ( self.end != self.size ) ) ) :
            start = 0
            size = None
            _values = self.start * [ 0 ] + self.__values + ( self.size - self.end ) * [ 0 ]
        return( values( _values, self.valueType, sep = self.__sep, start = start, size = size, 
                label = self.label ) )

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
    def values( self ) :

        return( self.__values )

    @values.setter
    def values( self, _values ) :

        if( self.valueType == standardsModule.types.float64Token ) :
            checker = float
        elif( self.valueType == standardsModule.types.float32Token ) :
            checker = float
        else :
            checker = int

        length = len( _values )
        self.__size = self.__start + length

        self.__values = [ checker( value ) for value in _values ]

    @property
    def valueType( self ) :

        return( self.__valueType )


    def offsetScaleValues( self, offset, scale ) :

        for i1, value in enumerate(self.__values) : self.__values[i1] = value * scale + offset

    def toString( self ) :

        strList = [ "%s" % value for value in self.__values ]
        return( self.sep.join( strList ) )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        def intValueFormatter( value, significantDigits = 0 ) :     # Ignores significantDigits.

            return( str( value ) )

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        valuesPerLine = kwargs.get( 'valuesPerLine', 100 )
        valueFormatter = kwargs.get( 'valueFormatter', floatToShortestString )
        significantDigits = kwargs.get( 'significantDigits', 15 )
        sep = kwargs.get( 'sep', None )
        outline = kwargs.get( 'outline', False )
        if( len( self ) < 14 ) : outline = False

        if( sep is None ) : sep = self.sep

        strList, line = [], None
        attributeStr = ''

        if( self.start != 0 ) : attributeStr += ' start="%d"' % self.start
        if( self.end != self.size ) : attributeStr += ' size="%d"' % self.size

        if( not( isinstance( sep, str ) ) ) : raise TypeError( 'sep must be a of type string' )
        if( len( sep ) != 1 ) : raise TypeError( 'sep string length must be 1, its length is %d' % len( sep ) )
        if( sep != ' ' ) :
            attributeStr += ' sep = "%s"' % sep
            sep += ' '
        if( self.valueType != standardsModule.types.float64Token ) : attributeStr += ' valueType="%s"' % self.valueType
        if( self.valueType == standardsModule.types.integer32Token ) : valueFormatter = intValueFormatter
        attributeStr += baseModule.xDataCoreMembers.attributesToXMLAttributeStr( self )
        XMLList = [ '%s<%s length="%d"%s>' % ( indent, self.moniker, len( self.__values ), attributeStr ) ]
        if( outline ) :                     # Logic above guarantees more than 14 elements in self.
            line = []
            for i1 in range( 6 ) : line.append( valueFormatter( self[i1], significantDigits = significantDigits ) )
            line.append( ' ... ' )
            for i1 in range( -6, 0 ) : line.append( valueFormatter( self[i1], significantDigits = significantDigits ) )
            strList.append( sep.join( line ) )
        else :
            dataToString = kwargs.get( 'dataToString', None )
            if( dataToString is None ) :
                for index, value in enumerate( self.__values ) :
                    if( ( index % valuesPerLine ) == 0 ) :
                        if( line is not None ) : strList.append( sep.join( line ) )
                        line = []
                    line.append( valueFormatter( value, significantDigits = significantDigits ) )
                if( ( line is not None ) and ( len( line ) > 0 ) ) : strList.append( sep.join( line ) )
            else :
                kwargs['valueFormatter'] = valueFormatter
                XMLList += kwargs['dataToString']( self, kwargs['dataToStringParent'], indent = indent2, **kwargs )
        if( len( strList ) > 0 ) : XMLList[-1] += strList.pop( 0 )

        for line in strList : XMLList.append( "%s%s" % ( indent2, line ) )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @staticmethod
    def parseXMLNode( element, xPath, linkData ) :

        from numericalFunctions import listOfDoubles_C

        attrs = { 'sep' : ' ', 'valueType' : standardsModule.types.float64Token, 'start' : 0, 'size' : None, 'label' : None }
        attributes = { 'length' : int, 'sep' : str, 'valueType' : str, 'start' : int, 'size' : int, 'label' : str }
        for key, item in element.items( ) :
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key]( item )

        if( element.text is None ) : element.text = ''
        if attrs['valueType'] == standardsModule.types.integer32Token:
            if attrs['sep'].isspace(): values1 = map(int, element.text.split())
            else: values1 = map(int, element.text.split(attrs['sep']))
        else:
# BRB: Need to check extraCharacters
            values1, extraCharacters = listOfDoubles_C.createFromString( element.text, sep = attrs['sep'] )
            values1 = [ value for value in values1 ]

        length = attrs.pop( 'length', len( values1 ) )
        if( length != len( values1 ) ) : raise Exception( 'length = %d != len( values1 ) = %d' % ( length, len( values1 ) ) )

        return( values( values1, **attrs ) )

    @staticmethod
    def parseXMLString( XMLString ) :

        from xml.etree import cElementTree

        return( values.parseXMLNode( cElementTree.fromstring( XMLString ), xPath=[], linkData={} ) )

    @staticmethod
    def valuesWithTrimmedZeros( _values, valueType = standardsModule.types.float64Token, sep = ' ' ) :

        length = len( _values )
        start, end = 0, length - 1
        while( end >= start ) :
            if( _values[end] != 0 ) : break
            end -= 1
        end += 1
        if( end != 0 ) :
            for value in _values :
                if( value != 0 ) : break
                start += 1
        if( ( start != 0 ) or ( end != length ) ) : _values = _values[start:end]
        return( values( _values, valueType = valueType, sep = sep, start = start, size = length ) )
