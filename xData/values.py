# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__metaclass__ = type

import copy

from numericalFunctions import pointwiseXY_C as pointwiseXY_CModule
floatToShortestString = pointwiseXY_CModule.floatToShortestString

from numericalFunctions import listOfDoubles_C as listOfDoubles_CModule

from . import base as baseModule
from . import standards as standardsModule

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
        writeLength = kwargs.get( 'writeLengthAttribute', False )
        sep = kwargs.get( 'sep', None )
        outline = kwargs.get( 'outline', False )
        if( len( self ) < 14 ) : outline = False

        if( sep is None ) : sep = self.sep

        strList, line = [], None
        attributeStr = ''

        if( writeLength ) : attributeStr += ' length="%d"' % len( self.__values )
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

        HDF_opts = kwargs.get("HDF_opts")
        if HDF_opts is not None and len(self) >= HDF_opts['minLength']:

            if HDF_opts['flatten']:
                if self.valueType == standardsModule.types.integer32Token:
                    datasetName = 'iData'
                else:
                    datasetName = 'dData'
                flatData = HDF_opts[datasetName]
                startLength = len(flatData)
                flatData.extend(list(self))
                XMLList = ['%s<%s href="HDF#/%s" offset="%d" count="%d"%s/>' %
                        (indent, self.moniker, datasetName, startLength, len(self), attributeStr)]
            else:
                datasetName = "data%d" % HDF_opts['index']
                HDF_opts['index'] += 1
                HDF_opts['h5file'].create_dataset(datasetName, data=list(self), **HDF_opts['compression'])
                XMLList = ['%s<%s href="HDF#/%s"%s/>' % (indent, self.moniker, datasetName, attributeStr)]
            return XMLList

        XMLList = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
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

        attrs = { 'sep' : ' ', 'valueType' : standardsModule.types.float64Token, 'start' : 0, 'size' : None, 'label' : None }
        attributes = { 'length' : int, 'sep' : str, 'valueType' : str, 'start' : int, 'size' : int, 'label' : str }
        for key, item in list( element.items( ) ) :
            if( key in ('href','indices') ): continue   # handled below
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key]( item )

        href = element.get("href")
        if href is not None:
            HDF = linkData['HDF']
            indices = element.get("indices")
            s1,s2 = None, None
            if indices:
                s1,s2 = map(int,indices.split(','))
            datasetName = href.split('#')[-1]
            if datasetName.lstrip('/') in HDF:
                data = HDF[datasetName.lstrip('/')]
                return values( data[s1:s2], **attrs )
            else:
                data = HDF['h5File'].get( datasetName )
                return values( data[s1:s2], **attrs )

        if( element.text is None ) : element.text = ''
        if attrs['valueType'] == standardsModule.types.integer32Token:
            if( attrs['sep'].isspace( ) ) :
                values1 = list( map( int, element.text.split( ) ) )
            else :
                values1 = list( map( int, element.text.split(attrs['sep'] ) ) )
        else:
# BRB: Need to check extraCharacters
            values1, extraCharacters = listOfDoubles_CModule.createFromString( element.text, sep = attrs['sep'] )
            values1 = [ value for value in values1 ]

        length = attrs.pop( 'length', len( values1 ) )
        if( length != len( values1 ) ) : raise Exception( 'length = %d != len( values1 ) = %d' % ( length, len( values1 ) ) )

        return( values( values1, **attrs ) )

    @staticmethod
    def parseXMLString( XMLString ) :

        from xml.etree import cElementTree

        return( values.parseXMLNode( cElementTree.fromstring( XMLString ), xPath=[], linkData={} ) )

    @staticmethod
    def sortAndThin( _values, rel_tol = 0.0, abs_tol = 0.0 ) :
        """
        This function sorts the list of _values into an ascending order and then thins the list. 
        If there are 2 or fewer points in _values, nothing is done.
        """

        __values = [ float( value ) for value in _values ]
        __values.sort( )

        length = len( __values )
        ___values = __values[:1]

        if( length < 2 ) : return( ___values )

        x1 = ___values[0]
        for i1 in range( 1, length ) :
            x2 = __values[i1]
            deltaMax = max( abs_tol, rel_tol * min( abs( x1 ), abs( x2 ) ) )
            if( x2 - x1 <= deltaMax ) : continue
            ___values.append( x2 )
            x1 = x2

        ___values[-1] = __values[-1]                # If last point of __values did not get added, replace last point with it, otherwise this does nothing.
            
        return( ___values )

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
