# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

__metaclass__ = type

from . import values as valuesModule
from . import base as baseModule
from . import standards as standardsModule

symmetryNoneToken = 'none'
symmetryLowerToken = 'lower'
symmetryUpperToken = 'upper'
symmetryOptions = ( symmetryNoneToken, symmetryLowerToken, symmetryUpperToken )

permutationPlusToken = '+'
permutationMinusToken = '-'
permutationOptions = ( permutationPlusToken, permutationMinusToken )

storageRowToken = 'row-major'
storageColumnToken = 'column-major'

class arrayBase( baseModule.xDataCoreMembers ) :

    moniker = 'array'

    def __init__( self, shape = None, symmetry = None, storageOrder = storageRowToken, 
                offset = None, permutation = permutationPlusToken,
                index = None, label = None ) :

        if( not( isinstance( self.compression, str ) ) ) : raise TypeError( 'compression must be a string' )
        baseModule.xDataCoreMembers.__init__( self, self.moniker, index = index, label = label )

        shape = tuple( int( value ) for value in shape )
        if( len( shape ) == 0 ) : raise ValueError( 'shape must contain at least one value' )
        if( min( shape ) <= 0 ) : raise ValueError( 'illegal shape "%s": lengths must all be greater than 0' % str(shape) )
        self.__shape = shape
        if( self.dimension > 3 ) : raise Exception( 'Currently, dimension = %d > 3 not supported' % len( self ) )

        if( not( isinstance( symmetry, str ) ) ) : raise TypeError( 'symmetry must be a string' )
        if( symmetry not in symmetryOptions ) :
                raise ValueError( 'invalid symmetry = "%s"' % symmetry )
        self.__symmetry = symmetry
        if( symmetry != symmetryNoneToken ) :
            for length in shape :
                if( length != shape[0] ) : raise ValueError( 'a symmetrical array must be "square": shape = %s' % shape )

        if( permutation != permutationPlusToken ) : raise TypeError( 'currently, only "%s" permutation is supported' % permutationPlusToken )
        self.__permutation = permutationPlusToken

        if( not( isinstance( storageOrder, str ) ) ) : raise TypeError( 'storageOrder must be a string' )
        if( storageOrder not in [ storageRowToken, storageColumnToken ] ) :
                raise ValueError( 'invalid storageOrder = "%s"' % storageOrder )
        self.__storageOrder = storageOrder

        if( offset is not None ) :
            offset = [ int( value ) for value in offset ]
            if( len( offset ) != len( shape ) ) : raise ValueError( 'offset must contain one value for each dimension' )
            if( min( offset ) < 0 ) : raise ValueError( 'offsets must be non-negative: %s' % offset )
        self.__offset = offset

    def __len__( self ) :

        return( self.dimension )

    @property
    def dimension( self ) :

        return( len( self.__shape ) )

    @property
    def shape( self ) :

        return( self.__shape )

    @property
    def size( self ) :

        size = 1
        for length in self.__shape : size *= length
        return( size )

    @property
    def symmetry( self ) :

        return( self.__symmetry )

    @property
    def permutation( self ) :

        return( self.__permutation )

    @property
    def compression( self ) :

        return( self.compression )

    @property
    def storageOrder( self ) :

        return( self.__storageOrder )

    @property
    def offset( self ) :

        return( self.__offset )

    def offsetScaleValues( self, offset, scale ):
        """Modify every element in the array: multiply by scale and add offset."""

        self.values.offsetScaleValues( offset, scale )

    def attributesToXMLAttributeStr( self ) :

        attributeStr = ' shape="%s"' % ','.join( [ "%d" % length for length in self.shape ] )
        if( self.compression != full.compression ) : attributeStr += ' compression="%s"' % self.compression
        if( self.symmetry != symmetryNoneToken ) : attributeStr += ' symmetry="%s"' % self.symmetry
        if( self.permutation != permutationPlusToken ) : attributeStr += ' permutation="%s"' % self.permutation
        if( self.offset is not None ) : attributeStr += ' offset="%s"' % ','.join( [ "%d" % offset for offset in self.offset ] )
        if( self.storageOrder != storageRowToken ) : attributeStr += ' storageOrder="%s"' % self.storageOrder
        attributeStr += baseModule.xDataCoreMembers.attributesToXMLAttributeStr( self )
        return( attributeStr )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    @classmethod
    def parseXMLNode( cls, xDataElement, xPath, linkData, **kwargs ) :

        xPath.append( xDataElement.tag )

        attributes = arrayBase.parseXMLNodeAttributes( xDataElement )
        compression = attributes.pop( 'compression' )
        numberOfValues = { full.compression : [ 1 ], diagonal.compression : [ 1, 2 ], flattened.compression : [ 3 ],
                embedded.compression : [ -1 ] }[compression]
        if( ( numberOfValues[0] != -1 ) and ( len( xDataElement ) not in numberOfValues ) ) :
                raise Exception( '%s array expects %s sub-elements: got %d' % ( cls.compression, numberOfValues, len( xDataElement ) ) )
        shape = attributes.pop( 'shape' )

        valuesDict = {}
        if( compression != embedded.compression ) :
            values = [ valuesModule.values.parseXMLNode( valuesElements, xPath, linkData ) for valuesElements in xDataElement ]
            for value in values :
                label = value.label
                if( value.label is None ) : label = 'data'
                valuesDict[label] = value

        if( compression == full.compression ) :
            array1 = full( shape, valuesDict['data'], **attributes )
        elif( compression == diagonal.compression ) :
            if( 'startingIndices' not in valuesDict ) : valuesDict['startingIndices'] = None
            array1 = diagonal( shape, valuesDict['data'], valuesDict['startingIndices'], **attributes )
        elif( compression == flattened.compression ) :
            array1 = flattened( shape, valuesDict['data'], valuesDict['starts'], valuesDict['lengths'], **attributes )
        elif( compression == embedded.compression ) :
            array1 = embedded( shape, **attributes )
            for subArrayElement in xDataElement :
                array2 = arrayBase.parseXMLNode( subArrayElement, xPath, linkData )
                array1.addArray( array2 )
        else :
            raise TypeError( 'Unsupported array type = "%s"' % compression )

        xPath.pop( )
        return( array1 )

    @classmethod
    def parseXMLString( cls, XMLString ) :

        from xml.etree import cElementTree

        return( cls.parseXMLNode( cElementTree.fromstring( XMLString ), xPath=[], linkData={} ) )

    @staticmethod
    def parseXMLNodeAttributes( xDataElement ) :

        attributes = {  'shape'         : ( None, str ),
                        'symmetry'      : ( symmetryNoneToken, str ),
                        'permutation'   : ( permutationPlusToken, str ),
                        'storageOrder'  : ( storageRowToken, str ),
                        'offset'        : ( None, str ),
                        'compression'   : ( None, str ),
                        'index'         : ( None, int ),
                        'label'         : ( None, str ) }
        attrs = {}
        for key, item in list( attributes.items( ) ) : attrs[key] = item[0]

        for key, item in list( xDataElement.items( ) ) :
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key][1]( item )
        if( attrs['shape'] is None ) : raise ValueError( 'shape attribute is missing from array' )
        attrs['shape'] = attrs['shape'].split( ',' )
        if( attrs['offset'] is not None ) : attrs['offset'] = attrs['offset'].split( ',' )
        if( attrs['compression'] is None ) : attrs['compression'] = full.compression
        return( attrs )

    @staticmethod
    def indicesToFlatIndices( shape, indices ) :

        flatIndices = []
        for index in indices :
            if( len( index ) != len( shape ) ) :
                raise Exception( 'len( index ) = %d != len( shape ) = %d' % ( len( index ), len( shape ) ) )
            length, flatIndex = 1, 0
            for i1, i2 in enumerate( index ) :
                flatIndex *= length
                flatIndex += i2
                length = shape[i1]
            flatIndices.append( flatIndex )
        return( flatIndices )

    @staticmethod
    def flatIndicesToIndices( shape, flatIndices ) :

        indices = []
        products = [ 1 ]
        for s1 in shape : products.append( s1 * products[-1] )
        del products[-1]
        products.reverse( )
        for i1 in flatIndices :
            index = []
            for product in products :
                i2, i1 = divmod( i1, product )
                index.append( i2 )
            indices.append( index )
        return( indices )

class full( arrayBase ) :

    compression = 'full'
    ancestryMembers = baseModule.xDataCoreMembers.ancestryMembers + ( 'values', )

    def __init__( self, shape = None, data = None, symmetry = symmetryNoneToken, storageOrder = storageRowToken, 
                offset = None, permutation = permutationPlusToken,
                index = None, label = None ) :

        arrayBase.__init__( self, shape, symmetry = symmetry, storageOrder = storageOrder, 
                offset = offset, permutation = permutation,
                index = index, label = label )

        if( not( isinstance( data, valuesModule.values ) ) ) : data = valuesModule.values( data )

        if( symmetry == symmetryNoneToken ) :
            size = self.size
        else :                  # Will be a 'square' array. Checked in arrayBase.
            length, size = self.shape[0], 1
            for i1 in range( self.dimension ) : size = size * ( length + i1 ) / ( i1 + 1 )
        if( size != len( data ) ) : raise ValueError( 'shape requires %d values while data has %d values' % ( size, len( data ) ) )

        self.values = data
        self.values.setAncestor( self )

    def constructArray( self ) :

        import numpy

        if( self.symmetry == symmetryNoneToken ) :
            array1 = numpy.array( [ value for value in self.values ] )

        elif len(self.shape) == 2:
            array1 = numpy.zeros( self.shape )
            if( self.symmetry == symmetryLowerToken ) :
                array1[numpy.tril_indices(self.shape[0])] = list(self.values)
                il = numpy.tril_indices(self.shape[0], -1)
                iu = (il[1],il[0])
                array1[iu] = array1[il]
            elif( self.symmetry == symmetryUpperToken ) :
                array1[numpy.triu_indices(self.shape[0])] = list(self.values)
                iu = numpy.triu_indices(self.shape[0], 1)
                il = (iu[1],iu[0])
                array1[il] = array1[iu]

        else :
            import itertools

            array1 = numpy.zeros( self.size )
            dimension = self.dimension
            length = self.shape[0]
            indexRange = range( len( self ) )
            indices = dimension * [ 0 ]
            indexChange = dimension - 1
            mode = ( ( self.symmetry == symmetryUpperToken ) and ( self.storageOrder == storageRowToken ) ) or \
                     ( self.symmetry == symmetryLowerToken ) and ( self.storageOrder == storageColumnToken )
            for value in self.values :
                permutations = itertools.permutations( indices )
                for permutation in permutations :
                    index = permutation[0]
                    for p1 in permutation[1:] : index = length * index + p1
                    array1[index] = value
                if( mode ) :
                    for i1 in indexRange :
                        indices[i1] += 1
                        if( indices[i1] < length ) : break
                    for i2 in indexRange :
                        if( i1 == i2 ) : break
                        indices[i2] = indices[i1]
                else :
                    indexChange += 1
                    if( indexChange == dimension ) :
                        indexChange -= 1
                        value = indices[indexChange]
                        for i1 in indexRange :
                            if( i1 == ( dimension - 1 ) ) : break
                            if( indices[indexChange-1] > value ) : break
                            indices[indexChange] = 0
                            indexChange -= 1
                    indices[indexChange] += 1

        order = { storageRowToken : 'C', storageColumnToken : 'F' }[self.storageOrder]
        return( array1.reshape( self.shape, order = order ) )

    def copy( self ) :

        return( full( self.shape, self.values.copy( ), 
                symmetry = self.symmetry, storageOrder = self.storageOrder, 
                offset = self.offset, permutation = self.permutation,
                index = self.index, label = self.label ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributesStr = self.attributesToXMLAttributeStr( )
        XMLList = [ '%s<%s%s>' % ( indent, self.moniker, attributesStr ) ]
        XMLList += self.values.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

class diagonal( arrayBase ) :

    compression = 'diagonal'
    ancestryMembers = baseModule.xDataCoreMembers.ancestryMembers + ( 'values', )

    def __init__( self, shape = None, data = None, startingIndices = None, symmetry = symmetryNoneToken, storageOrder = storageRowToken, 
                offset = None, permutation = permutationPlusToken,
                index = None, label = None ) :

        arrayBase.__init__( self, shape, symmetry, storageOrder = storageOrder,
                offset = offset, permutation = permutation,
                index = index, label = label )

        dimension = self.dimension
        if( not( isinstance( data, valuesModule.values ) ) ) : data = valuesModule.values( data )
        if( startingIndices is None ) :
            self.startingIndicesOriginal = None
            startingIndices = dimension * [ 0 ]
        else :
            if( not( isinstance( startingIndices, valuesModule.values ) ) ) :
                startingIndices = valuesModule.values( startingIndices, valueType = standardsModule.types.integer32Token )
            if( startingIndices.valueType not in [ standardsModule.types.integer32Token ] ) : raise TypeError( 'startingIndices must be a list of integers' )
            self.startingIndicesOriginal = startingIndices
            self.startingIndicesOriginal.label = 'startingIndices'

        if( ( len( startingIndices ) == 0 ) or ( ( len( startingIndices ) % dimension ) != 0 ) ) :
            raise ValueError( 'lenght of startingIndices = %d must be a multiple of dimension = %d' %
                    ( len( startingIndices ), dimension ) )
        startingIndices = [ value for value in startingIndices ]
        if( min( startingIndices ) < 0 ) : raise ValueError( 'negative starting index not allowed' )
        self.startingIndices = []
        size = 0
        while( len( startingIndices ) > 0 ) :
            self.startingIndices.append( startingIndices[:dimension] )
            startingIndices = startingIndices[dimension:]
            startingIndex = self.startingIndices[-1]
            offset = self.shape[0]
            for i1, length in enumerate( self.shape ) :
                offset = min( offset, length - startingIndex[i1] )
            if( offset < 0 ) : raise ValueError( 'starting index must be less than length: %s' % startingIndex )
            size += offset
        if( size != len( data ) ) : raise ValueError( 'shape requires %d values while data has %d values' % ( size, len( data ) ) )

        self.values = data
        self.values.setAncestor( self )

    def constructArray( self ) :

        import numpy
        import itertools

        valuesIndex = 0
        array1 = numpy.zeros( self.size )
        shape = self.shape
        range1 = range( self.dimension )
        for startIndex in self.startingIndices :
            si = [ index for index in startIndex ]
            moreToDo = True
            while( moreToDo ) :
                value = self.values[valuesIndex]
                valuesIndex += 1
                if( self.symmetry == symmetryNoneToken ) :
                    permutations = [ si ]
                else :
                    permutations = itertools.permutations( si )
                for permutation in permutations :
                    scale, flatIndex = 1, 0
                    for i1, index in enumerate( permutation ) :
                        flatIndex = flatIndex * scale + index
                        scale = shape[i1]
                    array1[flatIndex] = value
                for i1 in range1 :
                    si[i1] += 1
                    if( si[i1] >= shape[i1] ) : moreToDo = False

        order = { storageRowToken : 'C', storageColumnToken : 'F' }[self.storageOrder]
        return( array1.reshape( self.shape, order = order ) )

    def copy( self ) :

        startingIndicesOriginal = self.startingIndicesOriginal
        if( startingIndicesOriginal is not None ) : startingIndicesOriginal = self.startingIndicesOriginal.copy( )
        return( diagonal( self.shape, self.values.copy( ), startingIndicesOriginal, 
                symmetry = self.symmetry, storageOrder = self.storageOrder,
                offset = self.offset, permutation = self.permutation,
                index = self.index, label = self.label ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributesStr = self.attributesToXMLAttributeStr( )
        XMLList = [ '%s<%s%s>' % ( indent, self.moniker, attributesStr ) ]
        if( self.startingIndicesOriginal is not None ) :
            XMLList += self.startingIndicesOriginal.toXMLList( indent2, **kwargs )
        XMLList += self.values.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

class flattened( arrayBase ) :

    compression = 'flattened'
    ancestryMembers = baseModule.xDataCoreMembers.ancestryMembers + ( 'values', 'lengths', 'starts' )

    def __init__( self, shape = None, data = None, starts = None, lengths = None, symmetry = symmetryNoneToken, 
                storageOrder = storageRowToken,
                offset = None, permutation = permutationPlusToken,
                index = None, label = None, 
                dataToString = None ) :

        arrayBase.__init__( self, shape, symmetry, storageOrder = storageOrder,
                offset = offset, permutation = permutation,
                index = index, label = label )

        self.dataToString = dataToString

        if( not( isinstance( data, valuesModule.values ) ) ) : data = valuesModule.values( data )
        if( not( isinstance( starts, valuesModule.values ) ) ) : 
                starts = valuesModule.values( starts, valueType = standardsModule.types.integer32Token )
        if( not( isinstance( lengths, valuesModule.values ) ) ) :
                lengths = valuesModule.values( lengths, valueType = standardsModule.types.integer32Token )

        if( len( starts ) != len( lengths ) ) : raise ValueError( 'length of starts = %d must equal length of lengths = %d' %
                ( len( starts ), len( lengths ) ) )

        size = len( data )
        length = 0
        for i1 in lengths : length += i1
        if( size != length ) : raise ValueError( 'number of data = %d and sum of length = %d differ' % ( size, length ) )

        indexPriorEnd = -1
        for i1, start in enumerate( starts ) :
            if( start < 0 ) : raise ValueError( 'negative start (=%d) not allowed' % start )
            if( start < indexPriorEnd ) : raise ValueError( 'data overlap: prior index end = %d current start = %d' % ( indexPriorEnd, start ) )
            length = lengths[i1]
            if( length < 0 ) : raise ValueError( 'negative length (=%d) not allowed' % length )
            indexPriorEnd = start + length
            if( indexPriorEnd > self.size ) :
                    raise ValueError( 'data beyond array boundary: indexPriorEnd = %d, size = %d', ( indexPriorEnd, self.size ) )

        self.starts = starts
        starts.label = 'starts'
        self.starts.setAncestor( self )

        self.lengths = lengths
        lengths.label = 'lengths'
        self.lengths.setAncestor( self )

        self.values = data
        self.values.setAncestor( self )

    def constructArray( self ) :

        import numpy

        index = 0
        array1 = numpy.zeros( self.size )
        for i1, start in enumerate( self.starts ) :
            length = self.lengths[i1]
            for i2 in range( length ) :
                array1[start+i2] = self.values[index]
                index += 1

        order = { storageRowToken : 'C', storageColumnToken : 'F' }[self.storageOrder]
        array1 = array1.reshape( self.shape, order = order )
        if self.symmetry == symmetryLowerToken:
            array1 = numpy.tril(array1) + numpy.tril(array1, -1).T
        elif self.symmetry == symmetryUpperToken:
            array1 = numpy.triu(array1) + numpy.triu(array1, -1).T
        return array1

    def copy( self ) :

        return( flattened( self.shape, self.values.copy( ), self.starts.copy( ), self.lengths.copy( ), 
                symmetry = self.symmetry, storageOrder = self.storageOrder,
                offset = self.offset, permutation = self.permutation,
                index = self.index, label = self.label ) )

    @staticmethod
    def fromNumpyArray( array, symmetry = symmetryNoneToken, nzeroes = 4 ):
        """
        Generate a sparse flattened array that represents an arbitrary numpy array.
        Only supports 'full' or 'lower-symmetric' matrices, with row-major data storage.
        :param array: input numpy array
        :param symmetry: allowed values are 'none' or 'lower'
        :param nzeroes: how many zeroes to allow before adding a new 'start' and 'length' 
        :return: 
        """

        starts, lengths, sparseData = [], [], []

        def helper( data, offset = 0 ):
            idx = 0
            end = len(data)
            while idx < end:
                if data[idx] != 0:
                    stop = idx+1
                    while stop < end:
                        if data[stop] != 0:
                            stop += 1
                        elif any(data[stop:stop + nzeroes]):
                            for i in range(nzeroes):
                                if stop + i < end and data[stop + i] != 0: stop += i
                        else:
                            break
                    starts.append( idx + offset )
                    lengths.append( stop - idx )
                    sparseData.extend( data[idx:stop] )
                    idx = stop
                idx += 1

        if symmetry == symmetryNoneToken:
            helper( array.flatten() )
        elif symmetry == symmetryLowerToken:
            rows, cols = array.shape
            for row in range(rows):
                dat = array[row][:row+1]
                helper( dat, offset = row * cols )
        else:
            raise NotImplementedError("Symmetry = '%s'" % symmetry)

        return flattened(shape=array.shape, data=sparseData, starts=starts, lengths=lengths, symmetry=symmetry)


    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributesStr = self.attributesToXMLAttributeStr( )
        XMLList = [ '%s<%s%s>' % ( indent, self.moniker, attributesStr ) ]
        XMLList += self.starts.toXMLList( indent2, **kwargs )
        XMLList += self.lengths.toXMLList( indent2, **kwargs )
        XMLList += self.values.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

class embedded( arrayBase ) :

    compression = 'embedded'

    def __init__( self, shape = None, symmetry = symmetryNoneToken, storageOrder = storageRowToken,
                offset = None, permutation = permutationPlusToken,
                index = None, label = None ) :

        arrayBase.__init__( self, shape, symmetry, storageOrder = storageOrder, 
                offset = offset, permutation = permutation,
                index = index, label = label )

        self.arrays = []

    def addArray( self, array ) :

        if( not( isinstance( array, arrayBase ) ) ) : raise TypeError( 'variable not an array instance' )
        if( self.dimension < array.dimension ) : raise ValueError( 'cannot embedded array into a smaller dimensional array: %s %s' %
                ( self.dimension, array.dimension ) )

        shape = list( reversed( array.shape ) )
        if( array.offset is None ) : raise TypeError( 'embedded array must have offset defined' )
        offsets = list( reversed( array.offset ) )
        shapeOfParent = list( reversed( self.shape ) )
        for i1, offset in enumerate( offsets ) :
            if( ( offset + shape[i1] ) > shapeOfParent[i1] ) :
                    raise ValueError( 'child array outside of parent: %s %s %s' % ( self.shape, array.shape, array.offset ) )

        self.arrays.append( array )

    def constructArray( self ) :

        import numpy

        order = { storageRowToken : 'C', storageColumnToken : 'F' }[self.storageOrder]
        array1 = numpy.zeros( self.shape, order = order )

        for array in self.arrays :
            array2 = array.constructArray( )
            slice1 = [ slice( offset, offset + array2.shape[i1] ) for i1, offset in enumerate( array.offset ) ]
            array1[slice1] = array2

        return( array1 )

    def copy( self ) :

        array1 = embedded( self.shape, 
                symmetry = self.symmetry, storageOrder = self.storageOrder,
                offset = self.offset, permutation = self.permutation,
                index = self.index, label = self.label )
        for array in self.arrays : array1.addArray( array.copy( ) )
        return array1

    def offsetScaleValues( self, offset, scale ):
        """Modify every sub-array: multiply by scale and add offset."""

        for subarray in self.arrays: subarray.offsetScaleValues( offset, scale )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributesStr = self.attributesToXMLAttributeStr( )
        XMLList = [ '%s<%s%s>' % ( indent, self.moniker, attributesStr ) ]
        for array in self.arrays : XMLList += array.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )
