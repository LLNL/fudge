# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module containes all the classes for handling GNDS array containers.

This module contains the following classes:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | Symmetry                          | This class is an enum of the supported symmetries for an array.       |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Permutation                       | This class is an enum of the supported premutations for an array.     |
    +-----------------------------------+-----------------------------------------------------------------------+
    | ArrayBase                         | This class is the base class for the array classes.                   |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Full                              | This class represents an array with all cells containing data.        |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Diagonal                          | This class represents an array where the data are stored in the       |
    |                                   | compact GNDS diagaonal array format.                                  |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Flattened                         | This class represents an array where the data are stored in the       |
    |                                   | compact GNDS flattened array format.                                  |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Embedded                          | This class represents an array that contains other arrays.            |
    +-----------------------------------+-----------------------------------------------------------------------+
"""

import numpy

from LUPY import enums as LUPY_enumsModule

from . import enums as enumsModule
from . import values as valuesModule
from . import base as baseModule
from . import vector as vectorModule
from . import matrix as matrixModule

class Symmetry(LUPY_enumsModule.Enum):
    """
    This class is an enum of the supported symmetries for an array.
    """

    none = LUPY_enumsModule.auto()
    lower = LUPY_enumsModule.auto()
    upper = LUPY_enumsModule.auto()

class Permutation(LUPY_enumsModule.Enum):
    """
    This class is an enum of the supported premutations for an array.
    """

    plus = '+'
    minus = '-'

class ArrayBase( baseModule.XDataCoreMembers ) :
    """
    This class is the base class for the array classes.

    The following table list the primary members of this class:

    +-------------------+-----------------------------------------------------------------------------------------------+
    | Member            | Description                                                                                   |
    +===================+===============================================================================================+
    | shape             | This is the shape of the array.                                                               |
    +-------------------+-----------------------------------------------------------------------------------------------+
    | symmetry          | This is the symmtry (none, lower or upper) for the data in the array.                         |
    +-------------------+-----------------------------------------------------------------------------------------------+
    | storageOrder      | This is the storage order (row or column) for the data in the arrary.                         |
    +-------------------+-----------------------------------------------------------------------------------------------+
    | offset            | This is the offset index where the first data value in the array starts.                      |
    +-------------------+-----------------------------------------------------------------------------------------------+
    | permutation       | This is the permutation for the data in the array.                                            |
    +-------------------+-----------------------------------------------------------------------------------------------+
    | index             | This is index of the array inside of other xData function.                                    |
    +-------------------+-----------------------------------------------------------------------------------------------+
    | label             | This is the label for the array.                                                              |
    +-------------------+-----------------------------------------------------------------------------------------------+
    """

    moniker = 'array'

    def __init__( self, shape = None, symmetry = None, storageOrder=enumsModule.StorageOrder.rowMajor, 
                offset = None, permutation = Permutation.plus,
                index = None, label = None ) :
        """
        :param shape:               This is the shape of the array.
        :param symmetry:            This is the symmtry (none, lower or upper) for the data in the array.
        :param storageOrder:        This is the storage order (row or column) for the data in the arrary.
        :param offset:              This is the offset index where the first data value in the array starts.
        :param permutation:         This is the permutation for the data in the array.
        :param index:               This is index of the array inside of other xData function.
        :param label:               This is the label for the array.
        """

        if( not( isinstance( self.compression, str ) ) ) : raise TypeError( 'compression must be a string' )
        baseModule.XDataCoreMembers.__init__(self, index=index, label=label)

        shape = tuple( int( value ) for value in shape )
        if( len( shape ) == 0 ) : raise ValueError( 'shape must contain at least one value' )
        if( min( shape ) <= 0 ) : raise ValueError( 'illegal shape "%s": lengths must all be greater than 0' % str(shape) )
        self.__shape = shape
        if( self.dimension > 3 ) : raise Exception( 'Currently, dimension = %d > 3 not supported' % len( self ) )

        self.__symmetry = Symmetry.checkEnumOrString(symmetry)
        if symmetry != Symmetry.none:
            for length in shape :
                if( length != shape[0] ) : raise ValueError( 'a symmetrical array must be "square": shape = %s' % shape )

        if permutation != Permutation.plus: raise TypeError('currently, only "%s" permutation is supported' % Permutation.plus)
        self.__permutation = permutation

        self.__storageOrder = storageOrder = enumsModule.StorageOrder.checkEnumOrString(storageOrder)

        if( offset is not None ) :
            offset = [ int( value ) for value in offset ]
            if( len( offset ) != len( shape ) ) : raise ValueError( 'offset must contain one value for each dimension' )
            if( min( offset ) < 0 ) : raise ValueError( 'offsets must be non-negative: %s' % offset )
        self.__offset = offset

    def __len__( self ) :
        """
        This method returns the dimension of the array.

        :returns:           A python int.
        """

        return( self.dimension )

    @property
    def dimension( self ) :
        """
        This method returns the dimension of the array.

        :returns:           A python int.
        """

        return( len( self.__shape ) )

    @property
    def shape( self ) :
        """
        This method returns a reference to the shape member.

        :returns:            A python tuple.
        """

        return( self.__shape )

    @property
    def size( self ) :
        """
        This method returns the number of cells in a full representation of *self*.

        :returns:            A python int.
        """

        size = 1
        for length in self.__shape : size *= length
        return( size )

    @property
    def symmetry( self ) :
        """
        This method returns the symmetry of *self*.

        :returns:            An instance of :py:class:`Symmetry`.
        """

        return( self.__symmetry )

    @property
    def permutation( self ) :
        """
        This method returns the permutation of *self*.

        :returns:            An instance of :py:class:`Permutation`.
        """

        return( self.__permutation )

    @property
    def compression( self ) :
        """
        This method returns the compression of *self*.

        :returns:            A python str.
        """

        return( self.compression )

    @property
    def storageOrder( self ) :
        """
        This method returns the storageOrder of *self*.

        :returns:            An instance of :py:class:`enumsModule.StorageOrder`.
        """

        return( self.__storageOrder )

    @property
    def offset( self ) :
        """
        This method returns the offset of *self*.

        :returns:            A python tuple.
        """

        return( self.__offset )

    def offsetScaleValues( self, offset, scale ):
        """
        This method modifies every cell in the array by multiply it by *scale* and adding *offset*.

        :param offset:      The offset to apply to each cell.
        :param scale:       The multiplicative scale factor to apply to each cell.
        """

        self.values.offsetScaleValues( offset, scale )

    def constructVector(self, secondIndex=None):
        """
        This method only works for 1d and 2d arrays. For a 1d array, its returns the data as an instancce of :py:class:`vectorModule.Vector`.
        For a 2d array, it returns the data for row *secondIndex* as an instancce of :py:class:`vectorModule.Vector`.

        :param secondIndex:     Index along the second dimension from which to extract 1d-array.

        :returns:               An instance of :py:class:`vectorModule.Vector`.
        """

        if self.dimension == 2:
            assert secondIndex is not None and isinstance(secondIndex, int), 'The method constructVector requires the first argument to be an integer.'

            if self.shape[1] <= secondIndex:
                return vectorModule.Vector()

            else:
                return vectorModule.Vector(values=self.constructArray()[:, secondIndex])

        elif self.dimension == 1:
            return vectorModule.Vector(values=self.constructArray())

        else:
            raise TypeError(f'The method constructVector is only defined for 1D or 2D arrays')

    def constructMatrix(self, thirdIndex=None):
        """
        This method only works for 3d array.
        This method returns the data at dimension *thirdIndex* as a :py:class:`matrixModule.Matrix` instance.

        :param thirdIndex:      Index along the third dimension from which to extract 2d-array.

        :returns:               An instance of :py:class:`matrixModule.Matrix`.
        """

        if self.dimension != 3:
            raise TypeError(f'The method constructMatrix is only defined for 3D arrays')

        assert thirdIndex is not None and isinstance(thirdIndex, int), 'The method constructMatrix requires the first argument to be an integer.'
        if self.shape[2] <= thirdIndex:
            return matrixModule.Matrix()

        else:
            return matrixModule.Matrix(values=self.constructArray()[:, :, thirdIndex])

    def attributesToXMLAttributeStr( self ) :
        """
        This method converts all the attributes of *self* into an XML attribute string.

        :returns:       A python string.
        """

        attributeStr = ' shape="%s"' % ','.join( [ "%d" % length for length in self.shape ] )
        if self.compression != Full.compression: attributeStr += ' compression="%s"' % self.compression
        if self.symmetry != Symmetry.none: attributeStr += ' symmetry="%s"' % self.symmetry
        if self.permutation != Permutation.plus: attributeStr += ' permutation="%s"' % self.permutation
        if self.offset is not None: attributeStr += ' offset="%s"' % ','.join( [ "%d" % offset for offset in self.offset ] )
        if self.storageOrder != enumsModule.StorageOrder.rowMajor: attributeStr += ' storageOrder="%s"' % self.storageOrder
        attributeStr += baseModule.XDataCoreMembers.attributesToXMLAttributeStr( self )
        return( attributeStr )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :returns:           An instance of *cls* representing *node*.
        """

        xPath.append(node.tag)

        attributes = ArrayBase.parseNodeAttributes(node, **kwargs)
        compression = attributes.pop('compression')
        numberOfValues = { Full.compression : [ 1 ], Diagonal.compression : [ 1, 2 ], Flattened.compression : [ 3 ],
                Embedded.compression : [ -1 ] }[compression]
        if numberOfValues[0] != -1 and len(node) not in numberOfValues:
                raise Exception('%s array expects %s sub-elements: got %d' % ( cls.compression, numberOfValues, len(node) ))
        shape = attributes.pop('shape')

        valuesDict = {}
        if compression != Embedded.compression:
            values = [ valuesModule.Values.parseNodeUsingClass(valuesElements, xPath, linkData, **kwargs) for valuesElements in node ]
            for value in values :
                label = value.label
                if value.label is None: label = 'data'
                valuesDict[label] = value

        if compression == Full.compression:
            array1 = Full(shape, valuesDict['data'], **attributes)
        elif compression == Diagonal.compression:
            if 'startingIndices' not in valuesDict: valuesDict['startingIndices'] = None
            array1 = Diagonal(shape, valuesDict['data'], valuesDict['startingIndices'], **attributes)
        elif compression == Flattened.compression:
            array1 = Flattened(shape, valuesDict['data'], valuesDict['starts'], valuesDict['lengths'], **attributes)
        elif compression == Embedded.compression:
            array1 = Embedded(shape, **attributes)
            for subArrayElement in node:
                array2 = ArrayBase.parseNodeUsingClass(subArrayElement, xPath, linkData, **kwargs)
                array1.addArray(array2)
        else :
            raise TypeError('Unsupported array type = "%s"' % compression)

        xPath.pop()

        return array1

    @staticmethod
    def parseNodeAttributes(node, **kwargs):
        """
        This static method parses the attributes for an array.

        :param node:        Node whose attributes are parsed.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :returns:           A python dict.
        """

        attributes = {  'shape'         : (None, str),
                        'symmetry'      : (Symmetry.none, str),
                        'permutation'   : (Permutation.plus, str),
                        'storageOrder'  : (enumsModule.StorageOrder.rowMajor, str),
                        'offset'        : (None, str),
                        'compression'   : (None, str),
                        'index'         : (None, int),
                        'label'         : (None, str) }
        attrs = {}
        for key, item in list(attributes.items()): attrs[key] = item[0]

        for key, item in list(node.items()):
            if key not in attributes: raise TypeError('Invalid attribute "%s"' % key)
            attrs[key] = attributes[key][1](item)
        if attrs['shape'] is None: raise ValueError('shape attribute is missing from array')
        attrs['shape'] = attrs['shape'].split(',')
        if attrs['offset'] is not None: attrs['offset'] = attrs['offset'].split(',')
        if attrs['compression'] is None: attrs['compression'] = Full.compression

        return attrs

    @staticmethod
    def indicesToFlatIndices( shape, indices ) :
        """
        This static method converts each index in *indices* into a flat index.

        :param shape:       The shape of the array.
        :param indices:     The list of indices within *shape*.

        :returns:           A python list.
        """

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
        """
        This static method converts each flat index in *flatIndices* into an index in an array with shape *shape*.

        :param shape:       The shape of the array.
        :param indices:     The list of flat indices.

        :returns:           A python list.
        """

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

class Full( ArrayBase ) :
    """
    This class store the array without any compression.

    The following table list the primary members of this class that are not specified in the primary members
    table in :py:class:`ArrayBase`:

    +-------------------+-----------------------------------------------------------------------------------------------+
    | Member            | Description                                                                                   |
    +===================+===============================================================================================+
    | data              | This is a :py:class:`valuesModule.Values` instance that has the data of the array.            |
    +-------------------+-----------------------------------------------------------------------------------------------+
    """

    compression = 'full'
    ancestryMembers = baseModule.XDataCoreMembers.ancestryMembers + ( 'values', )

    def __init__(self, shape=None, data=None, symmetry=Symmetry.none, storageOrder=enumsModule.StorageOrder.rowMajor, 
                offset=None, permutation=Permutation.plus, index=None, label=None):
        """
        :param shape:               This is the shape of the array.
        :param data:                This is a :py:class:`valuesModule.Values` instance that has the data of the array.
        :param symmetry:            This is the symmtry (none, lower or upper) for the data in the array.
        :param storageOrder:        This is the storage order (row or column) for the data in the arrary.
        :param offset:              This is the offset index where the first data value in the array starts.
        :param permutation:         This is the permutation for the data in the array.
        :param index:               This is index of the array inside of other xData function.
        :param label:               This is the label for the array.
        """

        ArrayBase.__init__( self, shape, symmetry = symmetry, storageOrder = storageOrder, 
                offset = offset, permutation = permutation,
                index = index, label = label )

        if( not( isinstance( data, valuesModule.Values ) ) ) : data = valuesModule.Values( data )

        if symmetry == Symmetry.none:
            size = self.size
        else :                  # Will be a 'square' array. Checked in ArrayBase.
            length, size = self.shape[0], 1
            for i1 in range( self.dimension ) : size = size * ( length + i1 ) / ( i1 + 1 )
        if( size != len( data ) ) : raise ValueError( 'shape requires %d values while data has %d values' % ( size, len( data ) ) )

        self.values = data
        self.values.setAncestor( self )

    def constructArray( self ) :
        """
        This method returns a numpy array that represents the data in *self* with all cells fulled.

        :returns:       An instance of :py:class:`numpy.array`.
        """

        import numpy

        if self.symmetry == Symmetry.none:
            array1 = numpy.array( [ value for value in self.values ] )

        elif len(self.shape) == 2:
            array1 = numpy.zeros( self.shape )
            if self.symmetry == Symmetry.lower:
                array1[numpy.tril_indices(self.shape[0])] = list(self.values)
                il = numpy.tril_indices(self.shape[0], -1)
                iu = (il[1],il[0])
                array1[iu] = array1[il]
            elif self.symmetry == Symmetry.upper:
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
            mode = ( ( self.symmetry == Symmetry.upper ) and ( self.storageOrder == enumsModule.StorageOrder.rowMajor ) ) or \
                     ( self.symmetry == Symmetry.lower ) and ( self.storageOrder == enumsModule.StorageOrder.columnMajor )
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

        order = { enumsModule.StorageOrder.rowMajor : 'C', enumsModule.StorageOrder.columnMajor : 'F' }[self.storageOrder]
        return( array1.reshape( self.shape, order = order ) )

    def copy( self ) :
        """
        This method returns a copy of *self*.

        :returns:       An instance of :py:class:`Full`.
        """

        return( Full( self.shape, self.values.copy( ), 
                symmetry = self.symmetry, storageOrder = self.storageOrder, 
                offset = self.offset, permutation = self.permutation,
                index = self.index, label = self.label ) )

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        attributesStr = self.attributesToXMLAttributeStr()
        XML_strList = [ '%s<%s%s>' % ( indent, self.moniker, attributesStr ) ]
        XML_strList += self.values.toXML_strList(indent2, **kwargs)
        XML_strList[-1] += '</%s>' % self.moniker

        return XML_strList

class Diagonal( ArrayBase ) :
    """
    This class store the array with GNDS diagonal compression where ohly some of the diagonals have non-zero data. 
    Multiple diagonals can be stored. If more than the main diagonal is represented, the starting point of each 
    diagaonal by me specified by the *startingIndices* member as per the GNDS specifications.

    The following table list the primary members of this class that are not specified in the primary members
    table in :py:class:`ArrayBase`:

    +-------------------+-----------------------------------------------------------------------------------------------+
    | Member            | Description                                                                                   |
    +===================+===============================================================================================+
    | data              | This is a :py:class:`valuesModule.Values` instance that has the data of the array.            |
    +-------------------+-----------------------------------------------------------------------------------------------+
    | startingIndices   | This is a :py:class:`valuesModule.Values` instance specifies where each diagonal starts.      |
    +-------------------+-----------------------------------------------------------------------------------------------+
    """

    compression = 'diagonal'
    ancestryMembers = baseModule.XDataCoreMembers.ancestryMembers + ( 'values', )

    def __init__( self, shape = None, data = None, startingIndices = None, symmetry = Symmetry.none, storageOrder = enumsModule.StorageOrder.rowMajor, 
                offset = None, permutation = Permutation.plus,
                index = None, label = None ) :
        """
        :param shape:               This is the shape of the array.
        :param data:                This is a :py:class:`valuesModule.Values` instance that has the data of the array.
        :param startingIndices:     This is a :py:class:`valuesModule.Values`
        :param symmetry:            This is the symmtry (none, lower or upper) for the data in the array.
        :param storageOrder:        This is the storage order (row or column) for the data in the arrary.
        :param offset:              This is the offset index where the first data value in the array starts.
        :param permutation:         This is the permutation for the data in the array.
        :param index:               This is index of the array inside of other xData function.
        :param label:               This is the label for the array.
        """

        ArrayBase.__init__( self, shape, symmetry, storageOrder = storageOrder,
                offset = offset, permutation = permutation,
                index = index, label = label )

        dimension = self.dimension
        if( not( isinstance( data, valuesModule.Values ) ) ) : data = valuesModule.Values( data )
        if( startingIndices is None ) :
            self.startingIndicesOriginal = None
            startingIndices = dimension * [ 0 ]
        else :
            if not isinstance(startingIndices, valuesModule.Values):
                startingIndices = valuesModule.Values(startingIndices, valueType=enumsModule.ValueType.integer32)
            if startingIndices.valueType not in [enumsModule.ValueType.integer32]:
                raise TypeError('startingIndices must be a list of integers.')
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
        """
        This method returns a numpy array that represents the data in *self* with all cells fulled.

        :returns:       An instance of :py:class:`numpy.array`.
        """

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
                if self.symmetry == Symmetry.none:
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

        order = { enumsModule.StorageOrder.rowMajor: 'C', enumsModule.StorageOrder.columnMajor: 'F' }[self.storageOrder]
        return( array1.reshape( self.shape, order = order ) )

    def copy( self ) :
        """
        This method returns a copy of *self*.

        :returns:       An instance of :py:class:`Diagonal`.
        """

        startingIndicesOriginal = self.startingIndicesOriginal
        if( startingIndicesOriginal is not None ) : startingIndicesOriginal = self.startingIndicesOriginal.copy( )
        return( Diagonal( self.shape, self.values.copy( ), startingIndicesOriginal, 
                symmetry = self.symmetry, storageOrder = self.storageOrder,
                offset = self.offset, permutation = self.permutation,
                index = self.index, label = self.label ) )

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributesStr = self.attributesToXMLAttributeStr()
        XML_strList = [ '%s<%s%s>' % ( indent, self.moniker, attributesStr ) ]
        if self.startingIndicesOriginal is not None:
            XML_strList += self.startingIndicesOriginal.toXML_strList(indent2, **kwargs)
        XML_strList += self.values.toXML_strList(indent2, **kwargs)
        XML_strList[-1] += '</%s>' % self.moniker

        return XML_strList

class Flattened( ArrayBase ) :
    """
    This class store the array with GNDS flattened compression. This array is good for storing an array where the
    non-zero parts of the array (called chunk herein) are sparse with many zero cells between each chunk.
    In addition to each chunk, the starting position and length of each chunk is needed. A starting position is
    is the flat index of the start of the chunk.

    The following table list the primary members of this class that are not specified in the primary members
    table in :py:class:`ArrayBase`:

    +-------------------+-----------------------------------------------------------------------------------------------+
    | Member            | Description                                                                                   |
    +===================+===============================================================================================+
    | data              | This is a :py:class:`valuesModule.Values` instance that has the data of the array.            |
    +-------------------+-----------------------------------------------------------------------------------------------+
    | starts            | This is a :py:class:`valuesModule.Values` instance that specifies where each chuck starts.    |
    +-------------------+-----------------------------------------------------------------------------------------------+
    | lengths           | This is a :py:class:`valuesModule.Values` instance that specifies the length of each chunk.   |
    +-------------------+-----------------------------------------------------------------------------------------------+
    """

    compression = 'flattened'
    ancestryMembers = baseModule.XDataCoreMembers.ancestryMembers + ( 'values', 'lengths', 'starts' )

    def __init__( self, shape = None, data = None, starts = None, lengths = None, symmetry = Symmetry.none, 
                storageOrder = enumsModule.StorageOrder.rowMajor,
                offset = None, permutation = Permutation.plus,
                index = None, label = None, 
                dataToString = None ) :
        """
        :param shape:               This is the shape of the array.
        :param data:                This is a :py:class:`valuesModule.Values` instance that has the data of the array.
        :param starts:              This is a :py:class:`valuesModule.Values` instance that contains the start flat index of each data chunk.
        :param lengths:             This is a :py:class:`valuesModule.Values` instance that contains the length of each data chunk.
        :param symmetry:            This is the symmtry (none, lower or upper) for the data in the array.
        :param storageOrder:        This is the storage order (row or column) for the data in the arrary.
        :param offset:              This is the offset index where the first data value in the array starts.
        :param permutation:         This is the permutation for the data in the array.
        :param index:               This is index of the array inside of other xData function.
        :param label:               This is the label for the array.
        """

        ArrayBase.__init__( self, shape, symmetry, storageOrder = storageOrder,
                offset = offset, permutation = permutation,
                index = index, label = label )

        self.dataToString = dataToString

        if( not( isinstance( data, valuesModule.Values ) ) ) : data = valuesModule.Values( data )
        if( not( isinstance( starts, valuesModule.Values ) ) ) : 
                starts = valuesModule.Values(starts, valueType=enumsModule.ValueType.integer32)
        if( not( isinstance( lengths, valuesModule.Values ) ) ) :
                lengths = valuesModule.Values(lengths, valueType=enumsModule.ValueType.integer32)

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
        """
        This method returns a numpy array that represents the data in *self* with all cells fulled.

        :returns:       An instance of :py:class:`numpy.array`.
        """

        import numpy

        index = 0
        array1 = numpy.zeros( self.size )
        for i1, start in enumerate( self.starts ) :
            length = self.lengths[i1]
            for i2 in range( length ) :
                array1[start+i2] = self.values[index]
                index += 1

        order = { enumsModule.StorageOrder.rowMajor: 'C', enumsModule.StorageOrder.columnMajor: 'F' }[self.storageOrder]
        array1 = array1.reshape( self.shape, order = order )
        if self.symmetry == Symmetry.lower:
            array1 = numpy.tril(array1) + numpy.tril(array1, -1).T
        elif self.symmetry == Symmetry.upper:
            array1 = numpy.triu(array1) + numpy.triu(array1, -1).T
        return array1

    def copy( self ) :
        """
        This method returns a copy of *self*.

        :returns:       An instance of :py:class:`Flattened`.
        """

        return( Flattened( self.shape, self.values.copy( ), self.starts.copy( ), self.lengths.copy( ), 
                symmetry = self.symmetry, storageOrder = self.storageOrder,
                offset = self.offset, permutation = self.permutation,
                index = self.index, label = self.label ) )

    @staticmethod
    def fromNumpyArray( array, symmetry = Symmetry.none, nzeroes = 4 ):
        """
        This static method generates a sparse flattened array that represents an arbitrary numpy array.
        Tjios method currently only supports 'full' or 'lower-symmetric' matrices, with row-major data storage.

        :param array:           Input numpy array to flatten.
        :param symmetry:        Allowed values are 'none' or 'lower'
        :param nzeroes:         How many zeroes to allow before adding a new 'start' and 'length' 

        :return:                An instance of :py:class:`Flattened`.
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

        if symmetry == Symmetry.none:
            helper( array.flatten() )
        elif symmetry == Symmetry.lower:
            rows, cols = array.shape
            for row in range(rows):
                dat = array[row][:row+1]
                helper( dat, offset = row * cols )
        else:
            raise NotImplementedError("Symmetry = '%s'" % symmetry)

        return Flattened(shape=array.shape, data=sparseData, starts=starts, lengths=lengths, symmetry=symmetry)


    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        attributesStr = self.attributesToXMLAttributeStr()
        XML_strList = [ '%s<%s%s>' % ( indent, self.moniker, attributesStr ) ]
        XML_strList += self.starts.toXML_strList(indent2, **kwargs)
        XML_strList += self.lengths.toXML_strList(indent2, **kwargs)
        XML_strList += self.values.toXML_strList(indent2, **kwargs)
        XML_strList[-1] += '</%s>' % self.moniker
        return XML_strList

class Embedded( ArrayBase ) :
    """
    This class store the array with GNDS embeeded compression where the array is stores as a list of sub-arrays.

    The primary members of this class are list in the primary members table in :py:class:`ArrayBase`:
    """

    compression = 'embedded'

    def __init__(self, shape=None, symmetry=Symmetry.none, storageOrder=enumsModule.StorageOrder.rowMajor,
                offset=None, permutation=Permutation.plus, index=None, label=None):
        """
        :param shape:               This is the shape of the array.
        :param storageOrder:        This is the storage order (row or column) for the data in the arrary.
        :param symmetry:            This is the symmtry (none, lower or upper) for the data in the array.
        :param offset:              This is the offset index where the first data value in the array starts.
        :param permutation:         This is the permutation for the data in the array.
        :param index:               This is index of the array inside of other xData function.
        :param label:               This is the label for the array.
        """

        ArrayBase.__init__( self, shape, symmetry, storageOrder = storageOrder, 
                offset = offset, permutation = permutation,
                index = index, label = label )

        self.arrays = []

    def addArray( self, array ) :
        """
        Thie method adds an :py:class:`ArrayBase` instance to *self*.

        :param array:       An :py:class:`ArrayBase` instance.
        """

        if( not( isinstance( array, ArrayBase ) ) ) : raise TypeError( 'variable not an array instance' )
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
        """
        This method returns a numpy array that represents the data in *self* with all cells fulled.

        :returns:       An instance of :py:class:`numpy.array`.
        """

        import numpy

        order = { enumsModule.StorageOrder.rowMajor: 'C', enumsModule.StorageOrder.columnMajor: 'F' }[self.storageOrder]
        array1 = numpy.zeros( self.shape, order = order )

        for array in self.arrays :
            array2 = array.constructArray( )
            slice1 = [ slice( offset, offset + array2.shape[i1] ) for i1, offset in enumerate( array.offset ) ]
            array1[slice1] = array2

        return( array1 )

    def copy( self ) :
        """
        This method returns a copy of *self*.

        :returns:       An instance of :py:class:`Embedded`.
        """

        array1 = Embedded( self.shape, 
                symmetry = self.symmetry, storageOrder = self.storageOrder,
                offset = self.offset, permutation = self.permutation,
                index = self.index, label = self.label )
        for array in self.arrays : array1.addArray( array.copy( ) )
        return array1

    def offsetScaleValues( self, offset, scale ):
        """
        This method class *offsetScaleValues* for each sub-array of *self*.

        :param offset:      The offset to apply to each cell of each sub-array.
        :param scale:       The multiplicative scale factor to apply to each cell of each sub-array.
        """

        for subarray in self.arrays: subarray.offsetScaleValues( offset, scale )

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        attributesStr = self.attributesToXMLAttributeStr()
        XML_strList = [ '%s<%s%s>' % ( indent, self.moniker, attributesStr ) ]
        for array in self.arrays: toXML_strList += array.toXML_strList(indent2, **kwargs)
        toXML_strList[-1] += '</%s>' % self.moniker
        return toXML_strList
