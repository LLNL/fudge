# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module containes the class for storing a mathematical matrix (i.e., 2d-array), and support several matrix operations on the matrix.

This module contains the following classes:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | Matrix                            | This class stores a mathematical matrix.                              |
    +-----------------------------------+-----------------------------------------------------------------------+
"""

import numpy

from . import vector as vectorModule

class Matrix():
    """
    This class stores a matrix (i.e., 2d-array) and can perform several matrix operations on the matrix.

    The following table list the primary members of this class:

    +-----------+---------------------------------------------------------------+
    | Member    | Description                                                   |
    +===========+===============================================================+
    | matrix    | The data of the matrix.                                       |
    +-----------+---------------------------------------------------------------+
    """

    def __init__(self, numberRows=0, numberColumns=0, values=None):
        """
        The input arguments are optional and the default behaviour is and instance of Matrix with an empty matrix.
        Either both or neither the numberRows and numberColumns arguments need to be speciified. If only the numberRows 
        and numberColumns arguments are provided, self.matrix is a 2-D numpy array with shape (numberRows, numberColumns). 
        If only the values argument is given, self.matriix is a 2-D numpy.array with these values as entries. If
        both the size and values arguments are given, the latter is expected to be a scalar.

        :param numberRows:      The number of rows of the matriix.
        :param numberColumns:   The number of columns of the matrix.
        :param values:          The values of the matrix.
        """

        if values is None and 0 in (numberRows, numberColumns):
            assert numberRows == numberColumns, f'When specified, the numberRows and numberColumns arguments should either be both zero or non-zero.'

            self.matrix = numpy.empty((0,0))

        elif numberRows > 0 and numberColumns > 0:
            self.matrix = numpy.zeros((numberRows, numberColumns))

            if isinstance(values, (int, float)):
                self.matrix += numpy.float64(values)

            elif values is not None:
                raise TypeError('Only scalar value of the named "values" argument is allowed to accompany non-zero instances of the numberRows and numberColumns arguments.')

        elif values is not None:
            if numberRows != 0 or numberColumns != 0:
                raise ValueError(f'When values is not None, numberRows (%s) and numberColumns (%s) must be zero.' % (numberRows, numberColumns))

            if isinstance(values, (list, tuple)):
                values = numpy.array(values, ndmin=2)

            if isinstance(values, numpy.ndarray):
                if values.ndim == 2:
                    self.matrix = values
                else:
                    raise TypeError(f'Invalid "values" named argument to {type(self)}. Numpy ndarrays should have a dimension of 2.')

            elif isinstance(values, Matrix):
                self.matrix = values.matrix.copy()

            else:
                raise TypeError(f'Invalid "values" named argument to {type(self)}. Either provide a list of lists, a 2-dimensional Numpy ndarray or another instance of {type(self)}.')

    def __len__(self):
        """
        This method returns the number of rows of *self*.

        :returns:       A python int.
        """

        return(len(self.matrix))

    def __getitem__(self, indexTuple):
        """
        This method returns the value of *self* at the location specified by the input *indexTuple*. If *indexTuple* is a number
        (e.g., called as self[1]) then a row is returned. If *indexTuple* is a tuple with two integers (e.g., called as self[1,2]), 
        then the value at that cell is returned. In general, any numpy indexing for a 2-d array is allowed.

        :param indexTuple:  Row or cell index of the data to return.

        :returns:           An instance of :py:class:`numpy.ndarray` or :py:class:`numpy.float64`.
        """

        return self.matrix[indexTuple]

    def __setitem__(self, indexTuple, newScalarValue):
        """
        This method replace the value at specified cell. 

        :param indexTuple:      Row or cell index of the data to replace.
        :param newScalarValue:  Replacement for value of self.matrix at given location.

        :returns:               A reference to self.
        """

        vectorModule.checkScalarValidity(newScalarValue, 'item assignment')

        self.matrix[indexTuple] = numpy.float64(newScalarValue)

        return self

    def __add__(self, otherMatrixOrScalar):
        """
        Thie method adds the matrix or scalar *otherMatrixOrScalar* to self.

        :param otherMatrixOrScalar: :py:class:`Matrix` or scalar to add.

        :returns:                   A new instance of :py:class:`Matrix`.
        """

        if len(self) == 0:
            if isinstance(otherMatrixOrScalar, Matrix):
                return Matrix(values=otherMatrixOrScalar.matrix)
            else:
                raise TypeError(f'Seconds argument is required to be of type {type(self)}.')

        else:
            if isinstance(otherMatrixOrScalar, Matrix):
                if len(otherMatrixOrScalar) == 0:
                    return Matrix(values=self.matrix)
                else:
                    return Matrix(values=self.matrix+otherMatrixOrScalar.matrix)

            elif vectorModule.checkScalarValidity(otherMatrixOrScalar):
                return (Matrix(values=self.matrix+otherMatrixOrScalar))

            else:
                raise TypeError(f'Adding a {type(otherMatrixOrScalar)} to a {type(self)} is not allowed.')

    def __iadd__(self, otherMatrixOrScalar):
        """
        Thie method adds a matrix or scalar to self (i.e., in-place addition).
`
        :param otherMatrixOrScalar:     :py:class:`Matrix` or scalar.

        :returns:                       A reference to self.
        """
        self = self.__add__(otherMatrixOrScalar)

        return self

    __radd__ = __add__

    def __sub__(self, otherMatrixOrScalar):
        """
        This method subtract the matrix or scalar given in the input argument from self.

        :param otherMatrixOrScalar:     :py:class:`Matrix` or scalar to subract.

        :returns:                       A new instance of :py:class:`Matrix`.
        """
        return self.__add__(-otherMatrixOrScalar)

    def __rsub__(self, otherMatrixOrScalar):
        """
        This method subtract self from the matrix or scalar given in the input argument.

        :param otherMatrixOrScalar:     :py:class:`Matrix` or scalar.

        :returns:                       A new instance of :py:class:`Matrix`.
        """
        return -self.__sub__(otherMatrixOrScalar)

    def __neg__(self):
        """
        This method returns the negations of all the value of self as a new :py:class:`Matrix` instance.

        :returns:           A new instance of :py:class:`Matrix`.
        """
        return Matrix(values=-self.matrix)

    def __mul__(self, scalarValue):
        """
        This method returns a new :py:class:`Matrix` instance that is *self* multiplied by a scalar.

        :param scalarValue:     Multiplication scalar.

        :returns:               A new instance of :py:class:`Matrix`.
        """
        vectorModule.checkScalarValidity(scalarValue, 'multiplication')

        return Matrix(values=scalarValue*self.matrix)

    __rmul__ = __mul__

    def __imul__(self, scalarValue):
        """
        This method multiplies all value in self by a scaler (i.e., in-place multiplication).

        :param scalarValue:     Multiplication scalar.

        :returns:               A reference to self.
        """
        self = self.__mul__(scalarValue)

        return self

    def __truediv__(self, scalarValue):
        """
        This method returns a new :py:class:`Matrix` instance with its values being the values of *self* divided the scalar *scalarValue*.

        :param scalarValue:     Division scalar.

        :returns:               A new instance of :py:class:`Matrix`.
        """
        vectorModule.checkScalarValidity(scalarValue, 'division')
        with numpy.errstate(divide='raise'):
            return Matrix(values=self.matrix/scalarValue)

    def __itruediv__(self, scalarValue):
        """
        Thie method divides all values of *self* by the scalar *scalarValue* (i.e., in-place division).

        :param scalarValue:     Division scalar.

        :returns:               A reference to self.
        """
        self = self.__truediv__(scalarValue)

        return self

    def __str__(self):
        """
        Thie method returns a string representation of self.

        :returns:           A python str)
        """
        return str(self.matrix)

    @property
    def shape(self):
        """
        This method returns the shape of *self*.

        :returns:           A tuple of two integers.
        """

        return self.matrix.shape

    @property
    def numberOfRows(self):
        """
        This method returns the number of rows of *self*. This is the same as len(self).

        :returns:           A python int.
        """
        return self.matrix.shape[0]

    @property
    def numberOfColumns(self):
        """
        This method returns the number of columns of *self*.

        :returns:           A python int.
        """
        return self.matrix.shape[1]

    def append(self, newRow):
        """
        This method appends a row to self, increasing its number of rows by 1.

        :param newRow:      New vector to append as a new row to self.
        """

        if isinstance(newRow, (list, tuple, numpy.ndarray)):
            newRow = vectorModule.Vector(newRow)

        if isinstance(newRow, vectorModule.Vector):
            if len(self) != 0 and len(newRow) != self.numberOfColumns:
                raise TypeError(f'Input {vectorModule.Vector} argument must be the same length as the number of columns.')
            rowToAppend = newRow.vector
        else:
            raise TypeError(f'Input argument of type {vectorModule.Vector}, {numpy.ndarray}, list or tuple expected.')

        if len(self) != 0:
            self.matrix = numpy.vstack([self.matrix, rowToAppend])
        else:
            self.matrix = numpy.array(rowToAppend, ndmin=2)

    def transpose(self):
        """
        This method returns the transpose of self.

        :returns:           A new instance of :py:class:`Matrix`.
        """
        return Matrix(values=numpy.transpose(self.matrix))

    def reverse(self):
        """
        This method returns a new :py:class:`Matrix` with its cells reversed along both the row and column dimensions.

        :returns:           A new instance of :py:class:`Matrix`.
        """
        return Matrix(values=numpy.flipud(numpy.fliplr(self.matrix)))
