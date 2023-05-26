# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import numpy

from . import vector as vectorModule

class Matrix():
    """Class to store mathematical matrix that perform several matrix operations."""

    def __init__(self, numberRows=0, numberColumns=0, values=None):
        """
        Initialize an instance of Matrix.

        The input arguments are optional and the default behaviour is and instance of Matrix with self.matrix = None.
        Either both or neither the numberRows and numberColumns arguments need to be speciified. If only the numberRows 
        and numberColumns arguments are provided, self.matrix is a 2-D numpy array with shape (numberRows, numberColumns). 
        If only the values argument is given, self.matriix is a 2-D numpy.array with these values as entries. If
        both the size and values arguments are given, the latter is expected to be a scalar.

        The optional input arguments are as follows:
        :numberRows: The number of rows of the matriix (default = 0).
        :numberColumns: The number of columns of the matrix (default = 0).
        :values: The values of the matrix (default = None).
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
        """Returns the number of rows of *self*."""

        return(len(self.matrix))

    def __getitem__(self, indexTuple):
        """
        Return value of self.matrix at location specified by the input tuple.

        :indexTuple: Row and column indices of the self.matrix value to return.
        :returns: A numpy.float64 value.
        """

        return self.matrix[indexTuple]

    def __setitem__(self, indexTuple, newScalarValue):
        """
        Replace value at specified location of self.matrix.

        :indexTuple: Row and column indices of the self.matrix value to replace.
        :newScalarValue: Replacement for value of self.matrix at given location.
        """

        vectorModule.checkScalarValidity(newScalarValue, 'item assignment')

        self.matrix[indexTuple] = numpy.float64(newScalarValue)

        return self

    def __add__(self, otherMatrixOrScalar):
        """
        Add the matrix or scalar given in the input argument to self.matrix.

        :otherMatrixOrScalar: Vector or scalar to add.
        :returns: A new instance of Matrix.
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
        In-place addition.
        """
        self = self.__add__(otherMatrixOrScalar)

        return self

    __radd__ = __add__

    def __sub__(self, otherMatrixOrScalar):
        """
        Subtract the matrix or scalar given in the input argument from self.matrix.

        :otherMatrixOrScalar: Matrix or scalar to subract.
        :returns: A new instance of Matrix.
        """
        return self.__add__(-otherMatrixOrScalar)

    def __rsub__(self, otherMatrixOrScalar):
        """
        Subtract self.matrix from the matrix given in the input argument.
        """
        return -self.__sub__(otherMatrixOrScalar)

    def __neg__(self):
        """
        Negation of self.matrix.

        :returns: A new instance of Matrix.
        """
        return Matrix(values=-self.matrix)

    def __mul__(self, scalarValue):
        """
        Multiply a scalar to self.matrix.

        :scalarValue: Multiplication scalar.
        :returns: A new instance of Matrix.
        """
        vectorModule.checkScalarValidity(scalarValue, 'multiplication')

        return Matrix(values=scalarValue*self.matrix)

    __rmul__ = __mul__

    def __imul__(self, scalarValue):
        """
        In-place multiplication

        :scalarValue: Multiplication scalar.
        """
        self = self.__mul__(scalarValue)

        return self

    def __truediv__(self, scalarValue):
        """
        True division of self.matrix by the input scalar argument.

        :scalarValue: Division scalar.
        :returns: A new instance of Matrix.
        """
        vectorModule.checkScalarValidity(scalarValue, 'division')
        with numpy.errstate(divide='raise'):
            return Matrix(values=self.matrix/scalarValue)

    def __itruediv__(self, scalarValue):
        """
        In-place true division.
        """
        self = self.__truediv__(scalarValue)

        return self

    def __str__(self):
        """
        Return a string representation of self.matrix.
        """
        return str(self.matrix)

    @property
    def shape(self):
        """
        Return the shape of self.matrix.
        """

        return self.matrix.shape

    @property
    def numberOfRows(self):
        """
        Return the number of rows of self.matrix. This is the same as len(self).
        """
        return self.matrix.shape[0]

    @property
    def numberOfColumns(self):
        """
        Return the number of columns of self.matrix.
        """
        return self.matrix.shape[1]

    def append(self, newRow):
        """
        Append row to self.matrix.

        :newRow: New vector to append as a new row to self.matrix.
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
        Transpose self.matrix.

        :returns: A new instance of Matrix
        """
        return Matrix(values=numpy.transpose(self.matrix))

    def reverse(self):
        """
        Reverse the order of elements long both the row and column dimensions.

        :returns: A new instance of Matrix
        """
        return Matrix(values=numpy.flipud(numpy.fliplr(self.matrix)))
