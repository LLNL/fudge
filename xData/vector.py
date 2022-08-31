# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import numpy

class Vector:
    """Class to store a mathematical vector that perform several vector operations."""

    def __init__(self, size=0, values=None):
        """
        Initialize an instance of Vector.

        The input arguments are optional and the default behaviour is an instance of Vector with self.vector = numpy.array([]).
        If only the size argument is provided, self.vector is a 1-D numpy array with length corresponding to this size
        argument. If only the values argument is given, self.vector is a 1-D numpy.array with these values as entries. If
        both the size and values arguments are given, the latter is expected to be a number.

        The optional input arguments are as follows:
        :size: The size of the  vector (default = 0).
        :values: The values of the vector (default = None).
        """

        if size == 0 and values is None:
            self.vector = numpy.array([])
   
        elif size > 0:
            if values is None:
                self.vector = numpy.zeros(size)
            else:
                if isinstance(values, (int, float)):
                    self.vector = numpy.full(size, numpy.float(values))
                else:
                    raise TypeError('Only an int or float value of the named "values" argument is allowed to accompany non-zero instances of the size argument.')

        else:
            if isinstance(values, (list, tuple)):
                values = numpy.array(values).astype(numpy.float)

            if isinstance(values, numpy.ndarray) and values.ndim == 1:
                self.vector = values.astype(numpy.float)
            elif isinstance(values, Vector):
                self.vector = values.vector.copy()
            else:
                raise TypeError(f'Optional second argument is of type {type(values)}. It is expected to be a list, tuple or 1-D numpy.ndarray.')
    
    def __add__(self, otherVectorOrScalar):
        """
        Add scalar or Vector given in the input argument to self.

        :otherVectorOrScalar: Vector or scalar to add.
        :returns: A new instance of Vector.
        """

        if len(self) == 0:
            if isinstance(otherVectorOrScalar, Vector):
                return Vector(values=otherVectorOrScalar.vector)
            raise TypeError(f'Seconds argument is required to be of type {type(self)}.')
        else:
            if isinstance(otherVectorOrScalar, Vector):
                if len(otherVectorOrScalar) == 0:
                    return Vector(values=self.vector)
                return Vector(values=self.vector+otherVectorOrScalar.vector)
            elif isinstance(otherVectorOrScalar, (list, tuple)):
                return Vector(values=self.vector+otherVectorOrScalar)
            else:
                raise TypeError(f'Adding a {type(otherVectorOrScalar)} to a {type(self)} is not allowed.')

    def __iadd__(self, otherVectorOrScalar):
        """
        In-place addition.
        """

        self = self.__add__(otherVectorOrScalar)
            
        return self

    __radd__ = __add__

    def __sub__(self, otherVectorOrScalar):
        """
        Subtract the vector or scalar given in the input argument from self.vector.

        :otherVectorOrScalar: Vector or scalar to subract.
        :returns: A new instance of Vector.
        """
        return self.__add__(-otherVectorOrScalar)

    def __rsub__(self, otherVectorOrScalar):
        """
        Subtract input argument from self.vector.
        """

        return -self.__sub__(otherVectorOrScalar)

    def __mul__(self, scalarValue):
        """
        Multiply self.vector by a scalar.

        :scalarValue: Multiplication scalar.
        :returns: A new instance of Vector.
        """
        checkScalarValidity(scalarValue, 'multiplication')

        return Vector(values=scalarValue*self.vector)

    __rmul__ = __mul__

    def __imul__(self, scalarValue):
        """
        In-place multiplication

        :scalarValue: Multiplication scalar.
        """
        self = self.__mul__(scalarValue)

        return self

    def __neg__(self):
        """
        Negation of self.vector.

        :returns: A new instance of Vector.
        """

        return Vector(values=-self.vector)

    def __truediv__(self, scalarValue):
        """
        True division of self.vector by the input scalar argument.

        :scalarValue: Division scalar.
        :returns: A new instance of Vector.
        """
        checkScalarValidity(scalarValue, 'division')
        with numpy.errstate(divide='raise'):
            return Vector(values=self.vector/scalarValue)

    def __itruediv__(self, scalarValue):
        """
        In-place true division.
        """
        self = self.__truediv__(scalarValue)

        return self

    def __getitem__(self, index):
        """
        Return value of self.vector at location specified by the input scalar argument.

        :index: Index of the self.vector value to return.
        :returns: A numpy.float value.
        """
        return self.vector[index]

    def __setitem__(self, index, newScalarValue):
        """
        Replace value at specified location of self.vector.

        :index: Index of the self.vector value to replace.
        :newScalarValue: Replacement for value of self.vector at given location.
        """
        checkScalarValidity(newScalarValue, 'item assignment')

        self.vector[index] = numpy.float(newScalarValue)
        return self

    def __len__(self):
        """
        Returns the length of self.vector.
        """

        return len(self.vector)

    def __str__(self):
        """
        Return a string representation of self.vector.
        """
        return str(self.vector)

    @property
    def sum(self):
        """
        Return the sum of self.vector.
        """
        return self.vector.sum()

    @property
    def isSorted(self):
        """
        Return a boolean indicating whether self.vector is sorted.
        """
        return all(numpy.sort(self.vector) == self.vector)

    def reverse(self):
        """
        Return self.vector in reversed order.
        """

        def copy(self):
            """
            Return a copy of self.
            """
            return Vector(values=self.vector)
        self.vector = self.vector[::-1]

def checkScalarValidity(scalarValue, operation):
    """
    Check validity of a scalar argument.

    :scalarValue: The scalar value to test.
    :operation: A string with the name of the operation/method in which the string argument is used.
    """

    if not isinstance(scalarValue, (int, float)):
        raise TypeError(f'Only {operation} with scalar value is allowed.')

    elif numpy.isinf(scalarValue):
        raise ValueError(f'Attempting {operation} by {numpy.inf} scalar value.')

    elif numpy.isnan(scalarValue):
        raise ValueError(f'Attempting {operation} by {numpy.nan} scalar value.')
