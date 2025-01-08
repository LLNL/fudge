# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module containes all the classes for handling GNDS axes and its child nodes.

This module contains the following classes:

    +---------------+-------------------------------------------------------------------------------------------+
    | Class         | Description                                                                               |
    +===============+===================+=======================================================================+
    | Vector        | This class is used to represent multi-group data and to perform some math operations on   |
    |               | the data.                                                                                 |
    +---------------+-------------------------------------------------------------------------------------------+
"""

import numpy

class Vector:
    """
    This class is used to represent multi-group data and to perform some math operations on the data. In essence,
    multi-group data are a set of numbers. All of the multi-group data can be added, subtracted, multiplied or divided 
    by a number (i.e., a scalar). In addition, two instances of :py:class:`Vector` can be added or subtracted.
    In this latter case, both instances must be the same length or at least one must have a length of 0.
    When both instances have the same length, each datum in one instance is added (substracted) to the 
    datum in the other instance at the same index. When one instances has length of 0, it is treated as
    having the same length as the other instance but with all its data being 0.0's.

    The following table list the primary members of this class:

    +-----------+---------------------------------------------------------------+
    | Member    | Description                                                   |
    +===========+===============================================================+
    | vector    | The list of floats.                                           |
    +-----------+---------------------------------------------------------------+
    """

    def __init__(self, size=0, values=None):
        """
        The input arguments are optional and the default behaviour is an instance of Vector with self.vector = numpy.array([]).
        If only the size argument is provided, self.vector is a 1-D numpy array with length corresponding to this size
        argument. If only the values argument is given, self.vector is a 1-D numpy.array with these values as entries. If
        both the size and values arguments are given, the latter is expected to be a number.

        :size:      The size of the  vector (default = 0).
        :values:    The values of the vector (default = None).
        """

        if size == 0 and values is None:
            self.vector = numpy.array([])
   
        elif size > 0:
            if values is None:
                self.vector = numpy.zeros(size)
            else:
                if isinstance(values, (int, float)):
                    self.vector = numpy.full(size, numpy.float64(values))
                else:
                    raise TypeError('Only an int or float value of the named "values" argument is allowed to accompany non-zero instances of the size argument.')

        else:
            if isinstance(values, (list, tuple)):
                values = numpy.array(values).astype(numpy.float64)

            if isinstance(values, numpy.ndarray) and values.ndim == 1:
                self.vector = values.astype(numpy.float64)
            elif isinstance(values, Vector):
                self.vector = values.vector.copy()
            else:
                raise TypeError(f'Optional second argument is of type {type(values)}. It is expected to be a list, tuple or 1-D numpy.ndarray.')
    
    def __add__(self, otherVectorOrScalar):
        """
        This method adds *otherVectorOrScalar* to *self*.

        :otherVectorOrScalar:   :py:class:`Vector` instance or a scalar.

        :returns:               A new instance of :py:class:`Vector`.
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
        In-place addition. See :py:func:`__add__`.

        :otherVectorOrScalar:   :py:class:`Vector` or scalar.

        :returns:               This method returns *self*.
        """

        self = self.__add__(otherVectorOrScalar)
            
        return self

    __radd__ = __add__

    def __sub__(self, otherVectorOrScalar):
        """
        This method subtracts *otherVectorOrScalar* from *self*.

        :otherVectorOrScalar:   :py:class:`Vector` or scalar.

        :returns:               A new instance of :py:class:`Vector`.
        """

        return self.__add__(-otherVectorOrScalar)

    def __rsub__(self, otherVectorOrScalar):
        """
        This method subtracts *self* from *otherVectorOrScalar*.

        :otherVectorOrScalar:   :py:class:`Vector` or scalar.

        :returns:               A new instance of :py:class:`Vector`.
        """

        return -self.__sub__(otherVectorOrScalar)

    def __mul__(self, scalarValue):
        """
        This method multiply *self* by scalar.

        :scalarValue:           Multiplication scalar.

        :returns:               A new instance of :py:class:`Vector`.
        """
        checkScalarValidity(scalarValue, 'multiplication')

        return Vector(values=scalarValue*self.vector)

    __rmul__ = __mul__

    def __imul__(self, scalarValue):
        """
        In-place multiplication. See :py:func:`__mul__'.

        :scalarValue:           Multiplication scalar.
        
        :returns:               This method returns *self*.
        """
        self = self.__mul__(scalarValue)

        return self

    def __neg__(self):
        """
        Thie method returns new instance of :py:class:`Vector` whose values of the negation of *self*'s values.

        :returns:               A new instance of :py:class:`Vector`.
        """

        return Vector(values=-self.vector)

    def __truediv__(self, scalarValue):
        """
        True division of *self* by *scalarValue*.

        :scalarValue:           Division scalar.

        :returns:               A new instance of :py:class:`Vector`.
        """

        checkScalarValidity(scalarValue, 'division')
        with numpy.errstate(divide='raise'):
            return Vector(values=self.vector/scalarValue)

    def __itruediv__(self, scalarValue):
        """
        In-place true division. See :py:func:`__truediv__`.

        :scalarValue:           Division scalar.

        :returns:               This method returns *self*.
        """

        self = self.__truediv__(scalarValue)

        return self

    def __getitem__(self, index):
        """
        This method return the value of *self* at index *index*.

        :index:         Index of the value to return.

        :returns:       A numpy.float64 value.
        """

        return self.vector[index]

    def __setitem__(self, index, newScalarValue):
        """
        Replace value at specified location of self.vector.

        :index: Index of the self.vector value to replace.
        :newScalarValue: Replacement for value of self.vector at given location.
        """
        checkScalarValidity(newScalarValue, 'item assignment')

        self.vector[index] = numpy.float64(newScalarValue)
        return self

    def __len__(self):
        """
        This method returns the number of values in *self*.

        :returns:       A python int.
        """

        return len(self.vector)

    def __str__(self):
        """
        This method returns a string representation of *self*'s values.

        :returns:       A python str.
        """

        return str(self.vector)

    @property
    def sum(self):
        """
        This method returns the sum of the *self*'s values.

        :returns:       A python float.
        """

        return self.vector.sum()

    @property
    def isSorted(self):
        """
        This methoe returns True if the values of *self* are sorted and False otherwise.

        :returns:       A python boolean.
        """

        return all(numpy.sort(self.vector) == self.vector)
    
    @property
    def size(self):
        """
        This method returns number of values in *self*.

        :returns:       A python int.
        """

        return self.vector.size

    def reverse(self):
        """
        This method reverses the ordering of the value in *self*.
        """

        def copy(self):
            """
            Return a copy of self. Why is this here?
            """

            return Vector(values=self.vector)

        self.vector = self.vector[::-1]

def checkScalarValidity(scalarValue, operation):
    """
    This function check if *scalarValue* is a valid number.

    :scalarValue:       The scalar to test.
    :operation:         A string with the name of the calling function used in the Exception string.
    """

    if not isinstance(scalarValue, (int, float)):
        raise TypeError(f'Only {operation} with scalar value is allowed.')

    elif numpy.isinf(scalarValue):
        raise ValueError(f'Attempting {operation} by {numpy.inf} scalar value.')

    elif numpy.isnan(scalarValue):
        raise ValueError(f'Attempting {operation} by {numpy.nan} scalar value.')
