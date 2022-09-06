# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import numpy

class ProductArray:

    def __init__(self, array=None):
        """
        Constructor. If the dimension of array is 2, then it is converted to a 3 dimensional array.

        :param array: A 2 or 3 dimensional numpy array containing the data.
        """

        if array is None:
            array = numpy.array([])
        elif not isinstance(array, numpy.ndarray):
            raise TypeError('Invalid array of type "%s".' % type(array))

        if array.ndim == 2:
            array.reshape(array.shape + (1,))

        if array.size != 0:
            if array.ndim != 3:
                raise Exception('Invalid length of %d of input array.' % len(array))

        self.__array = array

    def __len__(self):
        """Returns the size of the 3rd dimension of the array."""

        return len(self.__array)

    def __add__(self, other):
        """Returns a new **ProductArray** instance that is the sum of *self* and *other*."""

        if not isinstance(other, ProductArray):
            raise TypeError('Supported other of type "%s".' % type(other))

        if len(other) == 0:
            return ProductArray(self.__array)
        if len(self) == 0:
            return ProductArray(other.array)

        if self.__array.shape[:2] != other.array.shape[:2]:
            raise Exception('Shape of self (%s) and other (%s) not compatible.' % (self.__array.shape, other.array.shape))

        smaller = other.array
        larger = self.__array
        if smaller.shape[2] > larger.shape[2]:
            smaller, larget = larger, smaller

        n3 = smaller.shape[2]
        if n3 == larger.shape[2]:
            array = larger + smaller
        else:
            array = larger
            array[:,:,0:n3] += smaller

        return ProductArray(array)

    @property
    def array(self):
        """Returns a reference to data of self."""

        return self.__array

    @property
    def shape(self):
        """Returns the shape of the internally stored numpy array."""

        return self.__array.shape
