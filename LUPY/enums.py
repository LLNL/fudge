# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Module to store base enum class used by other class enums in **FUDGE**.
"""

import enum

def auto():
    """
    Function to wrap the :py:func:`enum.auto` function in the python **enum** module.

    :retuns:        The returned value of :py:func:`enum.auto`.
    """

    return enum.auto()

class Enum(enum.Enum):
    """
    This is a base enum class used within **FUDGE**. If :py:func:`enums.auto` is called for a enum, 
    its value is assigned from its name per the :py:func:`_generate_next_value_` method.
    """

    def _generate_next_value_(name, start, count, last_values):
        """
        For internal use per :py:class:`enum.Enum` design. That is, an :py:class:`enum.Enum` hook that allows the 
        :py:func:`enum.auto` to assign the name of the enum to its value.
        """

        return name

    def __str__(self):
        """Returns a string representation of *self* which is its *value* member."""
            
        return self.value

    @classmethod
    def fromString(cls, string):
        """
        Returns a enum for cls from the specified string representation.

        :param cls:             An :py:class:`Enum` derived class.
        :param string:          A string representing one of the allowed values of *cls*.

        :returns:               The enum of *cls* with value *string*.

        :raises TypeError:      If *string* is not a python str instance.
        :raises ValueError:     If no match is found in *cls*.
        """
    
        if not isinstance(string, str):
            raise TypeError('Argument must be a string not "%s".' % type(string))
            
        for entry in cls:
            if entry.value == string:
                return entry

        raise ValueError('Invalid string value "%s" for enum "%s.%s".' % (string, cls.__module__, cls.__name__))

    @classmethod
    def checkEnumOrString(cls, value):
        """
        Checks that *value* is a valid enum representation of *cls*.

        :param cls:             An :py:class:`Enum` derived class.
        :param value:           Can be a string representing an enum of *cls* and an enum of *cls*.

        :returns:               The enum of *cls* represented by *value*.

        :raises TypeError:      If *string* is not a python str instance.
        """

        if isinstance(value, str):
            value = cls.fromString(value)
        if value not in cls:
            raise TypeError('Invalid enum "%s" for class "%s".' % (type(value), type(cls)))

        return value
