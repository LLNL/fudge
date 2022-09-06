# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Module to store base enum classes.
"""

import enum

def auto():

    return enum.auto()

class Enum(enum.Enum):
    """
    This is a base enum class used within **FUDGE**. If enum.auto() is called for a enum, its value is assigned from its name.
    """

    def _generate_next_value_(name, start, count, last_values):
        """
        For internal use per **enum.Enum** design. That is, an enum.Enum hook that allows the enum.auto() to assign
        the name of the enum to its value.
        """

        return name

    def __str__(self):
        """Returns a string representation of *self*."""
            
        return self.value

    @classmethod
    def fromString(cls, string):
        """Returns a enum for cls from the specified string representation."""
    
        if not isinstance(string, str):
            raise TypeError('Argument must be a string not "%s".' % type(string))
            
        for entry in cls:
            if entry.value == string:
                return entry

        raise ValueError('Invalid string value "%s" for enum "%s.%s".' % (string, cls.__module__, cls.__name__))

    @classmethod
    def checkEnumOrString(cls, value):

        if isinstance(value, str):
            value = cls.fromString(value)
        if value not in cls:
            raise TypeError('Invalid enum "%s" for class "%s".' % (type(value), type(cls)))

        return value
