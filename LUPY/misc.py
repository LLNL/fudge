# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

'''
This module contains general miscellaneous functions that do not fit in other LUPY modules.
'''

def isString(value, message):
    '''
    This function check if *value* is of type python str instance. If it is, *value* is returned. Otherwise, a TypeError raise
    is execute with message *message*.

    :param value:           The instance to check if it is a python str instance.
    :param message:         The message passed to **TypeError**.

    :returns:               *value*.

    :raises TypeError:      If *value* is not a python str instance.
    '''

    if isinstance(value, str):
        return value

    raise TypeError(message)
