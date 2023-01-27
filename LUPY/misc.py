# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

'''
This module contains miscellaneous functions that do not fit in other modules.
'''

def isString(value, message):
    '''
    This function check if *value* is a type str instance. If it is, *value* is returned. Otherwise, a TypeError raise
    is execute with message *message*.
    '''

    if isinstance(value, str):
        return value

    raise TypeError(message)
