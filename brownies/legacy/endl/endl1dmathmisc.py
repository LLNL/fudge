# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains miscellaneous routines for supporting the endl1dmath class from module endl1dmathClasses.

Module's global variables::

    endl1d_repr_xFormat : Default format string used for converting an element of the data of
                          an endl1dmath object to a string.
"""

from fudge.core.math import fudgemath

endl1d_repr_xFormat = "%14.7e"

def get1dmathData( object, callingRoutine, msg ) :
    """If object is a python list, then it is returned. Else, if object is a sub-class
    of endl1dmath its data is return otherwise a raise is executed."""

    if( type( object ) == type( [] ) ) :
        data = object
    else :
        data = valid1dClassType( object, callingRoutine, msg ).data
    return( data )

def valid1dClassType( object, callingRoutine = "", msg = "" ) :
    """Returns the first argument, object, if it is a subclass of the endl1dmath class; else, triggers a raise."""

    from brownies.legacy.endl import endl1dmathClasses
    if( isinstance(object, endl1dmathClasses.endl1dmath)) : return(object)
    raise Exception( "\nError in %s: invalid type = %s for %s" % ( callingRoutine, type( object ), msg ) )

def check1dData( data ) :
    """Checks that data is a python list of numbers. Triggers a raise if it is not.  Data may be a
python list or an instance that is a subclass of endl1dmath."""

    data = get1dmathData( data, "check1dData", "" )
    i = 0
    for p in data :
        if( not( fudgemath.isNumber( p ) ) ) : raise Exception( "\nError in check1dData: data at index %d is not a number." % i )
        i += 1
