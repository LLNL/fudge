# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

"""
This module contains miscellaneous routines for supporting the endl1dmath class from module endl1dmathClasses.

Module's global variables::

    endl1d_repr_xFormat : Default format string used for converting an element of the data of
                          an endl1dmath object to a string.
"""

import types
import endl1dmathClasses
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

    if( isinstance( object, endl1dmathClasses.endl1dmath ) ) : return( object )
    raise Exception( "\nError in %s: invalid type = %s for %s" % ( callingRoutine, type( object ), msg ) )

def check1dData( data ) :
    """Checks that data is a python list of numbers. Triggers a raise if it is not.  Data may be a
python list or an instance that is a subclass of endl1dmath."""

    data = get1dmathData( data, "check1dData", "" )
    i = 0
    for p in data :
        if( not( fudgemath.isNumber( p ) ) ) : raise Exception( "\nError in check1dData: data at index %d is not a number." % i )
        i += 1
