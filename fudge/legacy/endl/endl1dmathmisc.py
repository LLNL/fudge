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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
