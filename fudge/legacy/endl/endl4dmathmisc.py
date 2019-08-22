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
This module contains miscellaneous routines for supporting the endl4dmath class from module endl4dmathClasses.

Module's global variables::

    endl4d_repr_tFormat : Default format string used for converting an t-element of the data of an endl4dmath object to a string.
    endl4d_repr_xFormat : Default format string used for converting an x-element of the data of an endl4dmath object to a string.
    endl4d_repr_yFormat : Default format string used for converting an y-element of the data of an endl4dmath object to a string.
    endl4d_repr_zFormat : Default format string used for converting an z-element of the data of an endl4dmath object to a string.
"""

import os
import types

from fudge.core import fudgemisc
from fudge.core.utilities import fudgeFileMisc, subprocessing
from fudge.core.math import fudgemath
import endl4dmathClasses

endl4d_repr_tFormat = "%14.7e"
endl4d_repr_xFormat = "%12.5e"
endl4d_repr_yFormat = "%12.5e"
endl4d_repr_zFormat = "%12.5e"

def get4dmathData( object, callingRoutine, msg ) :
    """If object is a python list, then it is returned. Else, if object is a sub-class
    of endl4dmath its data is return otherwise a raise is executed."""

    if( type( object ) == type( [] ) ) :
        data = object
    else :
        data = valid4dClassType( object, callingRoutine, msg ).data
    return( data )

def valid4dClassType( object, callingRoutine, msg ) :
    "Returns the first argument, object, if it is a subclass of the endl4dmath class; else, triggers a raise."

    if( isinstance( object, endl4dmathClasses.endl4dmath ) ) : return( object )
    raise Exception( "\nError in %s: invalid type = %s for %s" % ( callingRoutine, type( object ), msg ) )

def check4dData( data, allowNegativeT = False, allowZeroT = True, allowNegativeX = False, allowZeroX = True, allowSameY = False, allowNegativeY = False,
    allowZeroY = True, positiveZ = False, printWarning = True, printErrors = True, xCloseEps = None, maxAbsFloatValue = None ) :
    """Checks that data is a valid structure for endl4dmath data. Triggers a raise if it is not.  Data may be
    an instance that is a subclass of endl4dmath or a python list as required by endl4dmath. If allowZeroT is False, 
    a zero first t value will generate a raise.If allowNegativeT is true than the t data can have negative values 
    (many ENDL 4d-data have projectile energy as the t-data which must always be positive). For meanings of 
    allowNegativeX, allowZeroX, allowSameY, allowNegativeY, allowZeroY, positiveZ and xCloseEps see endl3dmathmisc.check3dData.
    If printErrors is 'True' then errors are printed, and if at least one error exist, a raise will be executed.  
    printWarning is just passed on to endl3dmathmisc.check3dData."""

    import endl3dmathmisc
    points = get4dmathData( data, "check4dData", "" )
    messages = []
    tPrior = None
    if( ( points[0][0] < 0. ) and ( not allowNegativeT ) ) :
        s = 'check4dData: zero x value at index = 0'
        messages.append( s )
        if( printErrors ) : fudgemisc.printWarning( '\n'.join( fudgemisc.checkMessagesToString( s, indentation = '    ' ) ) )
    elif( ( points[0][0] < 0. ) and ( not allowNegativeT ) ) :
        s = 'check4dData: negative t = %e' % points[0][0]
        messages.append( s )
        if( printErrors ) : fudgemisc.printWarning( '\n'.join( fudgemisc.checkMessagesToString( s, indentation = '    ' ) ) )
    i = 0
    for t, xyz in points :
        fudgemath.checkNumber( t, "check4dData: t[%d]" % i, messages = messages, indentation = '    ', printErrors = printErrors, \
            maxAbsFloatValue = maxAbsFloatValue )
        if( tPrior != None ) :
            if( t <= tPrior ) :
                s = 'check4dData: t value t[%d] >= t[%d] (%e >= %e)' % ( i-1, i, tPrior, t )
                messages.append( s )
                if( printErrors ) : fudgemisc.printWarning( '\n'.join( fudgemisc.checkMessagesToString( s, indentation = '    ' ) ) )
        tPrior = t
        messages3d = endl3dmathmisc.check3dData( xyz, allowNegativeX = allowNegativeX, allowZeroX = allowZeroX, allowSameY = allowSameY, 
            allowNegativeY = allowNegativeY, allowZeroY = allowZeroY, positiveZ = positiveZ, printWarning = printWarning, printErrors = False, 
            xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages3d ) > 0 ) :
            s = 'check4dData: t value = %e' % t
            messages.append( [ s, messages3d ] )
            if( printErrors ) : fudgemisc.printWarning( '\n'.join( fudgemisc.checkMessagesToString( [ messages[-1] ], indentation = '    ' ) ) )
        i += 1
    if( printErrors and len( messages ) > 0 ) : raise Exception( '      bad 4d data' )
    return( messages )

def plot3dFrom4d( data, i, options, xyzlog = 0 ) :
    "For internal use only."

    ly_etal = data[i][1]
    f = fudgeFileMisc.fudgeTempFile( )
    for y_etal in ly_etal :
        y = "%15.7e" % y_etal[0]
        for etal in y_etal[1] : f.write( "%s %14.6e %14.6e\n" % ( y, etal[0], etal[1] ) )
    f.close( )
    p = os.path.join( __file__.split( '/fudge/core/' )[0], "fudge", "vis", "gnuplot", "endl3dplot.py" )
    s = [ "python", p, 'xyzlog', str( xyzlog ) ] + options + [ f.getName( ) ]
    subprocessing.spawn( s )

def string4dData( d, i0 = 0, i1 = 1, i2 = 2, i3 = 3, fmt0 = None, fmt1 = None, fmt2 = None, fmt3 = None ) :
    """Returns a list of strings for data d.  D must be list[ number, list[ number, list[ number, number ] ] ]."""

    i = [ i0, i1, i2, i3 ]
    i.sort( )
    if ( i[0] != 0 ) or ( i[1] != 1 ) or ( i[2] != 2 ) or ( i[3] != 3 ) :
        raise Exception( "\nError in string4dData: %s" % "Invalid i0 = %s, i1 = %s, i2 = %s, and/or i3 = %s" % (`i0`,`i1`,`i2`,`i3`) )
    if( fmt0 == None ) : fmt0 = endl4d_repr_tFormat
    if( fmt1 == None ) : fmt1 = endl4d_repr_xFormat
    if( fmt2 == None ) : fmt2 = endl4d_repr_yFormat
    if( fmt3 == None ) : fmt3 = endl4d_repr_zFormat
    fmt = [ fmt0, fmt1, fmt2, fmt3 ]
    s = []
    sa = [ "", "", "", "" ]
    for d0 in d :
        sa[i0] = fmt[i0] % d0[0]
        for d1 in d0[1] :
            sa[i1] = fmt[i1] % d1[0]
            for d2 in d1[1] :
                sa[i2] = fmt[i2] % d2[0]
                sa[i3] = fmt[i3] % d2[1]
                s.append( ' '.join( sa ) )
    return s

