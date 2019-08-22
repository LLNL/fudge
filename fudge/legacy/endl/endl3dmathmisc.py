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
This module contains miscellaneous routines for supporting the endl3dmath class from module endl3dmathClasses.

Module's global variables::

    endl3d_repr_xFormat : Default format string used for converting an x-element of the data of an endl3dmath object to a string.
    endl3d_repr_yFormat : Default format string used for converting an y-element of the data of an endl3dmath object to a string.
    endl3d_repr_zFormat : Default format string used for converting an z-element of the data of an endl3dmath object to a string.
"""

import types
from fudge.core import fudgemisc
from fudge.core.math import fudgemath
import endl2dmathClasses
import endl2dmathmisc

endl3d_repr_xFormat = "%14.7e"
endl3d_repr_yFormat = "%12.5e"
endl3d_repr_zFormat = "%12.5e"
noExtrapolation = 'noExtrapolation'
flatExtrapolation = 'flatExtrapolation'

def get3dmathData( object, callingRoutine, msg ) :
    """If object is a python list, then it is returned. Else, if object is a sub-class
    of endl3dmath its data is return otherwise a raise is executed."""


    if( type( object ) == type( [] ) ) :
        data = object
    else :
        data = valid3dClassType( object, callingRoutine, msg ).data
    return( data )

def valid3dClassType( object, callingRoutine, msg ) :
    "Returns the first argument, object, if it is a subclass of the endl3dmath class; else, triggers a raise."
    import endl3dmathClasses
    if( isinstance( object, endl3dmathClasses.endl3dmath ) ) : return( object )
    raise Exception( "\nError in %s: invalid type = %s for %s" % ( callingRoutine, type( object ), msg ) )

def check3dData( data, allowNegativeX = False, allowZeroX = True, allowNegativeY = False, allowSameY = False, allowZeroY = True, positiveZ = False, 
    printWarning = True, printErrors = True, xCloseEps = None, maxAbsFloatValue = None ) :
    """Checks that data is a valid structure for endl3dmath data. Triggers a raise if it is not.  Data may be
    an instance that is a subclass of endl3dmath or a python list as required by endl3dmath. If allowZeroX is 
    False, a zero first x value will generate a raise. If allowNegativeX is true than the x data can have negative 
    values (ENDL 3d-data has projectile energy as the x-data which must always be > 0. ). See endl2dmathmisc.check2dData's 
    allowSameX, allowZeroX, positiveY and xCloseEps argument for meaning of allowSameY, allowZeroY, positiveZ and 
    xCloseEps respectively.  If printErrors is 'True' then errors are printed, and if at least one error exist, a 
    raise will be executed.  printWarning is just passed on to endl2dmathmisc.check2dData."""

    points = get3dmathData( data, "check3dData", "" )
    messages = []
    xPrior = None
    if( ( points[0][0] == 0. ) and ( not allowZeroX ) ) :
            s = 'check3dData: zero x value at index = 0'
            messages.append( s )
            if( printErrors ) : fudgemisc.printWarning( '\n'.join( fudgemisc.checkMessagesToString( s, indentation = '      ' ) ) )
    elif( ( points[0][0] < 0. ) and ( not allowNegativeX ) ) :
            s = 'check3dData: negative x = %e' % points[0][0]
            messages.append( s )
            if( printErrors ) : fudgemisc.printWarning( '\n'.join( fudgemisc.checkMessagesToString( s, indentation = '      ' ) ) )
    i = 0
    for x, yz in points :
        if( fudgemath.isNumber( x ) ) :
            fudgemath.checkNumber( x, "check3dData: x[%d]" % i, messages = messages, indentation = '      ', printErrors = printErrors, \
                maxAbsFloatValue = maxAbsFloatValue )
        else :
            s = 'check3dData: x = %s is not a number' % points[0][0]
            messages.append( s )
            if( printErrors ) : fudgemisc.printWarning( '\n'.join( fudgemisc.checkMessagesToString( s, indentation = '      ' ) ) )
        if( xPrior != None ) :
            if( x <= xPrior ) :
                s = 'check3dData: x value x[%d] >= x[%d] (%e >= %e)' % ( i-1, i, xPrior, x )
                messages.append( s )
                if( printErrors ) : fudgemisc.printWarning( '\n'.join( fudgemisc.checkMessagesToString( s, indentation = '      ' ) ) )
        ne, badXIndicies, messages2d = endl2dmathmisc.check2dData( yz, allowNegativeX = allowNegativeY, allowSameX = allowSameY, allowZeroX = allowZeroY, \
            positiveY = positiveZ, printWarning = printWarning, printErrors = False, xCloseEps = xCloseEps, maxAbsFloatValue = maxAbsFloatValue )
        if( len( messages2d ) > 0 ) :
            s = 'check3dData x value = %e' % x
            messages.append( [ s, messages2d ] )
            if( printErrors ) : fudgemisc.printWarning( '\n'.join( fudgemisc.checkMessagesToString( s, indentation = '      ' ) ) )
        xPrior = x
        i += 1
    return( messages )

def interpolate_XY( data, x, y ) :
    """For internal use only."""

    l = len( data )
    if( l < 1 ) : return( 0. )
    z = 0.
    if( data[0][0] <= x <= data[-1][0] ) :
        i = 0
        while( i < l ) :
            if( x <= data[i][0] ) : break
            i += 1
        if( x == data[i][0] ) :
            z = endl2dmathmisc.interpolate_X( data[i][1], y )
        else :
            d1 = data[i-1]
            d2 = data[i]
            z1 = endl2dmathmisc.interpolate_X( d1[1], y )
            z2 = endl2dmathmisc.interpolate_X( d2[1], y )
            x1 = d1[0]
            x2 = d2[0]
            z = ( z1 * ( x2 - x ) + z2 * ( x - x1 ) ) / ( x2 - x1 )
    return( z )

def interpolate3d( x, data, unitBase = False, extrapolation = noExtrapolation, endl2dmathObject = False ) :
    """Returns the interpolation of data (an endl3dmath object or suitable data)
    at x as a list of [y,z] pairs.  If x is outside the domain of data, then [] 
    is returned. If x is one of the x-values of data, then its yz-data is returned.
    Otherwise, for xl < x < xu where xl and xu are consecutive x-values in data 
    linear interpolation is performed on the yz data of xl and xu.  If unitBase is 
    True, then unit base interpolation is performed."""

    points = get3dmathData( data, "endl3dmathmisc.interpolate3d", "" )
    ubiData = []
    if( len( points ) > 0 ) :
        if( x < points[0][0] ) :
            if( extrapolation == flatExtrapolation ) :
                ubiData = [ [ y, z ] for y, z in points[0][1] ]
        elif( points[-1][0] < x ) :
            if( extrapolation == flatExtrapolation ) :
                ubiData = [ [ y, z ] for y, z in points[-1][1] ]
        else :
            i = 0
            for x_, yz in points :
                if( x <= x_  ) : break
                i += 1
            data_u = points[i]
            xu = data_u[0]
            u = data_u[1]
            if( xu == x ) :
                for y, z in u : ubiData.append( [ y, z ] )
            else :                                          # Note, if xl == xu or i is the last element in data we would not be here.
                data_l = points[i-1]
                xl = data_l[0]
                l = data_l[1]
                f = ( xu - x ) / ( xu - xl )
                g = 1. - f
                if( unitBase and ( ( l[0][0] != u[0][0] ) or ( l[-1][0] != u[-1][0] ) ) ) :
                    unitBase = True
                    yl = ( f * l[ 0][0] + g * u[ 0][0] )
                    yu = ( f * l[-1][0] + g * u[-1][0] )
                    s = ( yu - yl ) / ( l[-1][0] - l[0][0] )
                    lp = []
                    y0 = l[0][0]
                    for y, z in l : lp.append( [ s * ( y - y0 ) + yl, z / s ] )
                    l = lp
                    s = ( yu - yl ) / ( u[-1][0] - u[0][0] )
                    up = []
                    y0 = u[0][0]
                    for y, z in u : up.append( [ s * ( y - y0 ) + yl, z / s ] )
                    up[-1][0] = lp[-1][0]                   # Make sure the last x values are the same.
                    u = up
                l = endl2dmathClasses.endl2dmath( data = l, checkDataType = False )
                u = endl2dmathClasses.endl2dmath( data = u, checkDataType = False )
                if( not unitBase ) :
                    l = l.copyData( )
                    u = u.copyData( )
                    if( ( l.data[0][0] < u.data[0][0] ) and ( l.data[0][1] != 0. ) ) :
                        v = ( 1.0 - 1e-6 ) * u.data[0][0]
                        if( v == 0. ) : v = -1e-12
                        u.data.insert( 0, [ v, 0. ] )
                    elif( ( l.data[0][0] > u.data[0][0] ) and ( u.data[0][1] != 0. ) ) :
                        v = ( 1.0 - 1e-6 ) * l.data[0][0]
                        if( v == 0. ) : v = -1e-12
                        l.data.insert( 0, [ v, 0. ] )
                    if( ( l.data[-1][0] < u.data[-1][0] ) and ( l.data[-1][1] != 0. ) ) :
                        v = ( 1.0 + 1e-6 ) * l.data[-1][0]
                        if( v == 0. ) : v = 1e-12
                        l.data.append( [ v, 0. ] )
                    elif( ( l.data[-1][0] > u.data[-1][0] ) and ( u.data[-1][1] != 0. ) ) :
                        v = ( 1.0 + 1e-6 ) * u.data[-1][0]
                        if( v == 0. ) : v = 1e-12
                        u.data.append( [ v, 0. ] )
                ubiData = f * l + g * u
                ubiData = ubiData.data
    if( endl2dmathObject ) : ubiData = endl2dmathClasses.endl2dmath( ubiData )
    return( ubiData )
