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

import math 
from fudge.core.math import fudgemath

__metaclass__ = type

class equalProbableBinnedData :

    def __init__( self, data ) :
        """Constructor for the equalProbableBinnedData class."""

        self.data = data

    def __getitem__( self, i ) :
        """Returns the (i+1)^th element of self."""

        return( self.data[i] )

    def __setitem__( self, i, value ) :
        """Sets the (i+1)^th element of self to value."""

        self.data[i] = value

    def __len__( self ) :
        """Returns the number of data points in self."""

        return( len( self.data ) )

    def getData( self ) :
        """Returns the data of self."""

        return( self.data )

    def getEnd( self ) :
        """Returns the number of data points in self."""

        return( len( self ) )

    def getStart( self ) :
        """Always returns 0."""

        return( 0 )

def equalProbableBins( nBins, xy ) :

    S = fudgemath.runningYSum( xy )
    SMax = S.data[-1][1]
    dS = SMax / nBins
    i = 0
    yp1 = xy[i][1]
    x1, y1 = S.data[0]
    iSum = 0
    runningSum = 0.
    epbs = [ x1 ]
    xMax = max( abs( xy[0][0] ), abs( xy[-1][0] ) )
    for x2, y2 in S.data :
        yp2 = xy[i][1]
        exitLoop = False
        if( x1 != x2 ) :
            while( ( runningSum + dS ) <= y2 ) :
                iSum += 1
                if( iSum == nBins ) :
                    exitLoop = True
                    break
                runningSum = SMax * ( float( iSum ) / float( nBins ) )
                c = runningSum - y1
                if( c == 0. ) : continue
                if( yp1 == yp2 ) :
                    if( yp1 != 0. ) :
                        x = x1 + c / yp1
                    else :                                                      # Both yp1 and yp2 are zero.
                        x = None
                else :
                    a = ( yp2 - yp1 ) / ( x2 - x1 )
                    b = yp1
                    sqrtArgument = b * b + 2. * a * c
                    if( sqrtArgument < 0. ) :
                        if(b * b * 1e-12 < -sqrtArgument) : raise Exception( 'equalProbableBins: b^2 + 2 a c  = %e < 0. a = %e, b = %e c = %e' % ( sqrtArgument, a, b, c ) )
                        sqrtArgument = 0.
                    x = x1 + 2. * c / ( math.sqrt( sqrtArgument ) + b ) # c (a) should be -c (a/2)
                if( x != None ) :
                    if( ( abs( x ) < 3e-16 ) and ( 1e-8 * xMax > abs( x ) ) ) : x = 0.  # Special case for when x should probably be 0.
                    epbs.append( x )
        if( exitLoop ) : break
        i += 1
        x1 = x2
        y1 = y2
        yp1 = yp2
    if( len( epbs ) != ( nBins + 1 ) ) : epbs.append( S.data[-1][0] )

    i = 0                               # Correct for case where data starts with more than one P(mu) = 0.
    for x, y in xy :
        if( y != 0. ) : break
        i += 1
    i -= 1
    if( i > 0 ) : epbs[0] = xy[i][0]

    i = 0                               # Correct for case where data ends with more than one P(mu) = 0.
    iEnd = 0
    for x, y in xy :
        if( y != 0. ) : iEnd = i
        i += 1
    iEnd += 1
    if( iEnd < len( xy ) ) : epbs[-1] = xy[iEnd][0]

    return( equalProbableBinnedData( epbs ) )
