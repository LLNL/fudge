# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import math 
from fudge.core.math import fudgemath


class EqualProbableBinnedData :

    def __init__( self, data ) :
        """Constructor for the EqualProbableBinnedData class."""

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
    domainMax = max( abs( xy[0][0] ), abs( xy[-1][0] ) )
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
                if( x is not None ) :
                    if( ( abs( x ) < 3e-16 ) and ( 1e-8 * domainMax > abs( x ) ) ) : x = 0.  # Special case for when x should probably be 0.
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

    return( EqualProbableBinnedData( epbs ) )
