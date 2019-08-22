# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
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

    return( equalProbableBinnedData( epbs ) )
