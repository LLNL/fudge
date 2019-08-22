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

"""
This module contains useful fudge math routines that do not fit into any other module.
"""

from pqu import PQU
from fudge.core.utilities import brb
try :
    import numpy
    numpyFloat64 = numpy.float64( 1. )
except :
    numpyFloat64 = 1.

__metaclass__ = type

def runningZSum( data, xLabel = None, yLabel = None, zLabel = None, normalize = False ) :
    """Returns the running sum of dy * z (normalized to 1 of normalize is True) for each x as an endl3dmath object.
    Data must be list of ( x, list of ( y, z ) )."""

    d3 = []
    for x_yz in data : d3.append( [ x_yz[0], runningYSum( x_yz[1], normalize = normalize ).data ] )
    from fudge.legacy.endl import endl3dmathClasses
    return endl3dmathClasses.endl3dmath( d3, xLabel = xLabel, yLabel = yLabel, zLabel = zLabel, checkDataType = 0 )

def runningYSum( data, normalize = False ) :
    """Returns the running sum of dx * y (normalized to 1 of normalize is True) as an endl2dmath object. 
    Data must be list of ( x, y )."""

    x1 = None
    runningSum = []
    for xy in data :
        x2 = xy[0]
        y2 = xy[1]
        if ( x1 == None ) :
            Sum = 0.
        else :
            Sum += 0.5 * ( y2 + y1 ) * ( x2 - x1 )
        runningSum.append( [ x2, Sum ] )
        x1 = x2
        y1 = y2
    if( normalize and ( Sum != 0. ) ) :
        for xy in runningSum : xy[1] /= Sum
    from fudge.legacy.endl import endl2dmathClasses
    return endl2dmathClasses.endl2dmath( runningSum, checkDataType = 0 )

def ZSum( data ) :
    """Returns the area under the curve z(y) for each x as an endl2dmath object. Data must be list of
    ( x, list of ( y, z ) )."""

    d2 = []
    for x_yz in data : d2.append( [ x_yz[0], YSum( x_yz[1] ) ] )
    from fudge.legacy.endl import endl2dmathClasses
    return endl2dmathClasses.endl2dmath( d2, checkDataType = 0 )

def YSum( data ) :
    "Returns the area under the curve y(x). Data must be list of list( x, y )."

    x1 = None
    for x2, y2 in data :
        if ( x1 == None ) :
            Sum = 0.
        else :
            Sum += ( y2 + y1 ) * ( x2 - x1 )
        x1 = x2
        y1 = y2
    return 0.5 * Sum

class fastSumOfManyAddends :
    """This class in designed to sum a lot of endl2dmath or fudge2dmath object together efficiently. For example,
    consider the list f2d of 100,000 fudge2dmath objects that are to be summed. One way to do this is as
    s = fudge2dmath( )
    for f in f2d : s = s + f

    In general, this is very inefficient and will take a long time. Using, this class as

    fs = fastSumOfManyAddends( )
    for f in f2d : fs.appendAddend( f )
    s = fs.returnSum( )

    is, in general, much more efficient (i.e., runs a lot faster) and it should never be less efficient.

    While this class was designed for endl2dmath and fudge2dmath objects, it should work for any object
    for which the '+' operation is defined."""

    def __init__( self ) :
        """Constructor for fastSumOfManyAddends."""

        self.clear( )

    def appendAddend( self, addend ) :
        """Adds addend to current sum efficiently."""

        n = len( self.list )
        for i in xrange( n ) :
            if( self.list[i] == None ) :
                self.list[i] = addend
                addend = None
                break
            else :
                addend = addend + self.list[i]
                self.list[i] = None
        if( addend != None ) : self.list.append( addend )

    def clear( self ) :
        """Clears currently summed data."""

        self.list = []

    def returnSum( self ) :
        """Returns the current sum of all addends appended."""

        s = None
        for l in self.list :
            if( l != None ) :
                if( s == None ) :
                    s = l
                else :
                    s = s + l
        return( s )

def getValue( n ) :

    if( isNumber( n ) ) : return( n )
    if( isinstance( n, PQU.PQU ) ) : return( n.getValue( ) )
    raise Exception( 'Invalue number object = %s' % brb.getType( n ) )
