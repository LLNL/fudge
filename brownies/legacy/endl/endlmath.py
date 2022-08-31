# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
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
    from brownies.legacy.endl import endl3dmathClasses
    return endl3dmathClasses.endl3dmath(d3, xLabel = xLabel, yLabel = yLabel, zLabel = zLabel, checkDataType = 0)

def runningYSum( data, normalize = False ) :
    """Returns the running sum of dx * y (normalized to 1 of normalize is True) as an endl2dmath object. 
    Data must be list of ( x, y )."""

    x1 = None
    runningSum = []
    for xy in data :
        x2 = xy[0]
        y2 = xy[1]
        if ( x1 is None ) :
            Sum = 0.
        else :
            Sum += 0.5 * ( y2 + y1 ) * ( x2 - x1 )
        runningSum.append( [ x2, Sum ] )
        x1 = x2
        y1 = y2
    if( normalize and ( Sum != 0. ) ) :
        for xy in runningSum : xy[1] /= Sum
    from brownies.legacy.endl import endl2dmathClasses
    return endl2dmathClasses.endl2dmath(runningSum, checkDataType = 0)

def ZSum( data ) :
    """Returns the area under the curve z(y) for each x as an endl2dmath object. Data must be list of
    ( x, list of ( y, z ) )."""

    d2 = []
    for x_yz in data : d2.append( [ x_yz[0], YSum( x_yz[1] ) ] )
    from brownies.legacy.endl import endl2dmathClasses
    return endl2dmathClasses.endl2dmath(d2, checkDataType = 0)

def YSum( data ) :
    """Returns the area under the curve y(x). Data must be list of list( x, y )."""

    x1 = None
    for x2, y2 in data :
        if ( x1 is None ) :
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
        for i in range( n ) :
            if( self.list[i] is None ) :
                self.list[i] = addend
                addend = None
                break
            else :
                addend = addend + self.list[i]
                self.list[i] = None
        if( addend is not None ) : self.list.append( addend )

    def clear( self ) :
        """Clears currently summed data."""

        self.list = []

    def returnSum( self ) :
        """Returns the current sum of all addends appended."""

        s = None
        for l in self.list :
            if( l is not None ) :
                if( s is None ) :
                    s = l
                else :
                    s = s + l
        return( s )

def getValue( n ) :

    if( isNumber( n ) ) : return( n )
    if( isinstance( n, PQU.PQU ) ) : return( n.getValue( ) )
    raise Exception( 'Invalue number object = %s' % brb.getType( n ) )
