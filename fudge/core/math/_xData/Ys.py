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
This module contains the class Ys. This class stores a list of float data called the y-values and
axes information for the y-values. A y-value can be any python object that can be coerced into a
float. For example, the following 

    ys0 = Ys( axes.defaultAxes( 1 ), data = [ 1, 2., '3.', '4', '5e5' ] )

produces the list containing the floats 1.0, 2.0, 3.0, 4.0 and 500000.0.

Relevant Ys Members:
    name            Description
----------------------------------------------------------------------------
    ys              The list of y-values which must all be a number that is 
                    coerced into a float.
    axes            An axes class of dimension 1, describing the y data 
                        (e.g., label, interpolation, units).
"""

import axes, XYs
from xData.ancestry import ancestry

monikerYs = 'Ys'

class Ys( ancestry ) :

    type = monikerYs
    xData = monikerYs
    moniker = monikerYs

    def __init__( self, axes_, data = None, index = None, value = None, parent = None, isPrimaryXData = None ) :

        attribute = None
        if( parent is not None ) :
            attribute = 'index'
            if( index is None ) : raise Exception( 'index argument cannot be None when parent is not None' )
        ancestry.__init__( self )

        axes.isValidAxes( axes_ )
        axes_.checkDimension( 1 )

        if( data is None ) : data = []
        self.ys = [ float( y ) for y in data ]
        self.axes = axes_.copy( )           # Even secondary must have axes
        self.index = index
        self.value = value
        if( isPrimaryXData is not None ) : self.isPrimaryXData = isPrimaryXData

    def __len__( self ) :
        """Returns the number of y-values."""

        return( len( self.ys ) )

    def __getitem__( self, index ) :
        """Returns the y-value at index."""

        return( self.ys[index] )

    def __setitem__( self, index, y ) :
        """Sets the y-value an index to y."""

        self.ys[index] = float( y )

    def __abs__( self ) :
        """Returns an instance for which each y-value is the absolution of each y-value of self."""

        absed = Ys( self.axes )
        for y in self : absed.ys.append( abs( y ) )
        return( absed )

    def __eq__( self, other ) :
        """Returns True if each element of self is equal to its corresponding element in other, and False otherwise."""

        if( len( self ) != len( other ) ) : return( False )
        for i1, y in enumerate( self ) :
            if( y != float( other[i1] ) ) : return( False )
        return( True )

    def __neg__( self ) :
        """Returns an instance for which each y-value is the negation of each y-value of self."""

        negged = Ys( self.axes )
        for y in self : negged.ys.append( -y )
        return( negged )

    def __iadd__( self, other ) :
        """Adds other to self. Other can be a number, or an instance of length 0 or len( self ) containing numbers."""

        try :
            v = float( other )
            for i1 in xrange( len( self ) ) : self.ys[i1] += v
        except :
            n1 = len( self )
            n2 = len( other )
            if( n2 == 0 ) : return( self )
            if( n1 != n2 ) :
                if( n1 != 0 ) : raise Exception( "len( self ) = %d != len( other ) = %d" % ( n1, n2 ) )
                raise Exception( "Other of type '%s' does not support interating" % type( other ) )
            ys = self.ys
            if( n1 == 0 ) : ys = n2 * [ 0. ]
            for i1 in xrange( n2 ) : ys[i1] += float( other[i1] )
            self.ys = ys
        return( self )

    def __add__( self, other ) :
        """Returns a new instances that is the sum of self and other (also see __iadd__)."""

        added = self.copy( )
        return( added.__iadd__( other ) )

    def __radd__( self, other ) :
        """See __add__."""

        return( self.__add__( other ) )

    def __isub__( self, other ) :
        """Substracts other from self. Other can be a number, or an instance of length 0 or len( self ) containing numbers."""

        try :
            v = float( other )
            for i1 in xrange( len( self ) ) : self.ys[i1] -= v
        except :
            n1 = len( self )
            n2 = len( other )
            if( n2 == 0 ) : return( self )
            if( n1 != n2 ) :
                if( n1 != 0 ) : raise Exception( "len( self ) = %d != len( other ) = %d" % ( n1, n2 ) )
                raise Exception( "Other of type '%s' does not support interating" % type( other ) )
            ys = self.ys
            if( n1 == 0 ) : ys = n2 * [ 0. ]
            for i1 in xrange( n2 ) : ys[i1] -= float( other[i1] )
            self.ys = ys
        return( self )

    def __sub__( self, other ) :
        """Returns a new instances that is the difference between self and other (also see __isub__)."""

        subbed = self.copy( )
        return( subbed.__isub__( other ) )

    def __rsub__( self, other ) :
        """See __sub__."""

        return( self.__sub__( other ).__neg__( ) )

    def __imul__( self, other ) :
        """Return self where each element of self is multiplied by other. Other must be an object that can be coerced into a float."""

        v = float( other )
        for i1 in xrange( len( self ) ) : self.ys[i1] *= v
        return( self )

    def __mul__( self, other ) :
        """Returns a new instance that is each element of self times other. (also see __rmul__)."""

        mulled = self.copy( )
        return( mulled.__imul__( other ) )

    def __rmul__( self, other ) :
        """See __mul__."""

        return( self.__mul__( other ) )

    def __idiv__( self, other ) :
        """Return self where each element of self is divided by other. Other must be an object that can be coerced into a float."""

        v = float( other )
        for i1 in xrange( len( self ) ) : self.ys[i1] /= v
        return( self )

    def __div__( self, other ) :
        """Returns a new instance that is of each element of self divided by other. (also see __rdiv__)."""

        dived = self.copy( )
        return( dived.__idiv__( other ) )

    def __rdiv__( self, other ) :
        """See __div__."""

        return( self.__div__( other ) )

    def append( self, y ) :
        """Appends y to the end of the list of y-values. Y can be any thing that can be coerced into a float."""

        self[ len( self ) ] = float( y )

    def copy( self, index = None, value = None, parent = None, axes = None, isPrimaryXData = None ) :
        """Returns a new Ys instance that is a copy of self."""

        if( index is None ) : index = self.index
        if( value is None ) : value = self.value
        if( axes is None ) : axes = self.axes
        n = Ys( axes, index = index, value = value, parent = parent, isPrimaryXData = isPrimaryXData )
        for y in self.ys : n.ys.append( y )
        return( n )

    def toString( self, format = "%15.8e", sep = " " ) :
        """Returns a string representation of self. The separator between y-values is 'sep'. Each y-value is converted to a string using 'format'."""

        s = [ format % y for y in self ]
        return( sep.join( s ) )

    def plot( self, indexMin = None, indexMax = None ) :
        """Plot the y-values versus their index."""

        if( indexMin is None ) : indexMin = 0
        if( indexMax is None ) : indexMax = len( self )
        xys = [ [ i1, y ] for i1, y in enumerate( self ) ]
        axesXYs = axes.defaultAxes( )
        axesXYs[1] = self.axes[0]
        xys = XYs.XYs( axesXYs, xys, 1e-3 )
        xys.plot( )

    @staticmethod
    def defaultAxes( labelsUnits = {} ) :

        return( axes.defaultAxes( dimension = 1, labelsUnits = labelsUnits ) )

if( __name__ == "__main__" ) :

    axesYs = axes.defaultAxes( 1 )
    axesYs[0] = axes.axis( 'xSec', 0, 'b' )

    ys0 = Ys( axesYs, data = [ 1, 2., '3.', '4', '5e5' ] )
    print '   ', ys0.toString( format = "%.1f" )

    print 'Negation testing'
    ys1 = -ys0
    print '   ', ys0.toString( format = "% .1f" )
    print '   ', ys1.toString( format = "% .1f" )
    ys1 = -ys1
    print '   ', ys1.toString( format = "% .1f" )

    print '\nAbsolution testing'
    for i1, y in enumerate( ys1 ) : ys1[i1] = (-1)**i1 * y
    print '   ', ys1.toString( format = "% .1f" )
    ys2 = abs( ys1 )
    print '   ', ys1.toString( format = "% .1f" )
    print '   ', ys2.toString( format = "% .1f" )

    print '\nMin/Max testing'
    print '   ', ys1.toString( format = "% .1f" )
    print '   ', min( ys1 ), max( ys1 )

    print '\nPlot testing'
    ys1 = Ys( axesYs, data = [ 1, 2. / 7., 4, 7. / 1.01 ] )
    print '   ', ys1.toString( )
    print '   ', ys1.toString( sep = ', ', format = "%.3f" )
#    ys1.plot( )
    print '\nAddition testing'
    ys2 = Ys( axesYs, data = range( 4 ) )
    print '   ', ys2.toString( sep = ', ', format = "%.3f" )

    print 'ys2 += 2.: __iadd__'
    ys2 += 2.
    print '   ', ys2.toString( sep = ', ', format = "%.3f" )
    ys2 += '2.5e2'
    print '   ', ys2.toString( sep = ', ', format = "%.3f" )

    print 'ys1 + ys2: __add__'
    ys3 = ys1 + ys2
    print '   ', ys1.toString( sep = ', ', format = "%.3f" )
    print '   ', ys2.toString( sep = ', ', format = "%.3f" )
    print '   ', ys3.toString( sep = ', ', format = "%.3f" )

    print 'ys4 = -256 + ys2: testing __radd__'
    ys4 = -256 + ys2
    print '   ', ys4.toString( sep = ', ', format = "%.3f" )

    print '\nSubstraction testing'
    print 'ys4 = ys2 - -256: __sub__'
    print '   ', ys4.toString( sep = ', ', format = "%.3f" )

    print 'ys4 = 256 - ys2: __rsub__'
    ys4 = 256 - ys2
    print '   ', ys4.toString( sep = ', ', format = "%.3f" )
    ys4 = '256' - ys2
    print '   ', ys4.toString( sep = ', ', format = "%.3f" )

    print 'ys3 = ys3 - ys2'
    ys3 = ys3 - ys2
    print '   ', ys1.toString( sep = ', ', format = "%.3f" )
    print '   ', ys3.toString( sep = ', ', format = "%.3f" )

    print "ys2 -= '2.5e2'"
    print '   ', ys2.toString( sep = ', ', format = "%.3f" )
    ys2 -= '2.5e2'
    print '   ', ys2.toString( sep = ', ', format = "%.3f" )
    ys2 -= 2
    print '   ', ys2.toString( sep = ', ', format = "%.3f" )

    print '\nMultiplication testing'
    print '   ', ys2.toString( sep = ', ', format = "%.3f" )
    print "ys2 *= '3'"
    ys2 *= '3'
    print '   ', ys2.toString( sep = ', ', format = "%.3f" )
    print "ys3 = y2 * 4"
    ys3 = ys2 * 4
    print '   ', ys3.toString( sep = ', ', format = "%.3f" )
    print "ys3 = 5 * y2"
    ys3 = 5 * ys2
    print '   ', ys3.toString( sep = ', ', format = "%.3f" )

    print '\nDivision testing'
    print '   ', ys2.toString( sep = ', ', format = "%.3f" )
    print "ys2 /= '3'"
    ys2 /= '3'
    print '   ', ys2.toString( sep = ', ', format = "%.3f" )
    print "ys3 = y2 / 4"
    ys3 = ys2 / 4
    print '   ', ys3.toString( sep = ', ', format = "%.3f" )
    print "ys3 = 5 / y2"
    ys3 = 5 / ys2
    print '   ', ys3.toString( sep = ', ', format = "%.3f" )

    print '\nEquality testing'
    data = [ 1, 2., '3.', '4', '5e5' ]
    ys0 = Ys( axesYs, data = data )
    print '   ', ys0 == data, ys0 == data[:4]
    data.append( 6 )
    print '   ', ys0 == data, ys0 == data[:5]
    data[1] += 3
    print '   ', ys0 == data, ys0 == data[:5]
