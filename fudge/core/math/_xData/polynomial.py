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

from xData.ancestry import ancestry
from pqu import PQU
import axes, XYs

__metaclass__ = type

monikerPolynomial = 'polynomial'

class polynomial( ancestry ) :

    moniker = monikerPolynomial
    tag = monikerPolynomial

    def __init__( self, axes, coefficients, parent = None ) :

        axes.checkDimension( 2 )
        ancestry.__init__( self )
        self.axes = axes.copy( parent = self )
        self.coefficients = []
        for c_l in coefficients : self.coefficients.append( float( c_l ) )

    def __len__( self ) :

        return( len( self.coefficients ) )

    def __getitem__( self, order ) :

        return( self.coefficients[order] )

    def __str__( self ) :

        return( ' '.join( [ "%g" % c_o for c_o in self.coefficients ] ) )

    def copy( self, parant = None ) :

        return( polynomial( self.axes, self.coefficients, parent = parent ) )

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQU.PQU( '1 ' + self.domainUnit( ) ).getValueAs( unitTo ) )

    def domainUnit( self ) :

        return( self.axes[0].getUnit( ) )

    def getValue( self, x ) :

        x, P = float( x ), 0.
        for c_o in reversed( self.coefficients ) : P = c_o + P * x
        return( P )

    def toPointwise_withLinearXYs( self, xMin, xMax, accuracy, biSectionMax = 8 ) :
        """Calculates the self's y-value at various x-values and the returns the results as a XYs.pointwiseXY instance."""

        def func( x, dummy ) : return( self.getValue( x ) )

        accuracy = min( .1, max( 1e-8, accuracy ) )
        if( len( self ) == 0 ) :
            return( None )
        elif( len( self ) < 3 ) :
            data = [ [ xMin, self.getValue( xMin ) ], [ xMax, self.getValue( xMax ) ] ]
            return( XYs.XYs( self.axes, data, accuracy = accuracy, biSectionMax = biSectionMax ) )
        else :
            n1, data = 4 * len( self ), []
            dx = ( xMax - xMin ) / n1
            for i in xrange( n1 ) :
                x = xMin + dx * i
                data.append( [ x, x ] )
            data.append( [ xMax, xMax ] )
        return( XYs.XYs( self.axes, data, accuracy = accuracy, biSectionMax = biSectionMax ).applyFunction( func, None, checkForRoots = 1 ) )

    def toXML( self, tag = 'xData', indent = '', incrementalIndent = '  ' ) :

        return( '\n'.join( self.toXMLList( tag = tag, indent = indent, incrementalIndent = incrementalIndent ) ) )

    def toXMLList( self, tag = None, indent = '', incrementalIndent = '  ' ) :

        if( tag is None ) :
            if( hasattr( self, 'tag' ) ) :
                tag = self.tag
            else :
                tag = 'xData'
        indent2 = indent + incrementalIndent
        xmlString = [ '%s<%s xData="%s" length="%s">' % ( indent, tag, self.tag, len( self ) ) ]
        xmlString += self.axes.toXMLList( indent = indent2 )
        xmlString.append( indent2 + '<data>' + ' '.join( [ '%g' % c_o for c_o in self.coefficients ] ) + '</data>' )
        xmlString[-1] += '</%s>' % tag
        return( xmlString )

    @classmethod
    def parseXMLNode( cls, polyElement, xPath=[], linkData={} ):
        xPath.append( polyElement.tag )
        axes_ = axes.parseXMLNode( polyElement[0], xPath )
        coefficients = map(float, polyElement[1].text.split())
        poly_ = cls(axes_, coefficients)
        xPath.pop()
        return poly_

if( __name__ == '__main__' ) :

    axes_ = axes.axes( )
    axes_[0] = axes.axis( 'energy_in', 0, 'eV', interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
    axes_[1] = axes.axis( 'multiplicity', 0, '' )
    poly = polynomial( axes_, [ 0.5, 0.3, -0.2, 0.1, 0.01 ] )
    pw = poly.toPointwise_withLinearXYs( .1, 11, 1e-3 )
    if( True ) :
        print len( poly )
        print poly
        print 'f(.5) =', poly.getValue( 0.5 )
        print poly.toXML( )
        print
        print pw.toString( )
# f(x) = 0.5 + x *( 0.3 + x * ( -0.2 + x * ( 0.1 + x * 0.01 ) ) )
    poly = polynomial( axes_, [ -2520., 3132., -1418., 293., -28., 1.] )
    pw = poly.toPointwise_withLinearXYs( .1, 11, 1e-3 )
    if( True ) :
        print
        print len( poly )
        print poly
        print 'f(.5) =', "%25.17e" % poly.getValue( 0.5 )
        print 'f(2) =', "%25.17e" % poly.getValue( 2. )
        print 'f(3) =', "%25.17e" % poly.getValue( 3. )
        print 'f(6) =', "%25.17e" % poly.getValue( 6. )
        print 'f(7) =', "%25.17e" % poly.getValue( 7. )
        print 'f(10) =', "%25.17e" % poly.getValue( 10. )
        print poly.toXML( )
        print
        pw.plot( )
        pw2 = poly.toPointwise_withLinearXYs( .1, 11, 1e-5, biSectionMax = 12 )
        print len( pw ), len( pw2 )
        pw.setSafeDivide( 1 )
        r = pw / pw2 - 1
        r.plot( )
# f(x) = -2520. + x * ( 3132. + x * ( -1418. + x * ( 293. + x * ( -28. + x  ) ) ) ) # or ( x - 2 ) * ( x - 3 ) * ( x - 6 ) * ( x - 7 ) * ( x - 10 )
