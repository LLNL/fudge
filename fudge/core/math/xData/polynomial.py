# <<BEGIN-copyright>>
# <<END-copyright>>

from fudge.core.ancestry import ancestry
import axes, XYs

__metaclass__ = type

monikerPolynomial = 'polynomial'

class polynomial( ancestry ) :

    moniker = monikerPolynomial
    tag = monikerPolynomial

    def __init__( self, axes, coefficients, parent = None ) :

        axes.checkDimension( 2 )
        ancestry.__init__( self, self.moniker, parent )
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

    def getValue( self, x ) :

        x, P = float( x ), 0.
        for c_o in reversed( self.coefficients ) : P = c_o + P * x
        return( P )

    def toPointwiseLinear( self, xMin, xMax, accuracy, biSectionMax = 8 ) :
        """Calculates the self's y-value at various x-values and the returns the results as a XYs.pointwiseXY instance."""

        def func( x, dummy ) : return( self.getValue( x ) )

        accuracy = min( .1, max( 1e-8, accuracy ) )
        if( len( self ) == 0 ) :
            return( None )
        elif( len( self ) < 3 ) :
            data = [ [ xMin, self.getValue( xMin ) ], [ xMax, self.getValue( xMax ) ] ]
        else :
            n, data = 4 * len( self ), []
            dx = ( xMax - xMin ) / n
            for i in xrange( 4 * len( self ) ) :
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
    def parseXMLNode( cls, polyElement, linkData={} ):
        axes_ = axes.parseXMLNode( polyElement[0] )
        coefficients = map(float, polyElement[1].text.split())
        return cls(axes_, coefficients)

if( __name__ == '__main__' ) :

    axes_ = axes.axes( )
    axes_[0] = axes.axis( 'energy_in', 0, 'eV', frame = axes.labToken, interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
    axes_[1] = axes.axis( 'multiplicity', 0, '', frame = axes.labToken )
    poly = polynomial( axes_, [ 0.5, 0.3, -0.2, 0.1, 0.01 ] )
    pw = poly.toPointwiseLinear( .1, 11, 1e-3 )
    if( True ) :
        print len( poly )
        print poly
        print 'f(.5) =', poly.getValue( 0.5 )
        print poly.toXML( )
        print
        print pw.toString( )
# f(x) = 0.5 + x *( 0.3 + x * ( -0.2 + x * ( 0.1 + x * 0.01 ) ) )
    poly = polynomial( axes_, [ -2520., 3132., -1418., 293., -28., 1.] )
    pw = poly.toPointwiseLinear( .1, 11, 1e-3 )
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
        pw2 = poly.toPointwiseLinear( .1, 11, 1e-5, biSectionMax = 12 )
        print len( pw ), len( pw2 )
        pw.setSafeDivide( 1 )
        r = pw / pw2 - 1
        r.plot( )
# f(x) = -2520. + x * ( 3132. + x * ( -1418. + x * ( 293. + x * ( -28. + x  ) ) ) ) # or ( x - 2 ) * ( x - 3 ) * ( x - 6 ) * ( x - 7 ) * ( x - 10 )
