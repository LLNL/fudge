# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from . import ancestry as ancestryModule
from . import standards as standardsModule
from . import values as valuesModule
from . import base as xDataBaseModule
from . import axes as axesModule
from . import XYs as XYsModule

class data( ancestryModule.ancestry ) :

    def __init__( self, values ) :

        if( not( isinstance( values, valuesModule.values ) ) ) : raise TypeError( 'Not a values type' )
        ancestryModule.ancestry.__init__( self )

        self.values = values

    def __len__( self ) :

        return( len( self.values ) )

    def copy( self ) :

        return( self.__class__( self.values ) )

    def offsetScaleData( self, offset, scale ) :

        self.values.offsetScaleValues( offset, scale )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        XMLStringList = '%s<%s>' % ( indent, self.moniker )
        XMLStringList += ''.join( self.values.toXMLList( '', **kwargs ) )
        XMLStringList += '</%s>' % self.moniker

        return( [ XMLStringList ] )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        data1 = cls( valuesModule.values.parseXMLNode( element[0], xPath, linkData ) )

        xPath.pop( )
        return( data1 )

class xs( data ) :

    moniker = 'xs'

class pdf( data ) :

    moniker = 'pdf'

class cdf( data ) :

    moniker = 'cdf'

class xs_pdf_cdf1d( xDataBaseModule.xDataFunctional ) :

    moniker = 'xs_pdf_cdf1d'
    dimension = 1

    def __init__( self, _xs, _pdf, _cdf, outerDomainValue = None, axes = None, interpolation = standardsModule.interpolation.linlinToken ) :

        xDataBaseModule.xDataFunctional.__init__( self, self.moniker, axes = axes, outerDomainValue = outerDomainValue )

        self.interpolation = interpolation

        if( not( isinstance( _xs, xs ) ) ) : raise TypeError( 'Invalid xs: "%s"' % type( _xs ) )
        self.__xs = _xs
        self.__xs.setAncestor( self )

        if( not( isinstance( _pdf, pdf ) ) ) : raise TypeError( 'Invalid pdf: "%s"' % type( _pdf ) )
        self.__pdf = _pdf
        self.__pdf.setAncestor( self )

        if( not( isinstance( _cdf, cdf ) ) ) : raise TypeError( 'Invalid cdf: "%s"' % type( _cdf ) )
        self.__cdf = _cdf
        self.__cdf.setAncestor( self )

        if( len( _xs ) != len( _pdf ) != len( _cdf ) ) :
            raise ValueError( 'lenghts not the same: %s %s %s' % ( len( _xs ), len( _pdf ), len( _cdf ) ) )

    def __len__( self ) :

        return( len( self.__xs ) )

    @property
    def xs( self ) :

        return( self.__xs )

    @property
    def pdf( self ) :

        return( self.__pdf )

    @property
    def cdf( self ) :

        return( self.__cdf )

    def convertUnits( self, unitMap ) :
        """See documentation for reactionSuite.convertUnits."""

        axes = self.axes
        if( axes is None ) : axes = self.ancestor.axes.copy( )      # BRB FIXME, the prior line should have worked.
        factors = axes.convertUnits( unitMap )
        if( factors[0] != 1.0 ) : self.pdf.offsetScaleData( 0.0, factors[0] )
        if( factors[1] != 1.0 ) : self.xs.offsetScaleData( 0.0, factors[1] )
        self.fixValuePerUnitChange( factors )

    def copy( self ) :

        return( self.__class__( self.xs.copy( ), self.pdf.copy( ), self.cdf.copy( ), 
                outerDomainValue = self.outerDomainValue, axes = self.axes, interpolation = self.interpolation ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        _pdf, _cdf = self.to_pdf_and_cdf( )
        return( _pdf.toPointwise_withLinearXYs( **kwargs ) )

    def to_pdf_and_cdf( self ) :
        """
        Returns two XYs1d instances. One for the pdf and one for the cdf. The interpolation is the same as self's interpolation.
        Ergo, the interpolation for the cdf is not correct as, in general, the correct cdf interpolation is not defined in GNDS.
        The units for the cdf will be correct but not the label.
        """

        _pdf = []
        _cdf = []
        for index, x1 in enumerate( self.__xs.values ) :
            _pdf.append( [ x1, self.__pdf.values[index] ] )
            _cdf.append( [ x1, self.__cdf.values[index] ] )

        linear = self.toLinearXYsClass( )

        axes = self.axes
        axes = self.ancestor.axes
        _pdf = linear( data = _pdf, axes = axes, interpolation = self.interpolation, outerDomainValue = self.outerDomainValue )
        if( axes is not None ) :
            axes = axes.copy( )
            axes[0].unit = axes[0].multiplyUnit( axes[1] )
        _cdf = linear( data = _cdf, axes = axes, interpolation = self.interpolation, outerDomainValue = self.outerDomainValue )

        return( _pdf, _cdf )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        xs_pdf_cdf1d_singleLine = kwargs.get( 'xs_pdf_cdf1d_singleLine', False )
        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent
        if( xs_pdf_cdf1d_singleLine ) : indent2 = ''

        attributeStr = xDataBaseModule.xDataFunctional.attributesToXMLAttributeStr( self )
        interpolationStr = ''
        if( self.interpolation != standardsModule.interpolation.linlinToken ) :
            interpolationStr = ' interpolation="%s"' % self.interpolation
        XMLStringList = [ '%s<%s%s%s>' % ( indent, self.moniker, attributeStr, interpolationStr ) ]
        if( self.isPrimaryXData( ) ) :
            if( self.axes is not None ) : XMLStringList += self.axes.toXMLList( indent = indent2, **kwargs )
        XMLStringList += self.xs.toXMLList( indent2, **kwargs )
        XMLStringList += self.pdf.toXMLList( indent2, **kwargs )
        XMLStringList += self.cdf.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        if( xs_pdf_cdf1d_singleLine ) : XMLStringList = [ ''.join( XMLStringList ) ]
        return( XMLStringList )

    def toLinearXYsClass( self ) :

        return( XYsModule.XYs1d )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData, axes = None ) :

        xPath.append( element.tag )

        values = { 'xs' : xs, 'pdf' : pdf, 'cdf' : cdf, 'axes' : axesModule.axes }
        data1 = {}
        outerDomainValue = element.get( 'outerDomainValue', None )
        interpolation = element.get( 'interpolation', standardsModule.interpolation.linlinToken )
        for child in element :
            data1[child.tag] = values[child.tag].parseXMLNode( child, xPath, linkData )

        axes = data1.pop( 'axes', None )
        data1 = cls( data1['xs'], data1['pdf'], data1['cdf'], outerDomainValue = outerDomainValue, interpolation = interpolation, axes = axes )

        xPath.pop( )
        return( data1 )

    @classmethod
    def parseXMLString( cls, XMLString ) :

        from xml.etree import cElementTree

        return( cls.parseXMLNode( cElementTree.fromstring( XMLString ), xPath=[], linkData={} ) )

    @classmethod
    def fromXYs( cls, xys, outerDomainValue = None ) :

        if( xys.interpolation not in [ standardsModule.interpolation.linlinToken, standardsModule.interpolation.flatToken ] ) :
            xys = xys.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 0, upperEps = 1e-8 )

        try :
            norm = xys.normalize( )
        except :
            norm = xys.copy( )
            for i1 in range( len( norm ) ) : norm[i1] = [ norm[i1][0], 1.0 ]
            if( len( norm ) > 2 ) : norm[ 0] = [ norm[ 0][0], 0.0 ]
            norm[-1] = [ norm[-1][0], 0.0 ]
            norm = norm.normalize( )

        _cdf = norm.runningIntegral( )
        _cdf[-1] = 1.0
        _xs, _pdf = norm.copyDataToXsAndYs( )
        _xs = xs( valuesModule.values( _xs ) )
        _pdf = pdf( valuesModule.values( _pdf ) )
        _cdf = cdf( valuesModule.values( _cdf ) )

        return( cls( _xs, _pdf, _cdf, outerDomainValue = outerDomainValue, interpolation = xys.interpolation ) )

if( __name__ == '__main__' ) :

    xys = XYsModule.XYs1d( [ [ 1, 0 ], [ 2, 1 ], [ 4, 1 ], [ 6, .5 ], [ 7, .1 ] ] )

    xs1 = xs_pdf_cdf1d.fromXYs( xys, 1e-5 )
    xsXML1 = xs1.toXML( )
    print( xsXML1 )

    xs2 = xs_pdf_cdf1d.parseXMLString( xsXML1 )
    xsXML2 = xs2.toXML( )

    print( xsXML1 == xsXML2 )
