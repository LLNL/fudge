# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from xData import enums as enumsModule
from . import values as valuesModule
from . import base as baseModule
from . import XYs1d as XYs1dModule

class Data(ancestryModule.AncestryIO):

    def __init__( self, values ) :

        if( not( isinstance( values, valuesModule.Values ) ) ) : raise TypeError( 'Not a values type' )
        ancestryModule.AncestryIO.__init__(self)

        self.values = values

    def __len__( self ) :

        return( len( self.values ) )

    def copy( self ) :

        return( self.__class__( self.values ) )

    def offsetScaleData( self, offset, scale ) :

        self.values.offsetScaleValues( offset, scale )

    def toXML_strList(self, indent = '', **kwargs):

        XMLStringList = '%s<%s>' % ( indent, self.moniker )
        XMLStringList += ''.join(self.values.toXML_strList('', **kwargs))
        XMLStringList += '</%s>' % self.moniker

        return [ XMLStringList ]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        data1 = cls(valuesModule.Values.parseNodeUsingClass(node[0], xPath, linkData, **kwargs))

        xPath.pop()

        return data1

class Xs( Data ) :

    moniker = 'xs'

class Pdf( Data ) :

    moniker = 'pdf'

class Cdf( Data ) :

    moniker = 'cdf'

class Xs_pdf_cdf1d( baseModule.XDataFunctional ) :

    moniker = 'xs_pdf_cdf1d'
    dimension = 1
    ancestryMembers = ( 'xs', 'pdf', 'cdf' )

    def __init__(self, _xs, _pdf, _cdf, axes=None, label=None, index=None, outerDomainValue=None, interpolation=enumsModule.Interpolation.linlin):

        baseModule.XDataFunctional.__init__(self, axes=axes, label=label, index=index, outerDomainValue=outerDomainValue)

        self.interpolation = enumsModule.Interpolation.checkEnumOrString(interpolation)

        if not isinstance(_xs, Xs): raise TypeError('Invalid Xs: "%s"' % type(_xs))
        self.__xs = _xs
        self.__xs.setAncestor(self)

        if not isinstance(_pdf, Pdf): raise TypeError('Invalid pdf: "%s"' % type(_pdf))
        self.__pdf = _pdf
        self.__pdf.setAncestor(self)

        if not isinstance(_cdf, Cdf): raise TypeError('Invalid cdf: "%s"' % type(_cdf))
        self.__cdf = _cdf
        self.__cdf.setAncestor(self)

        if len(_xs) != len(_pdf) != len(_cdf):
            raise ValueError('lenghts not the same: %s %s %s' % ( len(_xs), len(_pdf), len(_cdf) ))

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

    def toXML_strList(self, indent = '', **kwargs):

        xs_pdf_cdf1d_singleLine = kwargs.get('xs_pdf_cdf1d_singleLine', False)
        incrementalIndent = kwargs.get('incrementalIndent', '  ')
        indent2 = indent + incrementalIndent
        if xs_pdf_cdf1d_singleLine: indent2 = ''

        attributeStr = baseModule.XDataFunctional.attributesToXMLAttributeStr(self)
        interpolationStr = ''
        if self.interpolation != enumsModule.Interpolation.linlin:
            interpolationStr = ' interpolation="%s"' % self.interpolation
        XML_strList = [ '%s<%s%s%s>' % ( indent, self.moniker, attributeStr, interpolationStr ) ]
        if self.isPrimaryXData():
            if self.axes is not None: XML_strList += self.axes.toXML_strList(indent = indent2, **kwargs)
        XML_strList += self.xs.toXML_strList(indent2, **kwargs)
        XML_strList += self.pdf.toXML_strList(indent2, **kwargs)
        XML_strList += self.cdf.toXML_strList(indent2, **kwargs)
        XML_strList[-1] += '</%s>' % self.moniker

        if xs_pdf_cdf1d_singleLine: XML_strList = [ ''.join(XML_strList) ]

        return XML_strList

    def toLinearXYsClass( self ) :

        return( XYs1dModule.XYs1d )

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parses *node* into class *cls*.
        """

        attributes, extraAttributes = baseModule.XDataFunctional.parseBareNodeCommonAttributes(node, xPath, allowInterapolation=True) # parseBareNodeCommonAttributes adds to xPath.
        if len(extraAttributes) > 0: raise Exception('Invalid attributes: %s.' % ( ', '.join(list(extraAttributes.keys())) ))

        datas = {}
        children = { 'xs': Xs, 'pdf': Pdf, 'cdf': Cdf }
        for child in node:
            if child.tag in children: datas[child.tag] = children[child.tag].parseNodeUsingClass(child, xPath, linkData, **kwargs)

        axes = kwargs.get('axes', None)
        if axes is not None: attributes['axes'] = axes

        instance = cls(datas['xs'], datas['pdf'], datas['cdf'], **attributes)

        extraNodesList = baseModule.XDataFunctional.parseNodeStandardChildren(instance, node, xPath, linkData, **kwargs)
        extraNodes = {}
        for extraNode in extraNodesList: extraNodes[extraNode.tag] = extraNode
        extraNodes.pop('xs')
        extraNodes.pop('pdf')
        extraNodes.pop('cdf')

        if len(extraNodes) > 0: raise Exception('Invalid nodes: %s.' % (', '.join([extraNode.tag for extraNode in extraNodes])))

        xPath.pop()

        return instance

    @classmethod
    def fromXYs( cls, xys, outerDomainValue = None ) :

        if xys.interpolation not in [enumsModule.Interpolation.linlin, enumsModule.Interpolation.flat]:
            xys = xys.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 0, upperEps = 1e-8 )

        if outerDomainValue is None:
            outerDomainValue = xys.outerDomainValue

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
        _xs = Xs( valuesModule.Values( _xs ) )
        _pdf = Pdf( valuesModule.Values( _pdf ) )
        _cdf = Cdf( valuesModule.Values( _cdf ) )

        return( cls( _xs, _pdf, _cdf, outerDomainValue = outerDomainValue, interpolation = xys.interpolation ) )
