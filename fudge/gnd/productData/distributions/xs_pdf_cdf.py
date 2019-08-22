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

from xData import ancestry as ancestryModule
from xData import standards as standardsModule
from xData import values as valuesModule
from xData import base as xDataBaseModule

class data( ancestryModule.ancestry ) :

    def __init__( self, values ) :

        if( not( isinstance( values, valuesModule.values ) ) ) : raise TypeError( 'Not a values type' )
        ancestryModule.ancestry.__init__( self )

        self.values = values

    def __len__( self ) :

        return( len( self.values ) )

    def copy( self ) :

        return( self.__class__( self.values ) )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        indent2 = indent + incrementalIndent

        XMLStringList = '%s<%s>' % ( indent, self.moniker )
        XMLStringList += ''.join( self.values.toXMLList( '', **kwargs ) )
        XMLStringList += '</%s>' % self.moniker

        return( [ XMLStringList ] )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData, axes = None, **kwargs ) :

        xPath.append( element.tag )

        _data = cls( valuesModule.values.parseXMLNode( element[0], xPath, linkData ) )

        xPath.pop( )
        return( _data )

class xs( data ) :

    moniker = 'xs'

class pdf( data ) :

    moniker = 'pdf'

class cdf( data ) :

    moniker = 'cdf'

class xs_pdf_cdf1d( xDataBaseModule.xDataFunctional ) :

    moniker = 'xs_pdf_cdf1d'
    dimension = 1

    def __init__( self, _xs, _pdf, _cdf, value = None, axes = None, interpolation = standardsModule.interpolation.linlinToken ) :

        xDataBaseModule.xDataFunctional.__init__( self, self.moniker, axes = axes, value = value )

        self.interpolation = interpolation

        if( not( isinstance( _xs, xs ) ) ) : raise TypeError( 'Invalid xs: "%s"' % type( _xs ) )
        self.xs = _xs
        _xs.setAncestor( self )

        if( not( isinstance( _pdf, pdf ) ) ) : raise TypeError( 'Invalid pdf: "%s"' % type( _pdf ) )
        self.pdf = _pdf
        _pdf.setAncestor( self )

        if( not( isinstance( _cdf, cdf ) ) ) : raise TypeError( 'Invalid cdf: "%s"' % type( _cdf ) )
        self.cdf = _cdf
        _cdf.setAncestor( self )

        if( len( _xs ) != len( _pdf ) != len( _cdf ) ) :
            raise ValueError( 'lenghts not the same: %s %s %s' % ( len( _xs ), len( _pdf ), len( _cdf ) ) )

    def convertUnits( self, unitMap ) :
        "See documentation for reactionSuite.convertUnits."

        raise Exception( 'Need to implement' )

    def copy( self ) :

        return( self.__class__( self.xs.copy( ), self.pdf.copy( ), self.cdf.copy( ), 
                value = self.value, axes = self.axes, interpolation = self.interpolation ) )

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
        XMLStringList += self.xs.toXMLList( indent2, **kwargs )
        XMLStringList += self.pdf.toXMLList( indent2, **kwargs )
        XMLStringList += self.cdf.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        if( xs_pdf_cdf1d_singleLine ) : XMLStringList = [ ''.join( XMLStringList ) ]
        return( XMLStringList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData, axes = None, **kwargs ) :

        xPath.append( element.tag )

        values = { 'xs' : xs, 'pdf' : pdf, 'cdf' : cdf }
        data = {}
        value = element.get( 'value', None )
        for child in element :
            data[child.tag] = values[child.tag].parseXMLNode( child, xPath, linkData )
        _data = cls( data['xs'], data['pdf'], data['cdf'], value = value )

        xPath.pop( )
        return( _data )

    @classmethod
    def parseXMLString( cls, XMLString ) :

        from xml.etree import cElementTree

        value = None

        return( cls.parseXMLNode( cElementTree.fromstring( XMLString ), xPath=[], linkData={} ) )

    @classmethod
    def fromXYs( cls, xys, value = None ) :

        if( xys.interpolation not in [ standardsModule.interpolation.linlinToken, standardsModule.interpolation.flatToken ] ) :
            xys = xys.toPointwise_withLinearXYs( accuracy = 1e-3, lowerEps = 0, upperEps = 1e-8 )
        try :
            norm = xys.normalize( )
        except :
            norm = xys.copy( )
            norm[0] = [ norm[0][0], 1.0 ]
            norm = norm.normalize( )
            print '    WARNGING: xs_pdf_cdf1d.fromXYs; distribution with 0 norm.'
        _cdf = norm.runningIntegral( )
        _xs, _pdf = norm.copyDataToXsAndYs( )
        _xs = xs( valuesModule.values( _xs ) )
        _pdf = pdf( valuesModule.values( _pdf ) )
        _cdf = cdf( valuesModule.values( _cdf ) )

        return( cls( _xs, _pdf, _cdf, value = value, interpolation = xys.interpolation ) )

if( __name__ == '__main__' ) :

    from xData import XYs as XYsModule

    xys = XYsModule.XYs1d( [ [ 1, 0 ], [ 2, 1 ], [ 4, 1 ], [ 6, .5 ], [ 7, .1 ] ] )

    xs1 = xs_pdf_cdf1d.fromXYs( xys, 1e-5 )
    xsXML1 = xs1.toXML( )
    print xsXML1

    xs2 = xs_pdf_cdf1d.parseXMLString( xsXML1 )
    xsXML2 = xs2.toXML( )

    print xsXML1 == xsXML2
