# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import ancestry as ancestryModule
from xData import xDataArray as arrayModule

class values( ancestryModule.ancestry ) :

    moniker = 'values'

    def __init__( self, _values ) :

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( _values, ( tuple, list ) ) ) ) : raise TypeError( 'values must be a values instance.' )
        self.__values = tuple( value for value in _values )

    def __getitem__( self, index ) :

        return( self.__values[index] )

    def __iter__( self ) :

        n1 = len( self.__values )
        for i1 in range( n1 ) : yield self.__values[i1]

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList[-1] += ' '.join( [ '%s' % value for value in self.__values ] )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseXMLNode(cls, element, xPath, linkData):

        xPath.append(element.tag)
        values_ = cls( list( map( float, element.text.split() ) ) )
        xPath.pop()
        return values_

class covariance( ancestryModule.ancestry ) :

    moniker = 'covariance'

    def __init__( self, matrix ) :

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( matrix, arrayModule.diagonal ) ) ) : raise TypeError( 'Invalid matrix instance.' )
        self.__matrix = matrix

    def matrix( self ) :

        return( self.__matrix )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        if( 'valuesPerLine' not in kwargs ) : kwargs['valuesPerLine'] = 1000
        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__matrix.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):
        xPath.append(element.tag)
        child = element[0]
        if child.tag == arrayModule.arrayBase.moniker:
            _matrix = arrayModule.arrayBase.parseXMLNode(child, xPath, linkData)
        else:
            raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))
        covar_ = cls(_matrix)
        xPath.pop()
        return covar_

class uncertainty( ancestryModule.ancestry ) :

    moniker = 'uncertainty'

    def __init__( self, form ) :

        ancestryModule.ancestry.__init__( self )

        self.__form = form

    def form( self ) :

        return( self.__form )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__form.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData ):
        xPath.append(element.tag)
        child = element[0]
        if child.tag == covariance.moniker:
            form = covariance.parseXMLNode(child, xPath, linkData)
        else:
            raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))
        uncert_ = cls(form)
        xPath.pop()
        return uncert_
