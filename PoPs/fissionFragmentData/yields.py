# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import ancestry as ancestryModule
from xData import values as valuesModule
from xData import xDataArray as arrayModule

from .. import misc as miscModule
from .. import suite as suiteModule

class values( miscModule.classWithLabelKey ) :

    moniker = 'values'

    def __init__( self, _values ) :

        ancestryModule.ancestry.__init__( self )

        self.__values = tuple( float( value ) for value in _values )

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
    def parseXMLNodeAsClass(cls, element, xPath, linkData):

        xPath.append( element.tag )
        values_ = cls( map(float, element.text.split()) )
        xPath.pop()
        return values_

class covariance( ancestryModule.ancestry ) :

    moniker = 'covariance'

    def __init__( self, _matrix ) :

        ancestryModule.ancestry.__init__( self )

        self.__matrix = _matrix

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
    def parseXMLNodeAsClass( cls, element, xPath, linkData ):
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
    def parseXMLNodeAsClass( cls, element, xPath, linkData ):
        xPath.append(element.tag)
        child = element[0]
        if child.tag == covariance.moniker:
            form = covariance.parseXMLNodeAsClass(child, xPath, linkData)
        else:
            raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))
        uncert_ = cls(form)
        xPath.pop()
        return uncert_

class yields( ancestryModule.ancestry ) :

    moniker = 'yields'

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

        self.__values = None
        self.__uncertainty = None

    @property
    def values( self ) :

        return( self.__values )

    @values.setter
    def values( self, _values ) :

        if( not( isinstance( _values, values ) ) ) : raise TypeError( 'Invalid values instance.' )
        self.__values = _values
        self.__values.setAncestor( self )

    @property
    def uncertainty( self ) :

        return( self.__uncertainty )

    @uncertainty.setter
    def uncertainty( self, _uncertainty ) :

        if( not( isinstance( _uncertainty, uncertainty ) ) ) : raise TypeError( 'Invalid uncertainty instance.' )
        self.__uncertainty = _uncertainty
        self.__uncertainty.setAncestor( self )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        if( self.__values is not None ) : XMLStringList += self.__values.toXMLList( indent2, **kwargs )
        if( self.__uncertainty is not None ) : XMLStringList += self.__uncertainty.toXMLList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ):

        xPath.append(element.tag)
        for child in element:
            if child.tag == values.moniker:
                self.values = values.parseXMLNodeAsClass(child, xPath, linkData)
            elif child.tag == uncertainty.moniker:
                self.uncertainty = uncertainty.parseXMLNodeAsClass(child, xPath, linkData)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return (self)

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ):

        xPath.append(element.tag)
        self = cls()
        xPath.pop()
        self.parseXMLNode(element, xPath, linkData)

        return self
