# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from xData import link as linkModule
from xData import xDataArray as arrayModule

from .. import misc as miscModule
from . import nuclides as nuclidesModule

class NuclidesLink(linkModule.Link):

    moniker = 'nuclides'

class Values(miscModule.ClassWithLabelKey):

    moniker = 'values'

    def __init__( self, _values ) :

        miscModule.ClassWithLabelKey.__init__(self, 'label')

        self.__values = tuple( float( value ) for value in _values )

    def __getitem__( self, index ) :

        return( self.__values[index] )

    def __iter__( self ) :

        n1 = len( self.__values )
        for i1 in range( n1 ) : yield self.__values[i1]

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList[-1] += ' '.join( [ '%s' % value for value in self.__values ] )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append( element.tag )
        values_ = cls( map(float, element.text.split()) )
        xPath.pop()
        return values_

class Covariance( ancestryModule.AncestryIO):

    moniker = 'covariance'

    def __init__( self, _matrix ) :

        ancestryModule.AncestryIO.__init__(self)

        self.__matrix = _matrix

    def matrix( self ) :

        return( self.__matrix )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        if( 'valuesPerLine' not in kwargs ) : kwargs['valuesPerLine'] = 1000
        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__matrix.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        xPath.append(element.tag)
        child = element[0]
        if child.tag == arrayModule.ArrayBase.moniker:
            _matrix = arrayModule.ArrayBase.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        else:
            raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))
        covar_ = cls(_matrix)
        xPath.pop()
        return covar_

class Uncertainty(ancestryModule.AncestryIO):

    moniker = 'uncertainty'

    def __init__( self, form ) :

        ancestryModule.AncestryIO.__init__(self)

        self.__form = form

    def form( self ) :

        return( self.__form )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__form.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        xPath.append(element.tag)
        child = element[0]
        if child.tag == Covariance.moniker:
            form = Covariance.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        else:
            raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))
        uncert_ = cls(form)
        xPath.pop()
        return uncert_

class Yields(ancestryModule.AncestryIO):

    moniker = 'yields'

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__(self)

        self.__nuclides = None
        self.__values = None
        self.__uncertainty = None

    @property
    def nuclides(self):

        return self.__nuclides

    @nuclides.setter
    def nuclides(self, value):

        if not isinstance(value, (NuclidesLink, nuclidesModule.Nuclides)):
            raise TypeError( 'Invalid nuclides instance.' )
        self.__nuclides = value
        self.__nuclides.setAncestor(self)

    @property
    def values( self ) :

        return( self.__values )

    @values.setter
    def values( self, _values ) :

        if( not( isinstance( _values, Values ) ) ) : raise TypeError( 'Invalid values instance.' )
        self.__values = _values
        self.__values.setAncestor( self )

    @property
    def uncertainty( self ) :

        return( self.__uncertainty )

    @uncertainty.setter
    def uncertainty( self, _uncertainty ) :

        if( not( isinstance( _uncertainty, Uncertainty ) ) ) : raise TypeError( 'Invalid Uncertainty instance.' )
        self.__uncertainty = _uncertainty
        self.__uncertainty.setAncestor( self )

    def toXML_strList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        if( self.__nuclides is not None ) : XMLStringList += self.__nuclides.toXML_strList( indent2, **kwargs )
        if( self.__values is not None ) : XMLStringList += self.__values.toXML_strList( indent2, **kwargs )
        if( self.__uncertainty is not None ) : XMLStringList += self.__uncertainty.toXML_strList( indent2, **kwargs )
        XMLStringList[-1] += '</%s>' % self.moniker

        return( XMLStringList )

    def parseNode(self, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        for child in element:
            if child.tag == NuclidesLink.moniker:
                if child.get('href') is not None:
                    self.nuclides = NuclidesLink.parseNodeUsingClass(child, xPath, linkData, **kwargs)
                else:
                    self.nuclides = nuclidesModule.Nuclides.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif child.tag == Values.moniker:
                self.values = Values.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            elif child.tag == Uncertainty.moniker:
                self.uncertainty = Uncertainty.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            else:
                raise TypeError("Unexpected child node '%s' in %s" % (child.tag, element.tag))

        xPath.pop()
        return (self)

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        self = cls()
        xPath.pop()
        self.parseNode(element, xPath, linkData, **kwargs)

        return self
