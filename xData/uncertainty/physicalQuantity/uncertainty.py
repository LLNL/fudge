# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the base uncertainty class.
"""

import sys
import abc
import fractions

from pqu import PQU as PQUModule

from ... import ancestry as ancestryModule

class base( ancestryModule.ancestry ) :

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

class uncertainty( ancestryModule.ancestry ) :

    moniker = "uncertainty"

    def __init__( self, form ) :

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( form, base ) ) ) : raise TypeError( "Invalid uncertainty form." )
        self.__form = form
        self.__form.setAncestor( self )

    @property
    def form( self ) :

        return( self.__form )

    def copy( self ) :

        return( self.__class__( self.__form.copy( ) ) )

    def parentConvertingUnits( self, factors ) :

        self.__form.parentConvertingUnits( factors )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs )  ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLStringList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XMLStringList += self.__form.toXMLList( indent = indent2, **kwargs )
        XMLStringList[-1] += "</%s>" % self.moniker

        return( XMLStringList )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        from . import standard as standardModule

        xPath.append( element.tag )

        child = element[0]
        if( child.tag == standardModule.standard.moniker ) :
            form = standardModule.standard.parseXMLNodeAsClass( child, xPath, linkData )
        else :
            raise Exception( 'Invalid child element = "%s".' % child.tag )

        _uncertainty = uncertainty( form )

        xPath.pop( )
        return( _uncertainty )

class quantity( ancestryModule.ancestry ) :

    absolute = 'absolute'
    relative = 'relative'
    percent = 'percent'
    relations = ( absolute, relative, percent )

    def __init__( self, value, relation = absolute ) :

        ancestryModule.ancestry.__init__( self )

        self.value = value

        if( relation not in self.relations ) : raise TypeError( 'invalid relation = "%s"' % relation )
        self.__relation = relation

    @property
    def relation( self ) :

        return( self.__relation )

    @property
    def value( self ) :

        return( self.__value )

    @value.setter
    def value( self, value ) :

        if( not( isinstance( value, self.valueType ) ) ) : raise TypeError( 'Invalid value type must be a "%s"' % self.valueType )
        self.__value = value

    @property
    def valueType( self ) :

        return( self.__valueType )

    def copy( self ) :

        return( self.__class__( self.value, self.relation ) )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs )  ) )

    def toXMLList( self, indent = '', **kwargs ) :

        relation = ''
        if( self.relation != self.absolute ) : relation = ' relation="%s"' % self.relation

        XMLStringList = [ '%s<%s value="%s"%s/>' % ( indent, self.moniker, self.value, relation ) ]

        return( XMLStringList )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( element.tag )

        attributes = ( 'value', 'relation' )
        for attributeName in element.attrib :
            if( attributeName not in attributes ) : raise ValueError( 'attribute = "%s" not allowed' % attributeName )

        value = cls.toValueType( element.attrib['value'] )
        self = cls( value, element.get( 'relation', quantity.absolute ) )

        xPath.pop()
        return( self )

    @classmethod
    def parseXMLStringAsClass( cls, string ) :

        from xml.etree import cElementTree
        return( cls.parseXMLNodeAsClass( cElementTree.fromstring( string ), [], [] ) )

    @classmethod
    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )

class number( quantity ) :
    """
    This is an abstract base class for number quantities. This class adds the pqu and float methods.
    """

    def pqu( self, unit = None ) :
        """
        Returns a PQU instance of self's value in units of unit. If unit is None, self's unit is used.
        """

        parent = self.ancestor.ancestor.ancestor
        if( self.relation == self.absolute ) :
            pqu = PQUModule.PQU( self.value, parent.unit )
            if( unit is not None ) : pqu.convertToUnit( unit )
        else :
            pqu = self.value * parent.pqu( unit )
            if( self.relation == self.percent ) : pqu /= 100
        return( pqu )

    def float( self, unit ) :
        """
        Returns a float instance of self's value in units of unit.
        """

        if( not( isinstance( unit, ( str, PQUModule.PhysicalUnit ) ) ) ) : raise TypeError( 'unit argument must be a str or a PQU.PhysicalUnit.' )
        return( float( self.pqu( unit ) ) )

class double( number ) :
    """
    This class is used to represent a (uncertainty) quantity whose value must be a float.
    """

    moniker = 'double'
    __valueType = float

    def __init__( self, value, relation = quantity.absolute ) :

        if( isinstance( value, ( int, fractions.Fraction ) ) ) : value = self.valueType( value )
        number.__init__( self, value, relation )

    @property
    def valueType( self ) :

        return( self.__valueType )

    def parentConvertingUnits( self, factors ) :

        if( self.relation == self.absolute ) : self.value *= factors[0]

    @classmethod
    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )
