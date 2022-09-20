# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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

from LUPY import ancestry as ancestryModule
from pqu import PQU as PQUModule

class Base(ancestryModule.AncestryIO):

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__( self )

class Uncertainty(ancestryModule.AncestryIO):

    moniker = "uncertainty"

    def __init__( self, form ) :

        ancestryModule.AncestryIO.__init__( self )

        if( not( isinstance( form, Base ) ) ) : raise TypeError( "Invalid uncertainty form." )
        self.__form = form
        self.__form.setAncestor( self )

    @property
    def form( self ) :

        return( self.__form )

    def copy( self ) :

        return( self.__class__( self.__form.copy( ) ) )

    def parentConvertingUnits( self, factors ) :

        self.__form.parentConvertingUnits( factors )

    def toXML_strList(self, indent = '', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XML_strList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XML_strList += self.__form.toXML_strList(indent = indent2, **kwargs)
        XML_strList[-1] += "</%s>" % self.moniker

        return XML_strList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        from . import standard as standardModule

        xPath.append(node.tag)

        child = node[0]
        if child.tag == standardModule.Standard.moniker:
            form = standardModule.Standard.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        else:
            raise Exception( 'Invalid child node = "%s".' % child.tag )

        uncertainty1 = Uncertainty(form)

        xPath.pop( )
        return uncertainty1

class Quantity(ancestryModule.AncestryIO):

    absolute = 'absolute'
    relative = 'relative'
    percent = 'percent'
    relations = ( absolute, relative, percent )

    def __init__( self, value, relation = absolute ) :

        ancestryModule.AncestryIO.__init__( self )

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

    def toXML_strList( self, indent = '', **kwargs ) :

        relation = ''
        if( self.relation != self.absolute ) : relation = ' relation="%s"' % self.relation

        return [ '%s<%s value="%s"%s/>' % ( indent, self.moniker, self.value, relation ) ]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append(node.tag)

        attributes = ( 'value', 'relation' )
        for attributeName in node.attrib:
            if attributeName not in attributes: raise ValueError('Attribute = "%s" not allowed.' % attributeName)

        value = cls.toValueType(node.attrib['value'])
        instance = cls(value, node.get('relation', Quantity.absolute))

        xPath.pop()

        return instance

    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )

class Number( Quantity ) :
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

class Double( Number ) :
    """
    This class is used to represent a (uncertainty) quantity whose value must be a float.
    """

    moniker = 'double'
    __valueType = float

    def __init__( self, value, relation = Quantity.absolute ) :

        if( isinstance( value, ( int, fractions.Fraction ) ) ) : value = self.valueType( value )
        Number.__init__( self, value, relation )

    @property
    def valueType( self ) :

        return( self.__valueType )

    def parentConvertingUnits( self, factors ) :

        if( self.relation == self.absolute ) : self.value *= factors[0]

    @classmethod
    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )
