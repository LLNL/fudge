# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the base uncertainty class.

This module contains the following classes:

    +---------------------------+-----------------------------------------------------------------------------------+
    | Class                     | Description                                                                       |
    +===========================+===================================================================================+
    | Base                      | This is the base class for all instances that can be a the for of uncertainty for |
    |                           | an instance of :py:class:`Uncertainty`.                                           |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Uncertainty               | This class presents an uncertain for a physical quantity.                         |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Quantity                  | This is the class for the GNDS computerCode version attribute.                    |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Number                    | This is the class for the GNDS computerCode version attribute.                    |
    +---------------------------+-----------------------------------------------------------------------------------+
    | Double                    | This is the class for the GNDS computerCode version attribute.                    |
    +---------------------------+-----------------------------------------------------------------------------------+
"""

import sys
import abc
import fractions

from LUPY import ancestry as ancestryModule
from pqu import PQU as PQUModule

class Base(ancestryModule.AncestryIO):
    """
    This is the base class for all instances that can be a the for of uncertainty for an instance of :py:class:`Uncertainty`.
    """

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__( self )

class Uncertainty(ancestryModule.AncestryIO):
    """
    This class presents an uncertainty for a physical quantity.

    The following table list the primary members of this class:

    +-----------------------+---------------------------------------------------------------------------+
    | Member                | Description                                                               |
    +=======================+===========================================================================+
    | form                  | This is the uncertainly function.                                         |
    +-----------------------+---------------------------------------------------------------------------+
    """

    moniker = "uncertainty"

    def __init__( self, form ) :
        """
        :param form:        This is the uncertainly function for *this* uncertainty instance..
        """

        ancestryModule.AncestryIO.__init__( self )

        if( not( isinstance( form, Base ) ) ) : raise TypeError( "Invalid uncertainty form." )
        self.__form = form
        self.__form.setAncestor( self )

    @property
    def form( self ) :
        """
        This method returns the form member of *self*.

        :returns:           An instance of :py:class:`Base`.
        """

        return( self.__form )

    def copy( self ) :
        """
        This method returns a copy of *self*.

        :returns:           An instance of :py:class:`Uncertainty`.
        """

        return( self.__class__( self.__form.copy( ) ) )

    def parentConvertingUnits( self, factors ) :
        """
        This method is call by the parent with the conversion factors used to change units.

        :param factors:     A list of python floats.
        """

        self.__form.parentConvertingUnits( factors )

    def toXML_strList(self, indent = '', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XML_strList = [ '%s<%s>' % ( indent, self.moniker ) ]
        XML_strList += self.__form.toXML_strList(indent = indent2, **kwargs)
        XML_strList[-1] += "</%s>" % self.moniker

        return XML_strList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :returns:           An instance of *cls* representing *node*.
        """

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
    """
    This class is the base class for all values (currently only floats) that an uncdertainty may have.

    The following table list the primary members of this class:

    +-----------------------+---------------------------------------------------------------------------+
    | Member                | Description                                                               |
    +=======================+===========================================================================+
    | value                 | This is the value of the uncertainty.                                     |
    +-----------------------+---------------------------------------------------------------------------+
    | relative              | This member indicates whether *value* is **absolute**, **relative** or    |
    |                       | percent.                                                                  |
    +-----------------------+---------------------------------------------------------------------------+
    """

    absolute = 'absolute'
    relative = 'relative'
    percent = 'percent'
    relations = ( absolute, relative, percent )

    def __init__( self, value, relation = absolute ) :
        """
        :param value:           This is the value of the uncertainty.
        :param relation:        This member indicates whether *value* is **absolute**, **relative** or percent.
        """

        ancestryModule.AncestryIO.__init__( self )

        self.value = value

        if( relation not in self.relations ) : raise TypeError( 'invalid relation = "%s"' % relation )
        self.__relation = relation

    @property
    def relation( self ) :
        """
        Thie method returns the *relation* member of *self*.

        :returns:       A python str.
        """

        return( self.__relation )

    @property
    def value( self ) :
        """
        Thie method returns the *value* member of *self*.

        :returns:       A type of self.valueType.
        """

        return( self.__value )

    @value.setter
    def value( self, value ) :
        """
        This member sets the *value* member of *self* to *value.

        :param value:       Must of type self.valueType.
        """

        if( not( isinstance( value, self.valueType ) ) ) : raise TypeError( 'Invalid value type must be a "%s"' % self.valueType )
        self.__value = value

    @property
    def valueType( self ) :
        """
        This method returns the *valueType* member of *self*.
        """

        return( self.__valueType )

    def copy( self ) :
        """
        This method returns a copy of *self*.

        :returns:           An instance of the class of *self*.
        """

        return( self.__class__( self.value, self.relation ) )

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        relation = ''
        if self.relation != self.absolute:
            relation = ' relation="%s"' % self.relation

        return ['%s<%s value="%s"%s/>' % (indent, self.moniker, self.valueToString(), relation)]

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :returns:           An instance of *cls* representing *node*.
        """

        xPath.append(node.tag)

        attributes = ( 'value', 'relation' )
        for attributeName in node.attrib:
            if attributeName not in attributes: raise ValueError('Attribute = "%s" not allowed.' % attributeName)

        value = cls.toValueType(node.attrib['value'])
        instance = cls(value, node.get('relation', Quantity.absolute))

        xPath.pop()

        return instance

    def toValueType( cls, value ) :
        """
        This method returns *value* as a type of the valueType of *self*.

        :param value:   Instance to convert.

        :returns:       An instance of the type of the valueType of *self*.
        """

        return( cls.__valueType( value ) )

class Number( Quantity ) :
    """
    This is an abstract base class for number quantities. This class adds the pqu and float methods.
    """

    def pqu( self, unit = None ) :
        """
        This method returns a :py:class:`PQUModule.PQU` instance of *self*'s value in units of *unit*. If *unit* is None, self's unit is used.

        :param unit:    The unit for the returned :py:class:`PQUModule.PQU` instance.

        :returns:       An instance of :py:class:`PQUModule.PQU`.
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
        This method returns a float instance of self's value in units of unit.

        :param unit:    The unit of the returned value.

        :returns:       A python float.
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
        """This seems redundant as it is defined in the base class Quantity."""

        return( self.__valueType )

    def parentConvertingUnits( self, factors ) :
        """
        This method is call by the parent with the conversion factors used to change units.

        :param factors:     A list of python floats.
        """

        if( self.relation == self.absolute ) : self.value *= factors[0]

    @classmethod
    def toValueType( cls, value ) :
        """
        This method returns *value* as a type of the valueType of this class.

        :param value:   Instance to convert.

        :returns:       An instance of the type of the valueType of this class.
        """

        return( cls.__valueType( value ) )

    def valueToString(self, precision=12):
        """
        This method returns a string version of the value of  *self*. See the function :py:func:`PQUModule.floatToShortestString`.

        :param precision:       The precision of the returned string.

        :returns:               A python str.
        """

        return PQUModule.floatToShortestString(self.value, min(max( 0, precision), 17), keepPeriod=True)
