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

    __metaclass__ = abc.ABCMeta

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

        if( self.relation == self.absolte ) : self.__value *= factors[0]

    @classmethod
    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )
