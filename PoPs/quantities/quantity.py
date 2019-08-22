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
Contains the quantity class. This class is used to represent a (physical) quantity (e.g., mass, spin or halflife).
"""

import sys
import abc
import fractions

from pqu import PQU as PQUModule

from .. import misc as miscModule
from .. import suite as suiteModule

from xData.uncertainty.physicalQuantity import uncertainty as uncertaintyModule

__metaclass__ = type

class quantity( miscModule.classWithLabelKey ) :
    """
    This class is used to represent a (physical) quantity (e.g., mass, spin or halflife).
    A quantity has members label, value and can also have members unit, documentation, and uncertainty.
    """

    __metaclass__ = abc.ABCMeta

    moniker = 'quantity'
    __valueType = None
    baseUnit = None

    def __init__( self, label, value, unit, documentation = '' ) :

        miscModule.classWithLabelKey.__init__( self, label )

        self.value = value

        self.unit = unit

        if( not( isinstance( documentation, str ) ) ) : raise TypeError( 'documentation must be a string' )
        self.__documentation = documentation

        self.uncertainty = None

    def __str__( self ) :

        return( "%s %s" % ( self.value, self.unit ) )

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

    @property
    def unit( self ) :

        return( self.__unit )

    @unit.setter
    def unit( self, unit ) :

        if( isinstance( unit, str ) ) : unit = stringToPhysicalUnit( unit )
        if( not( isinstance( unit, PQUModule.PhysicalUnit ) ) ) : raise TypeError( 'unit must be a PQU.PhysicalUnit' )
        if( self.baseUnit is not None ) :
            if( not( self.baseUnit.isCompatible( unit ) ) ) : raise ValueError( 'unit "%s" not compatible with baseUnit "%s"' % ( unit, self.baseUnit ) )
        self.__unit = unit

    @property
    def uncertainty( self ) :

        return( self.__uncertainty )

    @uncertainty.setter
    def uncertainty( self, _uncertainty ) :

        if( _uncertainty is not None ) :
            if( not( isinstance( _uncertainty, uncertaintyModule.uncertainty ) ) ) : raise TypeError( 'Invalid uncertainty instance.' )

        self.__uncertainty = _uncertainty
        if( self.__uncertainty is not None ) : self.__uncertainty.setAncestor( self )

    @property
    def documentation( self ) :

        return( self.__documentation )

    def convertUnits( self, unitMap ) :

        unit, factor = PQUModule.convertUnits( self.unit, unitMap )
        if( abs( factor - 1 ) > ( 4 * sys.float_info.epsilon ) ) :
            NotImplementedError( 'Conversion of units for quantity of "%s" not implemented: from unit "" to unit "%s"' %
                    ( self.moniker, self.unit, unit ) )

    def copy( self ) :

        _quantity = self.__class__( self.label, self.value, self.unit, self.documentation )
        if( self.__uncertainty is not None ) : _quantity.uncertainty = self.__uncertainty.copy( )
        return( _quantity )

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs )  ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        indent3 = indent2 + kwargs.get( 'incrementalIndent', '  ' )

        attributes = ''
        if( not( self.unit.isDimensionless( ) ) ) : attributes += ' unit="%s"' % self.unit

        ending = '>'
        if( ( self.documentation == '' ) and ( self.uncertainty is None ) ) : ending = '/>'
        XMLStringList = [ '%s<%s label="%s" value="%s"%s%s' % 
                ( indent, self.moniker, self.label, self.value, attributes, ending ) ]

        if( ending == '>' ) :
            if( self.documentation != '' ) :
                XMLStringList.append( '%s<documentation>' % indent2 )
                XMLStringList.append( '%s%s</documentation>' % ( indent3, self.documentation ) )
            if( self.uncertainty is not None ) : XMLStringList += self.uncertainty.toXMLList( indent = indent2, **kwargs )
            XMLStringList[-1] += '</%s>' % self.moniker
        return( XMLStringList )

    def parseXMLNode( self, element, xPath, linkData ) :

        xPath.append( element.tag )

        documentation = ''
        for child in element :
            if( child.tag == 'uncertainty' ) :
                self.uncertainty = uncertaintyModule.uncertainty.parseXMLNodeAsClass( child, xPath, linkData )
#            elif( child.tag 'documentation' ) :
#                if( child.tag == 'documentation' ) : documentation = child.findtext( )
            else :
                raise ValueError( 'child element with tag "%s" not allowed' % child.tag )

        xPath.pop( )
        return( self )

    @classmethod
    def parseXMLNodeAsClass( cls, element, xPath, linkData ) :

        xPath.append( '%s[@label="%s"]' % ( element.tag, element.get( 'label' ) ) )

        attributes = ( 'label', 'value', 'unit' )
        for attributeName in element.attrib :
            if( attributeName not in attributes ) : raise ValueError( 'attribute = "%s" not allowed' % attributeName )

        value = cls.toValueType( element.attrib['value'] )
        unit = stringToPhysicalUnit( element.get( 'unit', '' ) )

        self = cls( element.attrib['label'], value, unit )
        xPath.pop()

        self.parseXMLNode( element, xPath, linkData )

        return( self )

    @classmethod
    def parseXMLStringAsClass( cls, string ) :

        from xml.etree import cElementTree
        return( cls.parseXMLNodeAsClass( cElementTree.fromstring( string ), [], [] ) )

    @classmethod
    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )

class string( quantity ) :
    """
    This is an abstract base class for string quantities.
    """

    moniker = 'string'
    __valueType = str

    @property
    def valueType( self ) :

        return( self.__valueType )

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

        pqu = PQUModule.PQU( self.value, self.unit )
        if( unit is not None ) : pqu.convertToUnit( unit )
        return( pqu )

    def float( self, unit ) :
        """
        Returns a float instance of self's value in units of unit.
        """

        if( not( isinstance( unit, ( str, PQUModule.PhysicalUnit ) ) ) ) : raise TypeError( 'unit argument must be a str or a PQU.PhysicalUnit.' )
        return( float( self.pqu( unit ) ) )

class integer( number ) :
    """
    This class is used to represent a (physical) quantity whose value must be an integer.
    """

    moniker = 'integer'
    __valueType = int

    @property
    def valueType( self ) :

        return( self.__valueType )

    @classmethod
    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )

class double( number ) :
    """
    This class is used to represent a (physical) quantity whose value must be a float.
    """

    moniker = 'double'
    __valueType = float

    def __init__( self, label, value, unit, documentation = '' ) :

        if( isinstance( value, ( int, fractions.Fraction ) ) ) : value = self.valueType( value )
        number.__init__( self, label, value, unit, documentation = documentation )

    @property
    def valueType( self ) :

        return( self.__valueType )

    def convertUnits( self, unitMap ) :

        unit, factor = PQUModule.convertUnits( self.unit, unitMap )
        unit = stringToPhysicalUnit( unit )
        self.value = self.value * factor
        self.unit = unit

    @classmethod
    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )

class fraction( number ) :
    """
    This class is used to represent a (physical) quantity whose value must be a rational number
    (e.g., a fraction like 1/2, 3, 5/6).
    """

    moniker = 'fraction'
    __valueType = fractions.Fraction

    def __init__( self, label, value, unit, documentation = '' ) :

        if( isinstance( value, ( int, str ) ) ) : value = self.valueType( value )
        number.__init__( self, label, value, unit, documentation = documentation )

    @property
    def valueType( self ) :

        return( self.__valueType )

    @classmethod
    def toValueType( cls, value ) :

        return( cls.__valueType( value ) )

class suite( suiteModule.suite ) :

    __metaclass__ = abc.ABCMeta

    _allowedClasses = []

    def __init__( self ) :

        suiteModule.suite.__init__( self, self._allowedClasses )

class numberSuite( suite ) :
    """
    This is an abstract base class for a number suite. This class adds the pqu and float methods.
    """

    __metaclass__ = abc.ABCMeta

    def pqu( self, unit = None ) :
        """
        Returns a PQU instance of self's recommended value in units of unit. If unit is None, self's unit is used.
        """

        return( self[0].pqu( unit = unit ) )

    def float( self, unit ) :
        """
        Returns a float instance of self's recommended value in units of unit.
        """

        return( self[0].float( unit ) )

def stringToPhysicalUnit( string ) :

    return( PQUModule._findUnit( string ) )
