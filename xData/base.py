# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the xData class that must be the base class for all xData classes.
"""

import abc

__metaclass__ = type

from pqu import PQU as PQUModule

from . import ancestry as ancestryModule
from . import standards as standardsModule

def getArguments( kwargs, arguments ) :

    for argument in arguments : arguments[argument] = kwargs.get( argument, arguments[argument] )
    return( arguments )

class xData( ancestryModule.ancestry ) :

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

class xDataCoreMembers( xData ) :

    def __init__( self, moniker, index = None, label = None ) :

        xData.__init__( self )

        if( index is not None ) : index = int( index )
        self.index = index
        self.label = label

    def attributesToXMLAttributeStr( self ) :

        attributeStr = ''
        if( self.index is not None ) : attributeStr += ' index="%s"' % self.index
        if( self.label is not None ) : attributeStr += ' label="%s"' % self.label
        return( attributeStr )

    @property
    def index( self ) :

        return( self.__index )

    @index.setter
    def index( self, value ) :

        if( value is not None ) : value = int( value )
        self.__index = value

    @property
    def label( self ) :

        return( self.__label )

    @label.setter
    def label( self, value ) :

        if( value is not None ) :
            if( not( isinstance( value, str ) ) ) : raise TypeError( 'label must be a string, got %s' % type(value) )
        self.__label = value

    @staticmethod
    def getArguments( kwargs, arguments ) :

        return( getArguments( kwargs, arguments ) )

class xDataFunctional( xDataCoreMembers ) :

    ancestryMembers = ( '', ) # 'axes', )       # BRB, ignore axes until it is improved. Probably should have each sub-function have its own axes copy.

    def __init__( self, moniker, axes, index = None, valueType = standardsModule.types.float64Token, outerDomainValue = None, label = None ) :

        xDataCoreMembers.__init__( self, moniker, index = index, label = label )

        if( valueType is None ) : valueType = standardsModule.types.float64Token
        if( not( isinstance( valueType, str ) ) ) : raise TypeError( 'valueType must be a string' )
        if( valueType != standardsModule.types.float64Token ) : raise TypeError( 'valueType = "%s" not supported' % valueType )
        self.valueType = valueType

        if( outerDomainValue is not None ) : outerDomainValue = float( outerDomainValue )
        self.__outerDomainValue = outerDomainValue

        self.axes = axes

        self.uncertainty = None

    @property
    @abc.abstractmethod
    def dimension( self ) :

        pass

    def attributesToXMLAttributeStr( self ) :

        attributeStr = xDataCoreMembers.attributesToXMLAttributeStr( self )
        if( self.outerDomainValue is not None ) : attributeStr += ' outerDomainValue="%s"' % PQUModule.floatToShortestString( self.outerDomainValue, 12 )
        return( attributeStr )

    @property
    def axes( self ) :
        """Returns self's axes."""

        return( self.__axes )

    @axes.setter
    def axes( self, _axes ) :

        from . import axes as axesModule

        if( _axes is not None ) : 
            if( not( isinstance( _axes, ( axesModule.axes, axesModule.referenceAxes ) ) ) ) :
                raise TypeError( 'axes in not an axes or referenceAxes instance' )
            if( len( _axes ) <= self.dimension ) : raise Exception( 'len( axes ) = %d != ( self.dimension + 1 ) = %d' % ( len( _axes ), ( self.dimension + 1 ) ) )
            _axes = _axes.copy( )
            _axes.setAncestor( self )
        self.__axes = _axes

    @property
    def outerDomainValue( self ) :

        return( self.__outerDomainValue )

    @outerDomainValue.setter
    def outerDomainValue( self, outerDomainValue ) :

        if( outerDomainValue is not None ) : outerDomainValue = float( outerDomainValue )
        self.__outerDomainValue = outerDomainValue

    def fixValuePerUnitChange( self, factors ) :

        if( self.__outerDomainValue is not None ) : self.__outerDomainValue *= factors[self.dimension+1]

    def getAxisIndexByIndexOrName( self, indexOrName ) :

        if( isinstance( indexOrName, int ) ) :
            return( indexOrName )
        elif( isinstance( indexOrName, str ) ) :
            for index, axis in enumerate( self.axes ) :
                if( axis.label == indexOrName ) : return( index )
        raise TypeError( 'argument must be an integer or a string' )

    def getAxisUnitSafely( self, index ) :

        if( self.axes is None ) : return( '' )
        return( self.axes[index].unit )

    def getPrimaryXData( self ) :

        ancestor = self
        while( not( ancestor.isPrimaryXData( ) ) ) : ancestor = ancestor.ancestor
        return( ancestor )

    def isPrimaryXData( self ) :
        """Returns False if self is contained in a higher dimension xDataFunctional and False otherwise."""

        ancestry = self.ancestor
        if( ancestry is None ) : return( True )
        return( not( isinstance( ancestry, xDataFunctional ) ) )

    @abc.abstractmethod
    def toXMLList( self ): return []

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

def isXDataFunctional( object ) :
    """Returns True if object is an instance of xData and False otherwise."""

    return( isinstance( object, xData ) )

def getDomainValue( value, unit, default ) :

    if( value is None ) : return( default )
    if( isinstance( value, ( str, PQUModule.PQU ) ) ) : return( PQUModule.PQU( value ).getValueAs( unit ) )
    return( value )

def getDomainLimits( self, domainMin, domainMax, unit ) :

    defaultMin, defaultMax = self.domainMin, self.domainMax
    return( getDomainValue( domainMin, unit, defaultMin ), getDomainValue( domainMax, unit, defaultMax ) )

def processUnits( unit1, unit2, operator ) :

    if( operator not in [ '*', '/' ] ) : raise ArithmeticError( 'unsupported unit operation "%s"' % operator )
    result = eval( 'PQUModule.PQU( 1, unit1 ) %s PQUModule.PQU( 1, unit2 )' % operator )
    return( result.getUnitSymbol( ) )

def getDomainValue2( domainValue ) :

    if( isinstance( domainValue, PQUModule.PQU ) ) : return( domainValue )
    if( isinstance( domainValue, str ) ) :
        try :
            return( float( domainValue ) )
        except :
            return( PQUModule.PQU( domainValue ) )
    return( float( domainValue ) )
