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
This module contains the xData class that must be the base class for all xData classes.
"""

__metaclass__ = type

import ancestry as ancestryModule
import standards as standardsModule
from pqu import PQU as PQUModule

class xData( ancestryModule.ancestry ) :

    def __init__( self, moniker ) :

        ancestryModule.ancestry.__init__( self )

class xDataCoreMembers( xData ) :

    def __init__( self, moniker, index = None, label = None ) :

        xData.__init__( self, moniker )

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

class xDataFunctional( xDataCoreMembers ) :

    def __init__( self, moniker, dimension, axes, index = None, valueType = standardsModule.types.float64Token, 
                value = None, label = None ) :

        import xData.axes as axesModule

        xDataCoreMembers.__init__( self, moniker, index = index, label = label )

        self.__dimension = int( dimension )

        if( axes is None ) : 
            self.__axes = None
        else :
            if( not( isinstance( axes, axesModule.axes ) ) ) : raise TypeError( 'axes in not an axes instance' )
            if( len( axes ) <= dimension ) : raise Exception( 'len( axes ) = %d != ( dimension + 1 ) = %d' % ( len( axes ), ( dimension + 1 ) ) )
            self.__axes = axes.copy( )
            self.__axes.setAncestor( self )

        if( valueType is None ) : valueType = standardsModule.types.float64Token
        if( not( isinstance( valueType, str ) ) ) : raise TypeError( 'valueType must be a string' )
        if( valueType != standardsModule.types.float64Token ) : raise TypeError( 'valueType = "%s" not supported' % valueType )
        self.valueType = valueType
        if( value is not None ) : value = float( value )
        self.__value = value

        self.uncertainties = []

    def addUncertainties( self, uncertainty ) :

        raise TypeError( 'uncertainty other than "None" currently not supported' )
        self.uncertainties.append( uncertainty )       # BRB need to check for valid uncertainties.

    def attributesToXMLAttributeStr( self ) :

        attributeStr = xDataCoreMembers.attributesToXMLAttributeStr( self )
        if( self.value is not None ) : attributeStr += ' value="%s"' % self.value
        return( attributeStr )

    @property
    def axes( self ) :
        """Returns self's axes."""

        if( self.__axes is not None ) : return self.__axes
        while( self.isPrimaryXData( ) ) : return( None )
        return( self.ancestor.axes )

    @property
    def dimension( self ) :
        """Returns self's xData dimension."""

        return( self.__dimension )

    @property
    def value( self ) :

        return( self.__value )

    @value.setter
    def value( self, value ) :

        if( value is not None ) : value = float( value )
        self.__value = value

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
        while( not( ancestor.isPrimaryXData( ) ) ) : ancestor = ancestor.getAncestor( )
        return( ancestor )

    def isPrimaryXData( self ) :
        """Returns False if self is contained in a higher dimension xDataFunctional and False otherwise."""

        ancestry = self.getAncestor( )
        if( ancestry is None ) : return( True )
        return( not( isinstance( ancestry, xDataFunctional ) ) )

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

    defaultMin, defaultMax = self.domain( )
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
