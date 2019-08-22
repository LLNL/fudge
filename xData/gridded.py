# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
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
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

__metaclass__ = type

import XYs
import standards as standardsModule
import base as baseModule
import axes as axesModule
import values as valuesModule
import array as arrayModule

class gridded( baseModule.xDataFunctional ) :

    moniker = 'gridded'

    def __init__( self, axes, array, copyArray = True,
            index = None, valueType = standardsModule.types.float64Token, value = None, label = None ) :

        for axis in axes :
            if( axis.index == 0 ) :
                if( not( isinstance( axis, axesModule.axis ) ) ) : raise Exception( 'dependent axis must not have a grid' )
            else :
                if( not( isinstance( axis, axesModule.grid ) ) ) : raise Exception( 'independent axis must have a grid' )
        if( not( isinstance( array, arrayModule.arrayBase ) ) ) : raise TypeError( 'array not value type' )

        baseModule.xDataFunctional.__init__( self, self.moniker, len( array ), axes, index = index, valueType = valueType,
                value = value, label = label )

        if( copyArray ) : array = array.copy( )
        self.array = array
        self.array.setAncestor( self )

    def __len__( self ) :
        """Returns the number of regions in self."""

        return( len( self.array ) )

    def copy( self ):

        return gridded( self.axes.copy(), self.array, copyArray=True, index=self.index, valueType=self.valueType,
                label=self.label )

    __deepcopy__ = __copy__ = copy

    def getGrid( self, index ) :

        if( 0 < index <= len( self ) ) : return( self.axes[index].grid )
        raise IndexError( 'index = %d out of range [1,%s]' % len( self ) )

    def getGridUnit( self, index ) :

        return( self.axes[index].unit )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = baseModule.xDataFunctional.attributesToXMLAttributeStr( self )
        XMLList = [ '%s<%s dimension="%s"%s>' % ( indent, self.moniker, self.dimension, attributeStr ) ]
        XMLList += self.axes.toXMLList( indent2, **kwargs )
        XMLList += self.array.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @classmethod
    def parseXMLNode( cls, element, xPath = [], linkData = {}, **kwargs ) :

        xPath.append( element.tag )

        axes = None
        if element.find('axes') is not None:
            axes = axesModule.axes.parseXMLNode( element.find('axes'), xPath )
        array = arrayModule.arrayBase.parseXMLNode( element[1], xPath, linkData )
        index = int( element.get('index') ) if 'index' in element.keys() else None
        value = float( element.get('value') ) if 'value' in element.keys() else None
        Grid = gridded( axes, array, copyArray = False, index=index, value=value, label=element.get('label') )
        xPath.pop( )
        return Grid
