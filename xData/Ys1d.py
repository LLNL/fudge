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
This module contains the Ys1d class. 
"""

__metaclass__ = type

import standards as standardsModule
import base as baseModule
import axes as axesModule
from xData import values as valuesModule
import XYs as XYsModule
from pqu import PQU as PQUModule

class Ys1d( baseModule.xDataFunctional ) :

    moniker = 'Ys1d'
    dimension = 1

    ancestryMembers = ( 'Ys', )

    def __init__( self, Ys, interpolation = standardsModule.interpolation.linlinToken, axes = None,
            index = None, valueType = standardsModule.types.float64Token, value = None, label = None ) :

        baseModule.xDataFunctional.__init__( self, self.moniker, axes, index = index, valueType = valueType,
                value = value, label = label )

        if( not( isinstance( interpolation, str ) ) ) : raise TypeError( 'interpolation must be a string' )
        self.interpolation = interpolation

        if( not( isinstance( Ys, valuesModule.values ) ) ) : raise TypeError( 'Ys must be an instance of values.values.' )
        self.__Ys = Ys
        self.__Ys.setAncestor( self )

    def __len__( self ) :

        return( len( self.__Ys ) )

    def __getitem__( self, index ) :

        return( self.__Ys[index] )

    def __add__( self, other ) :

        if( len( self.__Ys ) == 0 ) : return( other.copy( ) )
        if( self.Ys.size != other.Ys.size ) : raise Exception( 'self.Ys.size = %d != other.Ys.size = %d' % ( self.Ys.size, other.Ys.size ) )

        if( self.__Ys.start <= other.Ys.start ) :
            ys1d_1 = self.copy( )
            ys1d_2 = other
        else :
            ys1d_1 = other.copy( )
            ys1d_2 = self 

        offset = ys1d_2.Ys.start - ys1d_1.Ys.start
        values = [ y for y in ys1d_1.Ys.values ]
        for i1, y2 in enumerate( ys1d_2.Ys ) : values[i1+offset] += y2
        ys1d_1.Ys.values = values

        return( ys1d_1 )

    @property
    def Ys( self ) :

        return( self.__Ys )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        factors = self.axes.convertUnits( unitMap )
        yFactor = factors[0]
        self.__Ys = valuesModule.values( [ yFactor * value for value in self.__Ys ] )
        self.fixValuePerUnitChange( factors )

    def copy( self ) :

        axes = self.axes
        if( axes is not None ) : axes = axes.copy( )
        Ys = self.__class__( self.__Ys.copy( ), interpolation = self.interpolation, axes = axes, 
                index = self.index, value = self.value, label = self.label )
        return( Ys )

    __copy__ = copy
    __deepcopy__ = __copy__

    def evaluate( self, domainValue, extrapolation = standardsModule.noExtrapolationToken, epsilon = 0 ) :

        pass

    @property
    def domainMin( self ) :

        return( self.axes[1].domainMin )

    @property
    def domainMax( self ) :

        return( self.axes[1].domainMax )

    @property
    def domainUnit( self ) :

        return( self.axes[1].domainUnit )

    def domainUnitConversionFactor( self, unitTo ) :

        return( self.axes[1].domainUnitConversionFactor )

    @property
    def domainGrid( self ) :

        return( self.axes[1].domainGrid )

    @property
    def rangeMin( self ) :

        return( min( self.__Ys.values ) )

    @property
    def rangeMax( self ) :

        return( max( self.__Ys.values ) )

    @property
    def rangeUnit( self ) :

        return( self.getAxisUnitSafely( 0 ) )

    def rangeUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.rangeUnit ).getValueAs( unitTo ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :

        cls = kwargs.pop( 'cls', XYsModule.XYs1d )

        xys = [ self.domainGrid[self.__Ys.start:self.__Ys.end], self.__Ys ]
        return( cls( data = xys, dataForm = 'xsandys', interpolation = self.interpolation, axes = self.axes, index = self.index,
                    valueType = self.valueType, value = self.value, label = self.label ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )
        outline = kwargs.get( 'outline', False )
        if( len( self ) < 6 ) : outline = False

        attributeStr = baseModule.xDataFunctional.attributesToXMLAttributeStr( self )
        if( self.interpolation != standardsModule.interpolation.linlinToken ) :
            attributeStr += ' interpolation="%s"' % self.interpolation

        XMLList = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ] 
        if( self.isPrimaryXData( ) ) :
            if( self.axes is not None ) : XMLList += self.axes.toXMLList( indent2 )
        XMLList += self.__Ys.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker

        return( XMLList )

    @classmethod
    def parseXMLNode( cls, xDataElement, xPath, linkData, axes = None ) :
        """
        Translates XML Ys1d into a Ys1d instance.
        """

        xPath.append( xDataElement.tag )

        attrs = {      'interpolation' : standardsModule.interpolation.linlinToken, 'label' : None, 'index' : None, 'value' : None  }
        attributes = { 'interpolation' : str,                                       'label' : str,  'index' : int,  'value' : float }
        for key, item in xDataElement.items( ) :
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key]( item )

        axes = None
        for subElement in xDataElement :
            if( subElement.tag == axesModule.axes.moniker ) :
                axes = axesModule.axes.parseXMLNode( subElement, xPath, linkData )
            elif( subElement.tag == valuesModule.values.moniker ) :
                Ys = valuesModule.values.parseXMLNode( subElement, xPath, linkData )
            else :
                raise TypeError( 'sub-element "%s" not valid' % subElement.tag )

        ys1d = cls( Ys, axes = axes, **attrs )

        xPath.pop( )
        return( ys1d )

    @staticmethod
    def defaultAxes( labelsUnits = None ) :
        """
        :param labelsUnits: dictionary of form { 0 : ( 'dependent label',   'dependent unit' ),
                                                 1 : ( 'independent label', 'independent unit' ) }
        :return: new axes instance
        """

        return( axesModule.axes( rank = 2, labelsUnits = labelsUnits ) )
