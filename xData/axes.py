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

import string

import ancestry as ancestryModule
from xData.link import link
import values as valuesModule

noneGridToken = 'none'
pointsGridToken = 'points'
boundariesGridToken = 'boundaries'
parametersGridToken = 'parameters'
linkGridToken = 'link'

normalPDF = 'normal'
lognormalPDF = 'log-normal'

class axis( ancestryModule.ancestry ) :

    moniker = 'axis'

    def __init__( self, label, index, unit ) :
        """
        Returns a new instance of axis.
        """

        ancestryModule.ancestry.__init__( self )

        if( not( isinstance( label, str ) ) ) : raise Exception( 'label = "%s" is not a string' % label )
        self.__label = label.strip( )

        self.__index = int( index )

        self.setUnit( unit )

    def __str__( self ) :

        return( 'label="%s", index="%s", unit="%s"' % ( self.label, self.__index, self.unit ) )

    def __eq__( self, other ) :

# BRB: why is this method defined?
        return( isinstance( other, axis ) and ( self.label == other.label ) )

    def __ne__( self, other ) :

        return( not( self.__eq__( other ) ) )

    def copy( self, index = None ) :
        """Returns a new instance that is a copy of self."""

        if( index is None ) : index = self.__index
        return( axis( self.label, index, self.unit ) )

    __copy__ = copy
    __deepcopy__ = __copy__

    @property
    def index( self ) :

        return( self.__index )

    @property
    def label( self ) :

        return( self.__label )

    @property
    def unit( self ) :

        return( self.__unit )

    def plotLabel( self ) :

        label = self.label
        if( label == '' ) : label = 'unknown'
        if( self.unit != '' ) : label += ' (%s)' % self.unit
        return( label )

    def setUnit( self, unit ) :
        """Sets self's unit. Only checks that unit is a string. If unit is None, it is set to an empty string (i.e., '')."""

        if( unit is None ) : unit = ''
        if( not( isinstance( unit, str ) ) ) : raise TypeError( 'unit type "%s" is not a string' % type( unit ) )
        self.__unit = unit.strip( )

    def toXML( self, indent = '', **kwargs ) :

        XMLStr = '%s<%s index="%d" label="%s" unit="%s"/>' % ( indent, self.moniker, self.index, self.label, self.unit )
        return( XMLStr )

    def toXMLList( self, indent = '', **kwargs ) :

        return( [ self.toXML( indent = indent, **kwargs ) ] )

    @staticmethod
    def parseXMLNode( element, xPath ) :

        xPath.append( '%s[@index="%s"]' % ( axis.moniker, element.get( 'index' ) ) )

        _axis = axis( element.get( 'label' ), element.get( 'index' ), element.get( 'unit' ) )

        xPath.pop()
        return( _axis )

    def unitConversionFactor( self, newUnit ) :
        """Returns as a float the factor needed to convert self's unit to newUnit. If units are not compatible, a raise is executed."""

        from pqu import PQU
        return( PQU.PQU( 1., self.unit ).getValueAs( newUnit ) )

class grid( axis ) :

    moniker = 'grid'

    def __init__( self, label, index, unit, style, grid, uncertainty = None, pdf = normalPDF, interpolation = None ) :
        """
        Returns a new instance of grid.
        """

        axis.__init__( self, label, index, unit )

        if( style == linkGridToken ) :
            if( not isinstance( grid, link ) ):
                raise TypeError( "style = 'link' not consistent with grid type '%s'" % type( grid ) )
        else :
            if( style not in [ pointsGridToken, boundariesGridToken, parametersGridToken ] ) :
                raise ValueError( 'style = %s not supported' % style )
            if( not( isinstance( grid, valuesModule.values ) ) ) : TypeError( 'grid not values instance.' )

        self.__style = style
        self.grid = grid
        self.grid.setAncestor( self )

        self.interpolation = interpolation

            # BRB: uncertainty needs work.
        self.uncertainty = uncertainty
        if( not( isinstance( pdf, str ) ) ) : raise TypeError( 'pdf must be a string' )
        self.pdf = pdf

    @property
    def style( self ) :

        return( self.__style )

    def convertToUnit( self, unit ) :

        factor = self.unitConversionFactor( unit )
        self.setUnit( unit )
        self.grid = valuesModule.values( [ factor * value for value in self.grid ] )

    def copy( self, index = None ) :
        """Returns a new grid instance that is a copy of self."""

        if( index is None ) : index = self.__index
        return( grid( self.label, index, self.unit, self.style, self.grid, uncertainty = self.uncertainty, 
                pdf = self.pdf, interpolation = self.interpolation ) )

    __copy__ = copy
    __deepcopy__ = __copy__

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = ' style="%s"' % self.style
        if( self.interpolation is not None ) : attributeStr += ' interpolation="%s"' % self.interpolation
        if( self.uncertainty is not None ) :
            attributeStr += ' uncertainty="%s"' % self.uncertainty
            attributeStr += ' pdf="%s"' % self.pdf
        XMLStrList = [ '%s<%s index="%d" label="%s" unit="%s"%s>' % ( indent, self.moniker, self.index, self.label, self.unit, attributeStr ) ]
        XMLStrList += self.grid.toXMLList( indent2, **kwargs )
        XMLStrList[-1] += '</%s>' % self.moniker
        return( XMLStrList )

    @staticmethod
    def parseXMLNode( element, xPath ) :

        xPath.append( '%s[@index="%s"]' % ( grid.moniker, element.get( 'index' ) ) )

        style = element.get( 'style' )
        gridClass = {
            'link': link,
            'boundaries': valuesModule.values
        }.get( style )
        if( gridClass is None ) : raise Exception( "grid style '%s' not yet supported" % style )

        gridData = gridClass.parseXMLNode( element[0] )
        _grid = grid( element.get( 'label' ), element.get( 'index' ), element.get( 'unit' ), style, gridData )

        xPath.pop()
        return( _grid )

class axes( ancestryModule.ancestry ) :

    moniker = 'axes'

    def __init__( self, rank = None, labelsUnits = {} ) :
        """
        Constructor for ``axes`` class. For example::

            axes_ = axes( labelsUnits = { 0 : ( 'crossSection' , 'b' ), 1 : ( 'energy_in', 'eV' ) } )
        """

        ancestryModule.ancestry.__init__( self )

        if( rank is None ) :
            rank = 2
            if( len( labelsUnits ) > 0 ) : rank = len( labelsUnits )
        rank = int( rank )
        if( not( 0 < rank < 26 ) ) : raise Exception( 'rank = %d must be in the range [1, 25]' % rank )

        self.axes = []
        abcsOffset = string.ascii_lowercase.index( 'y' )
        for index in xrange( rank ) :
            label, unit = string.ascii_lowercase[abcsOffset-index], ''
            if( index in labelsUnits ) : label, unit = labelsUnits[index]
            self.axes.append( axis( label, index, unit ) )
            self.axes[-1].setAncestor( self, 'index' )

    def __eq__( self, other ) :

        if( isinstance( other, axes ) and ( len( self ) == len( other ) ) ) :
            for index, axis_ in enumerate( self.axes ) :
                if( axis_ != other[index] ) : return( False )
            return( True )
        return( False )

    def __ne__( self, other ) :

        return( not( self.__eq__( other ) ) )

    def __len__( self ) :

        return( len( self.axes ) )

    def __getitem__( self, index ) :

        return( self.axes[index] )

    def __setitem__( self, index, axisOrGrid ) :

        if( not( isinstance( axisOrGrid, ( axis, grid ) ) ) ) : raise TypeError( 'axisOrGrid is not an instance of axis or grid' )
        rank = len( self )
        index = int( index )
        if( index < 0 ) : index += rank 
        if( not( 0 <= index < rank ) ) : raise IndexError( "index = %s out of range for self of rank %s" % ( index, rank ) )
        self.axes[index] = axisOrGrid.copy( index = index )
        self.axes[index].setAncestor( self, 'index' )

    def __str__( self ) :

        l = [ str( axis ) for axis in self ]
        return( '\n'.join( l ) )

    def checkRank( self, rank ) :

        if( len( self ) != rank ) : raise Exception( "self's rank = %s != %s" % ( len( self ), rank ) )

    def copy( self ) :

        newAxes = axes( rank = len( self ) )
        for index, axis in enumerate( self ) : newAxes[index] = axis        # __setitem__ makes a copy, so no need to here.
        return( newAxes )

    __copy__ = copy
    __deepcopy__ = __copy__

    def toXML( self, indent = '', **kwargs ) :

        return( '\n'.join( self.toXMLList( indent = indent, **kwargs ) ) )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        XMLList = [ '%s<%s>' % ( indent, self.moniker ) ]
        xmlAxisStringList = []
        for axis in self : xmlAxisStringList.append( axis.toXML( indent = indent2, **kwargs ) )
        XMLList += reversed( xmlAxisStringList )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @staticmethod
    def parseXMLNode( axesElement, xPath ) :
        """Parse XML element with tag '<axes>'."""

        xPath.append( axesElement.tag )
        if( axesElement.tag == axes.moniker ) :
            axes_ = axes( rank = len( axesElement ) )
            for child in axesElement :
                childClass = {axis.moniker: axis, grid.moniker: grid}.get(child.tag)
                if childClass is None:
                    raise TypeError("Unexpected child element '%s' encountered in axes" % child.tag)
                index = child.get( "index" )
                axes_[index] = childClass.parseXMLNode( child, xPath )
        else :
            raise Exception( 'Invalid tag "%s" for axes' % ( axesElement.tag ) )
        xPath.pop()
        return( axes_ )

    @staticmethod
    def parseXMLString( axisString, xPath ) :

        from xml.etree import cElementTree
        return( axes.parseXMLNode( cElementTree.fromstring( axisString ), xPath = xPath ) )
