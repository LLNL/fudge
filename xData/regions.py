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

__metaclass__ = type

import abc

from pqu import PQU as PQUModule

import XYs as XYsModule
import series1d as series1dModule
import multiD_XYs as multiD_XYsModule
import standards as standardsModule
import base as baseModule
import axes as axesModule
import uncertainties as uncertaintiesModule

domainEpsilon = 1e-15

class regions( baseModule.xDataFunctional ) :
    """Abstract base class for regions."""

    __metaclass__ = abc.ABCMeta

    def __init__( self, axes = None, 
            index = None, valueType = standardsModule.types.float64Token, value = None, label = None ) :

        baseModule.xDataFunctional.__init__( self, self.moniker, self.dimension, axes, index = index, valueType = valueType,
                value = value, label = label )
        self.regions = []

    def __len__( self ) :
        """Returns the number of regions in self."""

        return( len( self.regions ) )

    def __getitem__( self, index ) :
        """Returns the (i-1)^th region of self."""

        return( self.regions[index] )

    def __setitem__( self, index, region ) :
        """
        Set region (the right-hand-side) to the :math:`{i-1}^{th}` region of self. Region must be an instance of 
        base.xDataFunctional with the same dimension as self.
        If :math:`i > 0`, the following must also be met:
        
            - all the prior regions must already exists,
            - region's minimum domain value must be equal to the prior region's maximum domain value.
        """

# BRB need to check axes.
        if( not( isinstance( region, self.allowedSubElements( ) ) ) ) :
            raise TypeError( 'Invalid class for insertion: %s' % region.__class__ )
        n1 = len( self )
        if( not( 0 <= index <= n1 ) ) : raise IndexError( 'Index = %s not in range 0 <= index <= %d' % ( index, n1 ) )
        region = region.copy( )
        region.setAncestor( self )
        if( len( self ) == 0 ) :
            self.regions.append( region )
        else :
            if( index > 0 ) :
                if( self.regions[index-1].domainMax( ) != region.domainMin( ) ) : raise Exception( "Prior region's domainMax %s != new region's domainMin = %s" \
                    % ( self.regions[index-1].domainMax( ), region.domainMin( ) ) )
            if( index < ( n1 - 1 ) ) :
                if( self.regions[index+1].domainMin( ) != region.domainMax( ) ) : raise Exception(  "Next region's domainMin %s != new region's domainMax = %s" \
                    % ( self.regions[index-1].domainMin( ), region.domainMax( ) ) )
            if( index == n1 ) :
                self.regions.append( region )
            else :
                self.regions[index] = region

    def __add__( self, other ) :

        self2, other2 = self.copyToCommonRegions( other )
        for i1, region in enumerate( self2 ) : self2[i1] = region + other2[i1]
        return( self2 )

    __radd__ = __add__

    def __sub__( self, other ) :

        self2, other2 = self.copyToCommonRegions( other )
        for i1, region in enumerate( self2 ) : self2[i1] = region - other2[i1]
        return( self2 )

    def append( self, region ) :

            self[len( self )] = region

    def prepend( self, region ) :
        """
        Adds region to the beginning of regions.
        """

        if( not( isinstance( region, self.allowedSubElements( ) ) ) ) :
            raise TypeError( 'Invalid class for insertion: %s' % region.__class__ )

        if( len( self ) > 0 ) :
            if( region.domainMax( ) != self[0].domainMin( ) ) :
                raise Exception( "Prepending region's domainMax %s != first region's domainMin = %s" \
                    % ( region.domainMax( ), self[0].domainMin( ) ) )
        region = region.copy( )
        region.setAncestor( self )
        self.regions.insert( 0, region )

    def copy( self, index = None, value = None, axes = None ) :
        # FIXME some of this should probably move to 'returnAsClass' method

        if( index is None ) : index = self.index
        if( value is None ) : value = self.value
        if( axes is None ) : axes = self.axes
        newRegions = self.__class__( axes = axes, index = index,
                    valueType = self.valueType, value = value, label = self.label )
        for child in self : newRegions.append( child.copy( ) )
        return( newRegions )

    def splitInTwo( self, domainValue, epsilon = domainEpsilon ) :
        """
        Splits the region containing domainValue into two regions.
        """

        for i1, region in enumerate( self ) :
            domainMin, domainMax = region.domain( )
            if( domainMin < domainValue < domainMax ) :
                r1, r2 = region.splitInTwo( domainValue, epsilon = domainEpsilon )
                self.regions[i1] = r2
                self.regions.insert( i1, r1 )
                return

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.domainUnit( ) ).getValueAs( unitTo ) )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( self.regions[0].domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( self.regions[-1].domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnit( self ) :

        return( self.getAxisUnitSafely( self.dimension ) )

    def getValue( self, unitTo ) :

        if( self.value is None ) : return( self.value )
        return( self.domainUnitConversionFactor( unitTo ) * self.value )

    def rangeMin( self, unitTo = None, asPQU = False ) :

        return( min( [ region.rangeMin( unitTo = unitTo, asPQU = asPQU ) for region in self ] ) )

    def rangeMax( self, unitTo = None, asPQU = False ) :

        return( max( [ region.rangeMax( unitTo = unitTo, asPQU = asPQU ) for region in self ] ) )

    def integrate( self, **limits ):
        """
        Integrate a piecewise function. Supports limits for each axis.
        Example:
        >regions.integrate( energy_in = ('1e-5 eV', '10 eV'), energy_out = ('1 keV', '10 keV') )

        :param limits: dictionary containing limits for each independent axis (keyed by axis label or index).
        If an independent axis is missing from the dictionary, integrate over the entire domain of that axis.

        :return: float or PQU
        """

        integral = 0
        for region in self:
            if( isinstance( region, ( XYsModule.XYs1d, series1dModule.series ) ) ) :
                Min, Max = None, None
                if self.axes[-1].label in limits:
                    Min, Max = limits.pop( self.axes[-1].label )
                elif self.axes[-1].index in limits:
                    Min, Max = limits.pop( self.axes[-1].index )

                integral += region.integrate( domainMin = Min, domainMax = Max )
            elif isinstance( region, multiD_XYsModule.XYsnd ):
                integral += region.integrate( **limits )
            else:
                raise TypeError( "Unsupported class for integration: %s" % type( region ) )
        return( integral )

    def toString( self ) :

        s1 = ''
        for i1, region in enumerate( self ) : s1 += 'region %s\n' % i1 + region.toString( )
        return( s1 )

    def toXMLList( self, indent = '', **kwargs ) :

        indent2 = indent + kwargs.get( 'incrementalIndent', '  ' )

        attributeStr = baseModule.xDataFunctional.attributesToXMLAttributeStr( self )
        XMLList = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
        if self.isPrimaryXData( ) :
            if( self.axes is not None ) : XMLList += self.axes.toXMLList( indent2, **kwargs )
        for region in self.regions : XMLList += region.toXMLList( indent2, **kwargs )
        if( self.uncertainties ) : XMLList += self.uncertainties.toXMLList( indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker
        return( XMLList )

    @classmethod
    def parseXMLNode( cls, element, xPath, linkData, axes = None, **kwargs ) :

        xPath.append( element.tag )

        allowedSubElements = cls.allowedSubElements( )

        if( element.find( 'axes' ) is not None ) :
            axes = axesModule.axes.parseXMLNode( element.find( 'axes' ), xPath, linkData )
        index = int( element.get( 'index' ) ) if 'index' in element.keys() else None
        value = float( element.get( 'value' ) ) if 'value' in element.keys() else None
        label = element.get( 'label', None )

        regions = cls( axes = axes, index = index, value = value, label = label )
        uncertainties = None

        for child in element :
            if( child.tag == 'axes' ) :
                continue
            elif( child.tag == 'uncertainties' ) :
                uncertainties = uncertaintiesModule.uncertainties.parseXMLNode( child, xPath, linkData )
                continue
            else :
                subElementClass = None
                for subElement in allowedSubElements :
                    if( subElement.moniker == child.tag ) :
                        subElementClass = subElement
                        break
                if( subElementClass is None ) : raise Exception( 'unknown sub-element "%s" in element "%s"' % ( child.tag, cls.moniker ) )
                regions.append( subElementClass.parseXMLNode( child, xPath, linkData, axes = axes, **kwargs ) )
# FIXME, Should set regions.uncertainties in method
        if( uncertainties is not None ) : regions.uncertainties = uncertainties

        xPath.pop( )
        return regions

class regions1d( regions ) :

    moniker = 'regions1d'
    dimension = 1

    def copyToCommonRegions( self, other, epsilon = domainEpsilon ) :
        """
        Returns two regions instances that are copies of self and other but with regions added as
        needed so that each has the same number of regions and the region boundaries align.
        """

        if( self.dimension != other.dimension ) : raise Exception( 'self.dimension = %s not equal to other.dimension = %s' % \
                ( self.dimension, other.dimension ) )

        self2 = self.copy( )
        other2 = other.copy( )
        if( isinstance( other2, XYsModule.XYs1d ) ) :
            temp = self.__class__( )
            temp.append( other2 )
            other2 = temp 
        elif( not( isinstance( other2, regions ) ) ) :
            raise NotImplementedError( 'object of instance "%s" not implemented' % other2.__class__ )

        if(   self2.domainMin( ) < other2.domainMin( ) ) :
            region1 = other2[0].copy( )
            region1.setData( [ [ self2.domainMin( ), 0 ], [ other2.domainMin( ), 0 ] ] )
            other2.prepend( region1 )
        elif( self2.domainMin( ) > other2.domainMin( ) ) :
            region1 = self2[0].copy( )
            region1.setData( [ [ other2.domainMin( ), 0 ], [ self2.domainMin( ), 0 ] ] )
            self2.prepend( region1 )

        if(   self2.domainMax( ) > other2.domainMax( ) ) :
            region1 = other2[-1].copy( )
            region1.setData( [ [ other2.domainMax( ), 0 ], [ self2.domainMax( ), 0 ] ] )
            other2.append( region1 )
        elif( self2.domainMax( ) < other2.domainMax( ) ) :
            region1 = self2[-1].copy( )
            region1.setData( [ [ self2.domainMax( ), 0 ], [ other2.domainMax( ), 0 ] ] )
            self2.append( region1 )

        boundaries = set( )
        for region in self2[1:] : boundaries.add( float( region.domainMin( ) ) )
        for region in other2[1:] : boundaries.add( float( region.domainMin( ) ) )
        boundaries = sorted( boundaries )

        count = 0
        priorBoundary = None
        boundariesToMove = []
        for boundary in boundaries :
            if( priorBoundary is not None ) :
                if( ( boundary - priorBoundary ) < epsilon * max( abs( boundary ), abs( priorBoundary ) ) ) :
                    boundariesToMove.append( [ boundary, priorBoundary ] )
                    boundary = priorBoundary
                    count += 1
                    if( count > 1 ) : raise Exception( 'more than one boundary within epsilon = %s of %s' % ( epsilon, priorBoundary ) )
                else :
                    count = 0
            priorBoundary = boundary

        for boundary, priorBoundary in boundariesToMove : boundaries.remove( boundary )

        def processBoundaries( regions_, boundaries, boundariesToMove ) :

            for region in regions_ :
                domainMin, domainMax = region.domain( )
                for boundary, priorBoundary in boundariesToMove :
                    if(   domainMin == boundary ) :
                        region.tweakDomain( domainMin = priorBoundary, epsilon = epsilon )
                    elif( domainMax == boundary ) :
                        region.tweakDomain( domainMax = priorBoundary, epsilon = epsilon )

            for boundary in boundaries :
                for region in regions_ :
                    if( region.domainMin( ) < boundary < region.domainMax( ) ) :
                        regions_.splitInTwo( boundary, epsilon = epsilon )
                        break

        processBoundaries( self2, boundaries, boundariesToMove )
        processBoundaries( other2, boundaries, boundariesToMove )
        return( self2, other2 )

    def normalize( self, insitu = False, dimension = 1 ) :
        """
        The dimension argument is currently ignored, but kept to be compatible with calling from XYsnd.normalize.
        """

        factor = 1.0 / self.integrate()
        if( insitu ): copy = self
        else: copy = self.copy()

        for region in copy:
            region *= factor
        return( copy )

    def toPointwise_withLinearXYs( self, accuracy, lowerEps, upperEps, removeOverAdjustedPoints = False, axes = None ) :
        """
        Converts the regions of self into a single ``XYs.XYs1d`` instance that has 'lin-lin' interpolation. At the
        boundary between two abutting regions, the x-values are the same, which is not allowed for an ``XYs.XYs1d`` instance.
        The arguments ``lowerEps``, ``upperEps`` are used to smear the x-values at a boundary as follows. Let :math:`(x_l, y_l)` and
        :math:`(x_u, y_u)` be the abutting points for two abutting regions. If :math:`y_l = y_u` then the point :math:`(x_u, y_u)` is removed.
        Otherwise, if( lowerEps > 0 ) the point :math:`(x_l, y_l)` is moved to :math:`x = x_l * ( 1 - lowerEps )` (or :math:`x = x_l * ( 1 + lowerEps )`
        if :math:`x_l < 0`) and the :math:`y` value is interpolated at :math:`x`. If :math:`x` is less than the x-value of the point below :math:`(x_l, y_l)`
        and ``removeOverAdjustedPoints`` is True then the point :math:`(x_l, y_l)` is removed; otherwise, a raise is executed. Similarly
        for upperEps and the point :math:`(x_u, y_u)`.
        """

        def getAdjustedX( x, eps ) :

            if( x == 0. ) :
                x_ = eps
            elif( x < 0 ) :
                x_ = x * ( 1. - eps )
            else :
                x_ = x * ( 1. + eps )
            return( x_ )

        if( len( self.regions ) == 0 ) : return( XYsModule.XYs1d( [], axes = self.axes ) )

        if( lowerEps < 0. ) : raise Exception( 'lowerEps = %s must >= 0.' % lowerEps )
        if( upperEps < 0. ) : raise Exception( 'upperEps = %s must >= 0.' % upperEps )
        if( ( lowerEps == 0. ) and ( upperEps == 0. ) ) : raise Exception( 'lowerEps and upperEps cannot both be 0.' )

        for region in self.regions : accuracy = max( accuracy, region.getAccuracy( ) )
        xys = []
        for iRegion, region in enumerate( self.regions ) :
            _region = region.changeInterpolation( standardsModule.interpolation.linlinToken, accuracy = accuracy, lowerEps = 2 * lowerEps, upperEps = 2 * upperEps )
            _region = _region.copyDataToXYs( )
            if( iRegion > 0 ) :
                x12, y12 = xys[-1]
                x21, y21 = _region[0]
                if( y12 == y21 ) :              # Remove first point of region as it is the same as the last point.
                    del _region[0]
                else :
                    if( lowerEps != 0. ) :
                        x11, y11 = xys[-2]
                        x = getAdjustedX( x12, -lowerEps )
                        if( x <= x11 ) :
                            if( removeOverAdjustedPoints ) :
                                del xys[-1]
                            else :
                                raise Exception( 'Adjustment at %s makes new x = %s >= prior x = %s; eps = %s' % ( x12, x, x11, lowerEps ) )
                        else :
                            xys[-1] = [ x, XYsModule.pointwiseXY_C.interpolatePoint( standardsModule.interpolation.linlinToken, x, x11, y11, x12, y12 ) ]
                    if( upperEps != 0. ) :
                        x22, y22 = _region[1]
                        x = getAdjustedX( x21, upperEps )
                        if( x >= x22 ) :
                            if( removeOverAdjustedPoints ) :
                                del _region[0]
                            else :
                                raise Exception( 'Adjustment at %s makes new x = %s >= next x = %s; eps = %s' % ( x21, x, x22, upperEps ) )
                        else :
                            _region[0] = [ x, XYsModule.pointwiseXY_C.interpolatePoint( standardsModule.interpolation.linlinToken, x, x21, y21, x22, y22 ) ]
            xys += _region
        pointwiseClass = self.toLinearXYsClass( )
        pointwise = pointwiseClass( data = xys, axes = self.axes, accuracy = accuracy, value = self.value ) # FIXME - need more work to insure all parameters are set properly.
        return( pointwise )

    def toLinearXYsClass( self ) :

        return( XYsModule.XYs1d )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYsModule.XYs1d, series1dModule.series ) )

class regions2d( regions ) :

    moniker = 'regions2d'
    dimension = 2

    def toLinearXYsClass( self ) :

        return( regions )

    @staticmethod
    def allowedSubElements( ) :

        return( ( multiD_XYsModule.XYs2d, ) )

class regions3d( regions ) :

    moniker = 'regions3d'
    dimension = 3

    def toLinearXYsClass( self ) :

        return( regions )

    @staticmethod
    def allowedSubElements( ) :

        return( ( multiD_XYsModule.XYs3d, ) )
