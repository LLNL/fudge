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

"""
    This module contains the class ``XYs``. This class treats a list of :math:`(x_i, y_i)` pairs as if they were the function :math:`y(x)`.
    That is, it is a numerical representation of :math:`f(x)`. The function :math:`y(x)` is also called a 1-dimensional or univariant
    function. As an example, let 
    
        .. math::
            y_1(x) = a + b * x

    and

        .. math::
            y_2(x) = c * \exp( - d * x )
            
    then, :math:`y_1` and :math:`y_2` can be added, subtracted, multiplied and divided, for example. Similiarly, two XYs instances can
    be added, subtracted, multiplied and divided. The only restriction on the :math:`(x_i, y_i)` pairs is that :math:`x_i < x_{i+1}`.

    The ``XYs`` class uses the class :py:class:`numericalFunctions.lib.pointwiseXY_C.pointwiseXY_C` as a base class.

    Members of an XYs instance of interest to most are:

        :accuracy:        see the base class :py:class:`numericalFunctions.lib.pointwiseXY_C.pointwiseXY_C`,
        :biSectionMax:    see the base class :py:class:`numericalFunctions.lib.pointwiseXY_C.pointwiseXY_C`,
        :infill:          see the base class :py:class:`numericalFunctions.lib.pointwiseXY_C.pointwiseXY_C`,
        :axes:            Description of the x and y data attributes (e.g., label, units).

    A ``XYs`` instance can be the ``XY`` instance of a 2-dimensional function (i.e., a ``multiD_XYs`` with
    dimension of 2). In this case, the ``XYs`` instances are secondary instances.
"""

"""
Notes:
    1) plot method still using fudge stuff.
    2) need list of values with unit. Used, for example, for method groupOneFunction.
"""

__metaclass__ = type

import base as baseModule
import axes as axesModule
import values as valuesModule
import standards as standardsModule
import uncertainties as uncertaintiesModule
from pqu import PQU

from numericalFunctions import pointwiseXY_C, pointwiseXY

defaultAccuracy = 1e-3
xAxisIndex = 1
yAxisIndex = 0
domainEpsilon = 1e-15

def return_pointwiseXY_AsXYs( self, C, units = {}, index = None, value = None, axes = None, template = None ) :

    if( index is None ) : index = self.index
    if( axes is None ) : axes = self.axes
    c = XYs( C, accuracy = 0.001, axes = axes, biSectionMax = 4, infill = 1, safeDivide = 0, 
            index = index, value = value )
    if( c.axes is not None ) :
        for k in units : c.axes[k].setUnit( units[k] )
    return( c )

def otherToSelfsUnits( self, other, checkXOnly = False ) :

    if( not( isinstance( other, XYs ) ) ) : raise TypeError( 'other instance not XYs instance' )
    if( ( self.axes is None ) and ( other.axes is None ) ) : return( other )
    yScale = 1
    xUnitSelf = PQU._getPhysicalUnit( self.axes[xAxisIndex].unit )
    xScale = xUnitSelf.conversionFactorTo( other.axes[xAxisIndex].unit )
    if( not( checkXOnly ) ) :
        yUnitSelf = PQU._getPhysicalUnit( self.axes[yAxisIndex].unit )
        yScale = yUnitSelf.conversionFactorTo( other.axes[yAxisIndex].unit )
    if( ( xScale != 1 ) or ( yScale != 1 ) ) : other = other.scaleOffsetXAndY( xScale = xScale, yScale = yScale )
    return( other )

def getValueAsUnit( unit, quantity ) :

    if( not( isinstance( quantity, PQU.PQU ) ) ) : raise TypeError( 'Quantity is not an instance of PQU.PQU' )
    return( quantity.getValueAs( unit ) )

def getValueAsAxis( axis, quantity ) :

    return( getValueAsUnit( axis.unit, quantity ) )

def getOtherAndUnit( self, other ) :
    """
    For internal use only. This function is used by the multiply and divide methods. Other must be a XYs instance or
    an object convertible to a PQU object.
    """

    if( isinstance( other, XYs ) ) :
        yUnit = other.getAxisUnitSafely( yAxisIndex )
        other = otherToSelfsUnits( self, other, checkXOnly = True )
    elif( isinstance( other, pointwiseXY ) ) :
        yUnit = ''
    else :                  # Must be an object convertible to a PQU instance.
        if( not( isinstance( other, PQU.PQU ) ) ) : other = PQU.PQU( other )
        yUnit = other.getUnitSymbol( )
        other = other.getValue( )
    return( other, yUnit )

def allow_XYsWithSameUnits_orNumberWithSameUnit( self, other ) :

    if( isinstance( other, XYs ) ) :
         other = otherToSelfsUnits( self, other )
    else :
        yUnit = self.getAxisUnitSafely( yAxisIndex )
        other = PQU.PQU( other, checkOrder = False ).getValueAs( yUnit )
    return( other )

class XYs( pointwiseXY, baseModule.xDataFunctional ) :

    moniker = 'XYs'
    mutableYUnit = True     # For __imul__ and __idiv__.

    def __init__( self, data = [], dataForm = "xys", interpolation = standardsModule.interpolation.linlinToken, axes = None,
            index = None, valueType = standardsModule.types.float64Token, value = None, label = None, 
            sep = ' ', accuracy = defaultAccuracy,
            initialSize = 10, overflowSize = 10, biSectionMax = 6, infill = True, safeDivide = False ) :
        """
        Constructor for XYs class. dataForm can be 'xys', 'xsandys' or 'list'.
        """

        baseModule.xDataFunctional.__init__( self, self.moniker, 1, axes, index = index, valueType = valueType,
                value = value, label = label )

        if( not( isinstance( interpolation, str ) ) ) : raise TypeError( 'interpolation must be a string' )
        self.interpolation = interpolation

        if( not( isinstance( sep, str ) ) ) : raise TypeError( 'sep must be of type str' )
        if( len( sep ) != 1 ) : raise TypeError( 'sep length must be 1 not %d' % len( sep ) )
        self.__sep = sep

        initialSize = max( initialSize, len( data ) )
        pointwiseXY.__init__( self, data = data, dataForm = dataForm, initialSize = initialSize, overflowSize = overflowSize,
            accuracy = accuracy, biSectionMax = biSectionMax, interpolation = interpolation, infill = infill, 
            safeDivide = safeDivide )

    def __abs__( self ) :

        return( self.returnAsClass( self, pointwiseXY.__abs__( self ) ) )

    def __neg__( self ) :

        return( self.returnAsClass( self, pointwiseXY.__neg__( self ) ) )

    def __add__( self, other ) :

        other_ = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )
        return( self.returnAsClass( self, pointwiseXY.__add__( self, other_ ) ) )

    __radd__ = __add__

    def __iadd__( self, other ) :

        other_ = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )
        pointwiseXY.__iadd__( self, other_ )
        return( self )

    def __sub__( self, other ) :

        other_ = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )
        return( self.returnAsClass( self, pointwiseXY.__sub__( self, other_ ) ) )

    def __rsub__( self, other ) :

        sub = self.__sub__( other )
        return( sub.__neg__( ) )

    def __isub__( self, other ) :

        other_ = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )
        pointwiseXY.__isub__( self, other_ )
        return( self )

    def __mul__( self, other ) :

        unit1 = self.getAxisUnitSafely( yAxisIndex )
        other, unit2 = getOtherAndUnit( self, other )
        unit = baseModule.processUnits( unit1, unit2, '*' )
        points = pointwiseXY.__mul__( self, other )
        return( return_pointwiseXY_AsXYs( self, points, units = { yAxisIndex : unit } ) )

    __rmul__ = __mul__

    def __imul__( self, other ) :

        other, otherUnit1 = getOtherAndUnit( self, other )
        if( not( self.mutableYUnit ) ) :
            if( otherUnit1 != '' ) : raise Exception( "Self's y-unit is immutable and other has unit of '%s'" % otherUnit1 )
        pointwiseXY.__imul__( self, other )
        if( otherUnit1 != '' ) : self.axes[yAxisIndex].setUnit( baseModule.processUnits( self.axes[yAxisIndex].unit, otherUnit1, '*' ) )
        return( self )

    def __div__( self, other ) :

        unit1 = self.getAxisUnitSafely( yAxisIndex )
        other, unit2 = getOtherAndUnit( self, other )
        unit = baseModule.processUnits( unit1, unit2, '/' )
        points = pointwiseXY.__div__( self, other )
        return( return_pointwiseXY_AsXYs( self, points, units = { yAxisIndex : unit } ) )

    def __rdiv__( self, other ) :

        unit2 = self.getAxisUnitSafely( yAxisIndex )
        other, unit1 = getOtherAndUnit( self, other )
        points = pointwiseXY.__rdiv__( self, other )
        unit = baseModule.processUnits( unit1, unit2, '/' )
        return( return_pointwiseXY_AsXYs( self, points, units = { yAxisIndex : unit } ) )

    def __idiv__( self, other ) :

        other, otherUnit1 = getOtherAndUnit( self, other )
        if( not( self.mutableYUnit ) ) :
            if( otherUnit1 != '' ) : raise Exception( "Self's y-unit is immutable and other has unit of '%s'" % otherUnit1 )
        pointwiseXY.__idiv__( self, other )
        if( otherUnit1 != '' ) : self.axes[yAxisIndex].setUnit( baseModule.processUnits( self.axes[yAxisIndex].unit, otherUnit1, '/' ) )
        return( self )

    def getitem_units( self, index ) :

        x, y = pointwiseXY.__getitem__( self, index )
        xUnit = self.getAxisUnitSafely( xAxisIndex )
        yUnit = self.getAxisUnitSafely( yAxisIndex )
        return( PQU.PQU( x, xUnit, checkOrder = False ), PQU.PQU( y, xUnit, checkOrder = False ) )

    def __setitem__( self, index, xy ) :

        if( len( xy ) != 2 ) : raise ValueError( 'right-hand-side must be list of length 2 not %s' % len( xy ) )
        pointwiseXY.__setitem__( self, index, xy )

    def setitem_units( self, index, xy ) :

        if( len( xy ) != 2 ) : raise ValueError( 'right-hand-side must be list of length 2 not length %s' % len( xy ) )
        xUnit = self.getAxisUnitSafely( xAxisIndex )
        yUnit = self.getAxisUnitSafely( yAxisIndex )
        xy = [ getValueAsAxis( xUnit, xy[0] ), getValueAsAxis( yUnit, xy[1] ) ]
        pointwiseXY.__setitem__( self, index, xy )

    def __getslice__( self, index1, index2 ) :

        return( self.returnAsClass( self, pointwiseXY.__getslice__( self, index1, index2 ) ) )

    def __setslice__( self, index1, index2, slice ) :

        slice = otherToSelfsUnits( self, slice )
        pointwiseXY.__setslice__( self, index1, index2, slice )

    @property
    def sep( self ) :

        return( self.__sep )

    def applyFunction( self, f, parameters, accuracy = -1, biSectionMax = -1, checkForRoots = False ) :
        """
        This method maps the function 'f' onto each y-value in an XYs object, returning a new XYs object.
        Additional points may be added to preserve the desired accuracy.

        For example, the following will transform all negative y-values in an XYs object to zero:
        >>> newXYs = XYs.applyFunction( lambda y,tmp : 0 if y < 0 else y, None )

        Extra parameters to the applied function are be passed via the 'parameters' argument.
        See XYs.XYs documentation for details on 'accuracy' and 'biSectionMax'.
        """

        return( self.returnAsClass( self, pointwiseXY.applyFunction( self, f, parameters, accuracy = accuracy, biSectionMax = biSectionMax, checkForRoots = checkForRoots ) ) )

    def changeInterpolation( self, interpolation, accuracy = None, lowerEps = 0, upperEps = 0, cls = None ) :

        if( interpolation != standardsModule.interpolation.linlinToken ) : raise ValueError( 'Only "%s" interpolation currently supported: not %s' %
                ( standardsModule.interpolation.linlinToken, interpolation ) )
        if( accuracy is None ) : accuracy = self.getAccuracy( )
        c1 = pointwiseXY.changeInterpolation( self, interpolation = interpolation, accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
        axes = self.axes
        c1 = return_pointwiseXY_AsXYs( self, c1, axes = axes, template = self, value = self.value )
        if( cls is None  ) : cls = self
        return( cls.returnAsClass( self, c1, axes = axes, interpolation = interpolation ) )

    def changeInterpolationIfNeeded( self, allowedInterpolations, accuracy = None, lowerEps = 0, upperEps = 0, cls = None ) :
        """
        If self's interpolation is one in list of allowedInterpolations, self is returned unchaged. Otherwise
        the returned instances is self's data converted to the interpolation allowedInterpolations[0].
        """

        for interpolation in allowedInterpolations :
            if( interpolation == self.interpolation ) : return( self )
        return( self.changeInterpolation( allowedInterpolations[0], accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps, cls = cls ) )

    def clip( self, rangeMin = None, rangeMax = None ) :

        if( rangeMin is None ) : rangeMin = self.rangeMin( )
        if( rangeMax is None ) : rangeMax = self.rangeMax( )
        return( self.returnAsClass( self, pointwiseXY.clip( self, rangeMin, rangeMax ) ) )

    def clip_units( self, rangeMin = None, rangeMax = None ) :

        unit = self.getAxisUnitSafely( yAxisIndex )
        rangeMin, rangeMax = baseModule.getDomainLimits( self, rangeMin, rangeMax, unit )
        return( self.returnAsClass( self, pointwiseXY.clip( self, rangeMin, rangeMax ) ) )

    def commonXGrid( self, others ) :
        """
        This method returns copies of self and others that are mapped to the same X-grid. That is, a union is made
        of all the x-values for self and others and all XYs-type instances are mapped to it. Others must be
        a list of instances of XYs. All XYs instances must have the same domain.
        """

# BRB
        domainMin, domainMax = self.domain( )      # Need to check units. Currently union is raising if not the same.
        grid = self
        for i1, other in enumerate( others ) :
            domainMinO, _domainMaxO = other.domain( )
            if( domainMin != domainMinO ) :
		raise Exception( "domainMin = %e != other's domainMin = %e for other index = %d" % ( domainMin, domainMinO, i1 ) )
            grid = grid.union( other, fillWithSelf = False )
        return( [ o1.union( grid, fillWithSelf = True ) for o1 in [ self ] + others ] )

    def convertAxisToUnit( self, indexOrName, newUnit ) :

        index = self.getAxisIndexByIndexOrName( indexOrName )
        axis = self.axes[index]
        factor = PQU.PQU( '1 ' + axis.unit ).getValueAs( newUnit )
        data = []
        xyIndex = 1 - index
        for xy in self :
            xy[xyIndex] *= factor
            data.append( xy )
        n = return_pointwiseXY_AsXYs( self, data, units = { index : newUnit }, template = self )
        return( self.returnAsClass( self, n, axes = n.axes ) )

    def copy( self, index = None, value = None, axes = None ) :

        xys = pointwiseXY.copy( self )
        if( axes is None ) : axes = self.axes
        return( self.returnAsClass( self, xys, index = index, value = value, axes = axes ) )

    __copy__ = copy
    __deepcopy__ = __copy__

    def copyDataToXYs( self, xUnitTo = None, yUnitTo = None ) :

        xScale, yScale = 1.0, 1.0
        if( xUnitTo is not None ) : xScale = PQU.PQU( '1 ' + self.axes[xAxisIndex].unit ).getValueAs( xUnitTo )
        if( yUnitTo is not None ) : yScale = PQU.PQU( '1 ' + self.axes[yAxisIndex].unit ).getValueAs( yUnitTo )
        return( pointwiseXY.copyDataToXYs( self, xScale = xScale, yScale = yScale ) )

    def copyDataToNestedLists( self, *units ) :

        xUnitTo, yUnitTo = None, None
        if( len( units ) > 0 ) : xUnitTo = units[0]
        if( len( units ) > 1 ) : yUnitTo = units[1]
        return( self.copyDataToXYs( xUnitTo = xUnitTo, yUnitTo = yUnitTo ) )

    def dullEdges( self, lowerEps = 0., upperEps = 0., positiveXOnly = 0 ) :

        d = pointwiseXY.dullEdges( self, lowerEps = lowerEps, upperEps = upperEps, positiveXOnly = positiveXOnly );
        return( self.returnAsClass( self, d ) )

    def getValue_units( self, x ) :

        x = PQU.PQU( x, checkOrder = False ).getValueAs( self.axes[xAxisIndex].unit )
        return( PQU.PQU( pointwiseXY.evaluate( self, x ), self.axes[yAxisIndex].unit, checkOrder = False ) )

    def getValue( self ) :

        return( self.value )

    def evaluate( self, x, unitTo = None ) :

        y = pointwiseXY.evaluate( self, x )
        if( unitTo is None ) : return( y )
        unit = self.getAxisUnitSafely( yAxisIndex )
        return( PQU.PQU( y, unit, checkOrder = False ).getValueAs( unitTo ) )

    def setValue_units( self, x, y ) :

        x = PQU.PQU( x, checkOrder = False ).getValueAs( self.axes[xAxisIndex].unit )
        y = PQU.PQU( y, checkOrder = False ).getValueAs( self.axes[yAxisIndex].unit )
        pointwiseXY.setValue( self, x, y )

    def mutualify( self, lowerEps1, upperEps1, positiveXOnly1, other, lowerEps2, upperEps2, positiveXOnly2 ) :
        '''
        .. note:: Need to check that x units are the same.
        '''
        m1, m2 = pointwiseXY.mutualify( self, lowerEps1, upperEps1, positiveXOnly1, other, lowerEps2, upperEps2, positiveXOnly2 )
        return( self.returnAsClass( self, m1 ), other.returnAsClass( other, m2 ) )

    def normalize( self, insitu = False, dimension = 1 ) :
        """
        The dimension argument is ignored. Only here to be compatable with calling from multiD_XYs.normalize.
        """

        xys = pointwiseXY.normalize( self )
        if( insitu ) :
            pointwiseXY.setData( self, xys )
            return( self )
        return( self.returnAsClass( self, xys ) )

    def thicken( self, sectionSubdivideMax = 1, dDomainMax = 0., fDomainMax = 1. ) :
        """
        .. note:: Need unit for dDomainMax.
        """

# BRB: unit for dDomainMax

        t = pointwiseXY.thicken( self, sectionSubdivideMax = sectionSubdivideMax, dDomainMax = dDomainMax, fDomainMax = fDomainMax )
        return( self.returnAsClass( self, t ) )

    def thin( self, accuracy ) :

        return( self.returnAsClass( self, pointwiseXY.thin( self, accuracy ) ) )

    def trim( self ) :

        return( self.returnAsClass( self, pointwiseXY.trim( self ) ) )

    def union( self, other, fillWithSelf = 1, trim = 0 ) :

        other = otherToSelfsUnits( self, other, checkXOnly = True )
        t = pointwiseXY.union( self, other, fillWithSelf = fillWithSelf, trim = trim  )
        return( self.returnAsClass( self, t ) )

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQU.PQU( '1 ' + self.domainUnit( ) ).getValueAs( unitTo ) )

    def domainMin( self, unitTo = None, asPQU = False ) :

        unit = self.getAxisUnitSafely( xAxisIndex )
        return( PQU.valueOrPQ( pointwiseXY.domainMin( self ), unitFrom = unit, unitTo = unitTo, asPQU = asPQU, checkOrder = False ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        unit = self.getAxisUnitSafely( xAxisIndex )
        return( PQU.valueOrPQ( pointwiseXY.domainMax( self ), unitFrom = unit, unitTo = unitTo, asPQU = asPQU, checkOrder = False ) )

    def domain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainGrid( self, unitTo = None ) :

        scale = self.domainUnitConversionFactor( unitTo )
        return( pointwiseXY.domainGrid( self, scale ) )

    def domainUnit( self ) :

        return( self.getAxisUnitSafely( xAxisIndex ) )

    def rangeMin( self, unitTo = None, asPQU = False ) :

        unit = self.getAxisUnitSafely( yAxisIndex )
        return( PQU.valueOrPQ( pointwiseXY.rangeMin( self ), unitFrom = unit, unitTo = unitTo, asPQU = asPQU, checkOrder = False ) )

    def rangeMax( self, unitTo = None, asPQU = False ) :

        unit = self.getAxisUnitSafely( yAxisIndex )
        return( PQU.valueOrPQ( pointwiseXY.rangeMax( self ), unitFrom = unit, unitTo = unitTo, asPQU = asPQU, checkOrder = False ) )

    def rangeUnit( self ) :

        return( self.getAxisUnitSafely( yAxisIndex ) )

    def domainSlice( self, domainMin = None, domainMax = None, fill = 1, dullEps = 0. ) :
        '''
        Returns a new instance with self sliced between ``domainMin`` and ``domainMax``.

        :param domainMin:   [optional] the lower x-value of the slice, default is domain minimum of self,
        :param domainMax:   [optional] the upper x-value of the slice, default is domain maximum of self,
        :param fill:        [optional] if True, points are added at domainMin and domainMax if they are not in self, 
                                       else only existing points in the range [domainMin, domainMax] are included.
        :param dullEps:     [optional] (Currently not implemented) the lower and upper points are dulled, default is 0.
        '''

        if( domainMin is None ) : domainMin = self.domainMin( )
        if( domainMax is None ) : domainMax = self.domainMax( )
        s = pointwiseXY.domainSlice( self, domainMin = domainMin, domainMax = domainMax, fill = fill, dullEps = dullEps )
        return( self.returnAsClass( self, s ) )

    def domainSlice_units( self, domainMin = None, domainMax = None, fill = 1, dullEps = 0. ) :
        """
        Same as domainSlice, only domainMin and domainMax must be PQU's.
        """

        unit = self.getAxisUnitSafely( xAxisIndex )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        slice = pointwiseXY.domainSlice( self, domainMin = domainMin, domainMax = domainMax, fill = fill, dullEps = dullEps )
        return( self.returnAsClass( self, slice ) )

    def __mod__( self, other ) : raise NotImplementedError( 'Currently, mod is not implemented' )

    def __pow__( self, other ) : raise NotImplementedError( 'Currently, pow is not implemented' )

    def __exp__( self, other ) : raise NotImplementedError( 'Currently, __exp__ is not implemented' )

    def convolute( self, func ) : 
        """
        Uses pointwiseXY's convolute function to do the grunt work.  
        """

        if( not( isinstance( func, pointwiseXY ) ) ) : raise TypeError( 'func argument of convolute() must be an instance of pointwiseXY' )
        return( self.returnAsClass( self, pointwiseXY.convolute( self, func, 0 ) ) )

    def group( self, xs, f2 = None, f3 = None, norm = None, asXYs = False ) :
        """
        The argument ``xs`` is a list of x-values with ``xs[i] < xs[i+1]``. This function calculates the integrals

        .. math::
            \int_{xs[i]}^{xs[i+1]} dx \; { f(x) \over n[i] }

        where :math:`f(x)` is defined as
            +------------------------------------+--------------------+
            | :math:`f(x)`                       |  conditional       |
            +====================================+====================+
            | :math:`f_1(x)`                     | for f2 = f3 = None |
            +------------------------------------+--------------------+
            | :math:`f_1(x) \, f_2(x)`           | for f3 = None      |
            +------------------------------------+--------------------+
            | :math:`f_1(x) \, f_3(x)`           | for f2 = None      |
            +------------------------------------+--------------------+
            | :math:`f_1(x) \, f_2(x) \, f_3(x)` | otherwise          |
            +------------------------------------+--------------------+
        :math:`f_1(x)` is self evaluated at `x` and `n[i]` is a normalization determined from the `norm` arugment as
            +---------------+---------------+
            | norm          | n[i]          |
            +===============+===============+
            | None          | 1             |
            +---------------+---------------+
            | 'dx'          | x[i+1] - x[i] |
            +---------------+---------------+
            | a python list | norm[i]       |
            +---------------+---------------+
        The arguments f2 (and f3) must be None or an XYs instance.

        If ``asXYs`` is ``False``, then ``len( xs ) - 1`` integrals are returned.
        If ``asXYs`` is ``True``, the last integral's values is appended to the end to make a list of length ``len( xs )``, and 
        an instance of class ``XYs`` is returned with the x-values from ``xs``, the y-values from the integrals and the 
        interpolation is 'flat'.

        Historical note: the word group comes from deterministic neutron transport (e.g., transport used to simulate nuclear reactors).
        """

        accuracy, yUnit = self.getAccuracy( ), PQU.PQU( 1, self.getAxisUnitSafely( yAxisIndex ) )
        if( f2 is None ) :
            if( f3 is None ) : 
                grouped = self.groupOneFunction( xs, norm = norm )
            else :
                grouped = self.groupTwoFunctions( xs, f3, norm = norm )
                yUnit = yUnit * PQU.PQU( 1,  f3.getAxisUnitSafely( yAxisIndex ) )
                accuracy = max( accuracy, f3.getAccuracy( ) )
        else :
            yUnit = yUnit * PQU.PQU( 1, f2.getAxisUnitSafely( yAxisIndex ) )
            accuracy = max( accuracy, f2.getAccuracy( ) )
            if( f3 is None ) :
                grouped = self.groupTwoFunctions( xs, f2, norm = norm )
            else :
                grouped = self.groupThreeFunctions( xs, f2, f3, norm = norm )
                yUnit = yUnit * PQU.PQU( 1, f3.getAxisUnitSafely( yAxisIndex ) )
                accuracy = max( accuracy, f3.getAccuracy( ) )
        if( norm is None ) :
            yUnit = PQU.PQU( 1, self.getAxisUnitSafely( xAxisIndex ) ) * yUnit
        elif( norm != 'dx' ) :
            pass                    # Need to add units to norm. That is, norm, grouped and xs should be an instance of Ys.
        if( asXYs ) :
            grouped.append( grouped[-1] )
            axes = axesModule.axes( labelsUnits = { 0 : [ self.axes[xAxisIndex].label, self.axes[xAxisIndex].unit ], 1 : [ "", yUnit.getUnitSymbol( ) ] } )
            grouped = XYs( [ xs, grouped ], dataForm = 'xsandys', accuracy = accuracy, 
                interpolation = standardsModule.interpolation.flatToken, axes = axes )
        return( grouped )

    def groupOneFunction( self, xs, norm = None ) :
        '''.. note:: Need unit of xs.'''

        return( pointwiseXY.groupOneFunction( self, xs, norm = norm ) )

    def groupTwoFunctions( self, xs, f2, norm = None ) :
        '''.. note:: Need unit of xs.'''

        return( pointwiseXY.groupTwoFunctions( self, xs, f2, norm = norm ) )

    def groupThreeFunctions( self, xs, f2, f3, norm = None ) :
        '''.. note:: Need unit of xs.'''

        return( pointwiseXY.groupThreeFunctions( self, xs, f2, f3, norm = norm ) )

    def integrate( self, domainMin = None, domainMax = None ) :
        """
        Definite integral of current ``XYs`` instance from ``domainMin`` to ``domainMax``:
        
        .. math::
            \int_{domainMin}^{domainMax} dx \; XYs(x)

        If ``domainMin`` or ``domainMax`` is unspecified, it is taken from the domain of the self.
        """

        unit = self.getAxisUnitSafely( xAxisIndex )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin( ) )
        domainMax = min( domainMax, self.domainMax( ) )
        return( PQU.PQU( pointwiseXY.integrate( self, domainMin = domainMin, domainMax = domainMax ), 
                baseModule.processUnits( unit, self.getAxisUnitSafely( yAxisIndex ), '*' ), checkOrder = False ) )

    def indefiniteIntegral( self, domainMin = None, domainMax = None ) :
        '''
        Indefinite integral of self:
        
        .. math::
            \int_0^x dx \; XYs(x)
            
        The new ``XYs`` instance is defined on the range of the old one and the units are wrong.
        '''

# BRB: I think this is just a running sum and may be implemented already. Needs units.
        myAxes = self.axes
        myData = [ [ self.domainMin( ), 0.0 ] ]
        for i in range( len( self ) - 1 ):
            domainMin = self[i][0]
            domainMax = self[i+1][0]
            myData.append( [ domainMax, myData[-1][1]+pointwiseXY.integrate( self, domainMin = domainMin, domainMax = domainMax ) ] )
        return XYs( myAxes, myData, 1e-6 )

    def integrateTwoFunctions( self, f2, domainMin = None, domainMax = None ) :
        """

        :param f2:
        :param domainMin:
        :param domainMax:
        :return:
        """
        if not isinstance(f2,XYs): raise TypeError("f2 must be an instance of an XYs")
        unit = self.getAxisUnitSafely( xAxisIndex )
        if f2.getAxisUnitSafely( xAxisIndex ) != unit: f2.convertAxisToUnit( xAxisIndex, unit )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin( ), f2.domainMin( ) )
        domainMax = min( domainMax, self.domainMax( ), f2.domainMax( ) )
        return( PQU.PQU( pointwiseXY.groupTwoFunctions( self, [ domainMin, domainMax ], f2 )[0],
                baseModule.processUnits( baseModule.processUnits( unit, self.getAxisUnitSafely( yAxisIndex ), '*' ), f2.getAxisUnitSafely( yAxisIndex ), '*' ), checkOrder = False ) )

    def integrateThreeFunctions( self, f2, f3, domainMin = None, domainMax = None ) :
        """

        :param f2:
        :param f3:
        :param domainMin:
        :param domainMax:
        :return:
        """
        if not isinstance(f2,XYs): raise TypeError("f2 must be an instance of an XYs")
        if not isinstance(f3,XYs): raise TypeError("f3 must be an instance of an XYs")
        unit = self.getAxisUnitSafely( xAxisIndex )
        if f2.getAxisUnitSafely( xAxisIndex ) != unit: f2.convertAxisToUnit( xAxisIndex, unit )
        if f3.getAxisUnitSafely( xAxisIndex ) != unit: f3.convertAxisToUnit( xAxisIndex, unit )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin( ), f2.domainMin( ), f3.domainMin( ) )
        domainMax = min( domainMax, self.domainMax( ), f2.domainMax( ), f3.domainMax( ) )
        return( pointwiseXY.groupThreeFunctions( self, [ domainMin, domainMax ], f2, f3 )[0] )

    def scaleDependent( self, value, insitu = False ) :

        xys = self
        if( not( insitu ) ) : xys = self.copy( )
        xys.scaleOffsetXAndY( yScale = value, insitu = True )
        return( xys )

    def splitInTwo( self, domainValue, epsilon = domainEpsilon ) : 
        """
        Splits self at domainValue into two XYs instances. Returns either None if domainValue not in
        domain or two XYs instances.
        """

        domainMin, domainMax = self.domain( )
        if( domainValue <= ( 1 + epsilon ) * domainMin ) : return( None )
        if( domainValue >= ( 1 - epsilon ) * domainMax ) : return( None )
        return( self.domainSlice( domainMax = domainValue ), self.domainSlice( domainMin = domainValue ) )

    def toPointwise_withLinearXYs( self, accuracy = None, lowerEps = 0, upperEps = 0, cls = None ) :

        if( accuracy is None ) : accuracy = self.getAccuracy( )
        return( self.changeInterpolation( standardsModule.interpolation.linlinToken, accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps, cls = cls ) )

    def toXMLList( self, indent = '', **kwargs ) :

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        oneLine = kwargs.get( 'oneLine', False )

        indent2 = indent + incrementalIndent
        if( oneLine ) : indent2 = ''

        attributeStr = baseModule.xDataFunctional.attributesToXMLAttributeStr( self )
        accuracy = self.getAccuracy( )
        if( accuracy != defaultAccuracy ) : attributeStr += ' accuracy="%s"' % accuracy
        if( self.interpolation != standardsModule.interpolation.linlinToken ) : attributeStr += ' interpolation="%s"' % self.interpolation

        XMLList = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ] 
        if( self.isPrimaryXData( ) ) : 
            if( self.axes is not None ) : XMLList += self.axes.toXMLList( indent = indent2, **kwargs )
        xys = []
        for x, y in self.copyDataToXYs( ) :
            xys.append( x )
            xys.append( y )
        XMLList += valuesModule.values( xys, valueType = self.valueType, sep = self.__sep ).toXMLList( indent2, **kwargs )
        if( self.uncertainties ) : XMLList += self.uncertainties.toXMLList( indent = indent2, **kwargs )
        XMLList[-1] += '</%s>' % self.moniker
        if( oneLine ) : return( [ ''.join( XMLList ) ] )
        return( XMLList )

    def tweakDomain( self, domainMin = None, domainMax = None, epsilon = domainEpsilon ) :

        if( len( self ) == 0 ) : return
        if( domainMin is not None ) : 
            x, y = self[0]
            if( abs( domainMin - x ) < epsilon * ( max( abs( domainMin ), abs( x ) ) ) ) : self[0] = domainMin, y
        if( domainMax is not None ) : 
            x, y = self[-1]
            if( abs( domainMax - x ) < epsilon * ( max( abs( domainMax ), abs( x ) ) ) ) : self[-1] = domainMax, y

    def plot( self, xylog = 0, domainMin = None, domainMax = None, rangeMin = None , rangeMax = None, title = '' ) :

        import subprocess, os
        from fudge.core.utilities import fudgeFileMisc
        from fudge.vis.gnuplot import plotbase

        def getUnitlessNumber( value, unit, default ) :

            if( value is None ) : return( default )
            if( not( isinstance( value, PQU.PQU ) ) ) : value = PQU.PQU( value )
            return( value.getValueAs( unit ) )

        if( self.axes is None ) :
            xUnit, yUnit = '', ''
            xLabel, yLabel = '', ''
        else :
            xUnit = self.getAxisUnitSafely( xAxisIndex )
            yUnit = self.getAxisUnitSafely( yAxisIndex )
            xLabel = self.axes[xAxisIndex].plotLabel( )
            yLabel = self.axes[yAxisIndex].plotLabel( )

        domainMin = getUnitlessNumber( domainMin, xUnit, self.domainMin( ) )
        domainMax = getUnitlessNumber( domainMax, xUnit, self.domainMax( ) )
        rangeMin = getUnitlessNumber( rangeMin, yUnit, self.rangeMin( ) )
        rangeMax = getUnitlessNumber( rangeMax, yUnit, self.rangeMax( ) )

        dt = plotbase.parsePlotOptions( domainMin, domainMax, rangeMin, rangeMax, xLabel, yLabel, title )
        f = fudgeFileMisc.fudgeTempFile( )
        f.write( self.toString( ) )
        f.close( )
        p = os.path.join( os.path.realpath( __file__ ).split( '/xData' )[0], "fudge", "vis", "gnuplot", "endl2dplot.py" )
        args = [ "python", p, 'xylog', str( xylog ) ] + dt + [ f.getName( ) ]
        subprocess.Popen( args )

    @classmethod
    def returnAsClass( cls, self, other, index = None, value = None, axes = None, interpolation = None ) :
        """
        Returns other as a class of cls. Other must be a sub-class pointwiseXY. cls must be a sub-class of XYs. 
        If index and axes are not specified, they are taken from self. The main use of this classmethod is for methods like __add__ where
        the addends may be a class derived from XYs. For example, the crossSection.pointwise class is derived from The XYs class.
        If two instances xSec1 and xSec2 of the crossSection.pointwise class are added together (i.e., xSec1 + xSec2) then, since
        the __add__ method used this classmethod, the returned instance will also be an instance of crossSection.pointwise class.
        """

        if( index is None ) : index = self.index
        if( value is None ) : value = self.value
        if( axes is None ) : axes = self.axes
        if( interpolation is None ) : interpolation = self.interpolation
        return( cls( data = other, accuracy = other.getAccuracy( ), interpolation = interpolation, axes = axes, 
                overflowSize = 10, biSectionMax = other.getBiSectionMax( ), infill = other.getInfill( ), 
                safeDivide = other.getSafeDivide( ), index = index, value = value ) )

    @classmethod
    def parseXMLNode( cls, xDataElement, xPath = [], linkData = {}, axes = None, **kwargs ) :
        """
        Translate an XYs XML element into the python XYs xData class.
        """

        xPath.append( xDataElement.tag )

        attrs = { 'accuracy' : defaultAccuracy, 'interpolation' : standardsModule.interpolation.linlinToken, 'label' : None, 'index' : None, 'value' : None, 
                'biSectionMax' : 6 }
        attributes = { 'accuracy' : float, 'interpolation' : str, 'label' : str, 'index' : int, 'value' : float }
        for key, item in xDataElement.items( ) :
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key]( item )

        values, uncertainties = None, None
        for subElement in xDataElement :
            if( subElement.tag == 'axes' ) :
                axes = axesModule.axes.parseXMLNode( subElement, xPath )
            elif( subElement.tag == 'values' ) :
                values = valuesModule.values.parseXMLNode( subElement, xPath )
            elif( subElement.tag == 'uncertainties' ) :
                uncertainties = uncertaintiesModule.uncertainties.parseXMLNode( subElement, xPath )
            else :
                raise TypeError( 'sub-element "%s" not valid' % subElement.tag )

        if( values is None ) : raise Exception( 'values element missing' )
        attrs['sep'] = values.sep

        xys = cls( data = values, dataForm = "list", axes = axes, **attrs )
        if uncertainties is not None: xys.uncertainties = uncertainties
        xPath.pop( )
        return( xys )

    @classmethod
    def parseXMLString( cls, XMLString ) :

        from xml.etree import cElementTree

        return( cls.parseXMLNode( cElementTree.fromstring( XMLString ) ) )

    @staticmethod
    def createFromFunction( axes, Xs, func, parameters, accuracy, biSectionMax, checkForRoots = False, infill = 1, safeDivide = 1 ) :
        """
        Given an ascending list of x-values (Xs) and a function (func), create a linear pointwise 
        representation of the function over the domain of Xs. There must be at least 2 values in Xs.
        See pointwiseXY_C.createFromFunction for all other arguments.
        """

        import math
        xys = pointwiseXY_C.createFromFunction( Xs, func, parameters, accuracy, biSectionMax, checkForRoots = checkForRoots, infill = infill, safeDivide = safeDivide )
        biSectionMax = max( 0, biSectionMax - math.log( len( xys ) / len( Xs ) ) / math.log( 2 ) )
        return( XYs( xys, axes = axes, accuracy = accuracy, infill = infill, safeDivide = safeDivide ) )

    @staticmethod
    def defaultAxes( labelsUnits = {} ) :

        return( axesModule.axes( rank = 2, labelsUnits = labelsUnits ) )
