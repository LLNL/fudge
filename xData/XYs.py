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
    This module contains the class ``XYs1d``. This class treats a list of :math:`(x_i, y_i)` pairs as if they were 
    the function :math:`y(x)`.  That is, it is a numerical representation of :math:`f(x)`. The function :math:`y(x)` 
    is also called a 1-dimensional or univariant function. As an example, let 
    
        .. math::
            y_1(x) = a + b * x

    and

        .. math::
            y_2(x) = c * \exp( - d * x )
            
    then, :math:`y_1` and :math:`y_2` can be added, subtracted, multiplied and divided, for example. Similiarly, 
    two XYs1d instances can be added, subtracted, multiplied and divided. The only restriction on the :math:`(x_i, y_i)` 
    pairs is that :math:`x_i < x_{i+1}`.

    The ``XYs1d`` class uses the class :py:class:`numericalFunctions.lib.pointwiseXY_C.pointwiseXY_C` as a base class.

    Members of an XYs1d instance of interest to most are:

        :infill:          see the base class :py:class:`numericalFunctions.lib.pointwiseXY_C.pointwiseXY_C`,
        :axes:            Description of the x and y data attributes (e.g., label, units).

    A ``XYs1d`` instance can be the ``XY`` instance of a 2-dimensional function (i.e., a ``XYs2d`` with
    dimension of 2). In this case, the ``XYs1d`` instances are secondary instances.
"""

"""
Notes:
    1) plot method still using fudge stuff.
    2) need list of values with unit. Used, for example, for method groupOneFunction.
"""

__metaclass__ = type

import math

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

def return_pointwiseXY_AsXYs( self, C, units = {}, index = None, value = None, axes = None ) :

    if( index is None ) : index = self.index
    if( axes is None ) : axes = self.axes
    c = XYs1d( C, axes = axes, infill = True, safeDivide = False, index = index, value = value )
    if( c.axes is not None ) :
        for k in units : c.axes[k].unit = units[k]
    return( c )

def otherToSelfsUnits( self, other, checkXOnly = False ) :

    if( not( isinstance( other, XYs1d ) ) ) : raise TypeError( 'other instance not XYs1d instance' )
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
    For internal use only. This function is used by the multiply and divide methods. Other must be a XYs1d instance or
    an object convertible to a PQU object.
    """

    if( isinstance( other, XYs1d ) ) :
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

    if( isinstance( other, XYs1d ) ) :
         other = otherToSelfsUnits( self, other )
    else :
        yUnit = self.getAxisUnitSafely( yAxisIndex )
        other = PQU.PQU( other, checkOrder = False ).getValueAs( yUnit )
    return( other )

class XYs1d( pointwiseXY, baseModule.xDataFunctional ) :

    moniker = 'XYs1d'
    dimension = 1
    mutableYUnit = True     # For __imul__ and __idiv__.

    def __init__( self, data, dataForm = "xys", interpolation = standardsModule.interpolation.linlinToken, axes = None,
            index = None, valueType = standardsModule.types.float64Token, value = None, label = None, 
            sep = ' ', initialSize = 10, overflowSize = 10, infill = True, safeDivide = False ) :
        """
        Constructor for XYs1d class. dataForm can be 'xys', 'xsandys' or 'list'.
        """

        baseModule.xDataFunctional.__init__( self, self.moniker, axes, index = index, valueType = valueType,
                value = value, label = label )

        if( not( isinstance( interpolation, str ) ) ) : raise TypeError( 'interpolation must be a string' )

        if( not( isinstance( sep, str ) ) ) : raise TypeError( 'sep must be of type str' )
        if( len( sep ) != 1 ) : raise TypeError( 'sep length must be 1 not %d' % len( sep ) )
        self.__sep = sep

        initialSize = max( initialSize, len( data ) )
        pointwiseXY.__init__( self, data = data, dataForm = dataForm, initialSize = initialSize, overflowSize = overflowSize,
            accuracy = defaultAccuracy, biSectionMax = 16, interpolation = interpolation, infill = infill, 
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
        if( otherUnit1 != '' ) : self.axes[yAxisIndex].unit = baseModule.processUnits( self.axes[yAxisIndex].unit, otherUnit1, '*' )
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
        if( otherUnit1 != '' ) : self.axes[yAxisIndex].unit = baseModule.processUnits( self.axes[yAxisIndex].unit, otherUnit1, '/' )
        return( self )

    def __setitem__( self, index, xy ) :

        if( len( xy ) != 2 ) : raise ValueError( 'right-hand-side must be list of length 2 not %s' % len( xy ) )
        pointwiseXY.__setitem__( self, index, xy )

    def __getslice__( self, index1, index2 ) :

        return( self.returnAsClass( self, pointwiseXY.__getslice__( self, index1, index2 ) ) )

    def __setslice__( self, index1, index2, slice ) :

        slice = otherToSelfsUnits( self, slice )
        pointwiseXY.__setslice__( self, index1, index2, slice )

    @property
    def interpolation( self ) :

        return( self.getInterpolation( ) )
#
# This property needs setInterpolation defined in pointwiseXY_C
#    @interpolation.setter
#    def interpolation( self, interpolation ) :
#
#        self.setInterpolation( interpolation )

    @property
    def sep( self ) :

        return( self.__sep )

    def applyFunction( self, f, parameters, accuracy = -1, biSectionMax = -1, checkForRoots = False ) :
        """
        This method maps the function 'f' onto each y-value in an XYs1d object, returning a new XYs1d object.
        Additional points may be added to preserve the desired accuracy.

        For example, the following will transform all negative y-values in an XYs1d object to zero:
        >>> newXYs = XYs1d.applyFunction( lambda y,tmp : 0 if y < 0 else y, None )

        Extra parameters to the applied function are be passed via the 'parameters' argument.
        See XYs.XYs1d documentation for details on 'accuracy' and 'biSectionMax'.
        """

        return( self.returnAsClass( self, pointwiseXY.applyFunction( self, f, parameters, accuracy = accuracy, biSectionMax = biSectionMax, checkForRoots = checkForRoots ) ) )

    def changeInterpolation( self, interpolation, accuracy, lowerEps = 0, upperEps = 0, cls = None ) :

        if( interpolation != standardsModule.interpolation.linlinToken ) : raise ValueError( 'Only "%s" interpolation currently supported: not %s' %
                ( standardsModule.interpolation.linlinToken, interpolation ) )
        c1 = pointwiseXY.changeInterpolation( self, interpolation = interpolation, accuracy = accuracy, lowerEps = lowerEps, 
                upperEps = upperEps )
        axes = self.axes
        c1 = return_pointwiseXY_AsXYs( self, c1, axes = axes, value = self.value )
        if( cls is None  ) : cls = self
        return( cls.returnAsClass( self, c1, axes = axes, interpolation = interpolation ) )

    def changeInterpolationIfNeeded( self, allowedInterpolations, accuracy, lowerEps = 0, upperEps = 0, cls = None ) :
        """
        If self's interpolation is one in list of allowedInterpolations, self is returned unchaged. Otherwise
        the returned instances is self's data converted to the interpolation allowedInterpolations[0].
        """

        for interpolation in allowedInterpolations :
            if( interpolation == self.interpolation ) : return( self )
        return( self.changeInterpolation( allowedInterpolations[0], accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps, cls = cls ) )

    def cloneToInterpolation( self, interpolation ) :

        c1 = pointwiseXY.cloneToInterpolation( self, interpolation )
        return( self.__class__.returnAsClass( self, c1, axes = self.axes, interpolation = interpolation ) )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        if( self.axes is None ) : print self.toXLink( )
        factors = self.axes.convertUnits( unitMap )
        if( factors[:2] != [ 1., 1. ] ) : self.scaleOffsetXAndY( xScale = factors[1], yScale = factors[0], insitu = True )
        self.fixValuePerUnitChange( factors )

    def clip( self, rangeMin = None, rangeMax = None ) :

        if( rangeMin is None ) : rangeMin = self.rangeMin
        if( rangeMax is None ) : rangeMax = self.rangeMax
        return( self.returnAsClass( self, pointwiseXY.clip( self, rangeMin, rangeMax ) ) )

    def commonDomainGrid( self, others ) :
        """
        This method returns copies of self and others that are mapped to the same X-grid. That is, a union is made
        of all the x-values for self and others and all XYs1d-type instances are mapped to it. Others must be
        a list of instances of XYs1d. All XYs1d instances must have the same domain.
        """

# BRB
        domainMin = self.domainMin                  # Need to check units. Currently union is raising if not the same.
        grid = self
        for i1, other in enumerate( others ) :
            domainMinO = other.domainMin
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
        n = return_pointwiseXY_AsXYs( self, data, units = { index : newUnit } )
        return( self.returnAsClass( self, n, axes = n.axes ) )

    def copy( self ) :

        xys = pointwiseXY.copy( self )
        axes = self.axes
        if( axes is not None ) : axes = axes.copy( )
        return( self.returnAsClass( self, xys, index = self.index, value = self.value, axes = axes ) )

    __copy__ = copy
    __deepcopy__ = __copy__

    def copyDataToNestedLists( self ) :

        return( self.copyDataToXYs( ) )

    def dullEdges( self, lowerEps = 0., upperEps = 0., positiveXOnly = 0 ) :

        d = pointwiseXY.dullEdges( self, lowerEps = lowerEps, upperEps = upperEps, positiveXOnly = positiveXOnly );
        return( self.returnAsClass( self, d ) )

    def evaluate( self, x ) :

        return( pointwiseXY.evaluate( self, x ) )

    def mutualify( self, lowerEps1, upperEps1, positiveXOnly1, other, lowerEps2, upperEps2, positiveXOnly2 ) :
        '''
        .. note:: Need to check that x units are the same.
        '''
        m1, m2 = pointwiseXY.mutualify( self, lowerEps1, upperEps1, positiveXOnly1, other, lowerEps2, upperEps2, positiveXOnly2 )
        return( self.returnAsClass( self, m1 ), other.returnAsClass( other, m2 ) )

    def normalize( self, insitu = False, dimension = 1 ) :
        """
        The dimension argument is ignored. Only here to be compatable with calling from XYsnd.normalize.
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

        t = pointwiseXY.thicken( self, sectionSubdivideMax = sectionSubdivideMax, dDomainMax = dDomainMax, fDomainMax = fDomainMax )
        return( self.returnAsClass( self, t ) )

    def thin( self, accuracy ) :

        return( self.returnAsClass( self, pointwiseXY.thin( self, accuracy ) ) )

    def thinToNumberOfPoints( self, maximumNumber, numberFraction = 0.95 ) :
        """
        Returns a tuple containing an accuracy and a new instance of self that hopefully has no more than 
        maximumNumber points.  The returned accuracy is the one passed to the thin method that obtained the 
        returned XYs1d instance.  No thinning is performed if the number of points in self is less than 6.
        A scan of accuracies is performed until the number of points is between minimumNumber and maximumNumber.
        minimumNumber is calculated as int( maximumNumber * numberFraction ) with the restriction that
        minimumNumber must be at least 2 less than maximumNumber. numberFraction is limited to be
        between 0.75 and 0.98. If accuracy is negative, self length is not greater than maximumNumber
        and a copy of self is returned as the XYs1d instance.

        Note, there is no guarantee that the returned instance will have less than maximumNumber points.
        The caller should check the returned instance's length.  For example, if self contains y-values 
        that oscillate between 1 and -1, no thinning is possible.
        """

        if( ( len( self ) > maximumNumber ) and ( len( self ) > 5 ) ) :
            numberFraction = min( .98, max( 0.75, numberFraction ) )
            minimumNumber = int( numberFraction * maximumNumber )
            if( minimumNumber > ( maximumNumber - 2 ) ) : minimumNumber = maximumNumber - 2

            accuracyMin = 1e-12
            accuracyMax = 1

            minSelf = self.thin( accuracyMin )
            maxSelf = self.thin( accuracyMax )
            for i1 in range( 10 ) :
                accuracyMid = math.sqrt( accuracyMin * accuracyMax )
                midSelf = self.thin( accuracyMid )
                if( len( midSelf ) < maximumNumber ) :
                    if( len( midSelf ) > minimumNumber ) : break
                    accuracyMax = accuracyMid
                else :
                    accuracyMin = accuracyMid
        else :
            accuracyMid = -1
            midSelf = self.copy( )

        return( accuracyMid, self.returnAsClass( self, midSelf ) )

    def trim( self ) :

        return( self.returnAsClass( self, pointwiseXY.trim( self ) ) )

    def union( self, other, fillWithSelf = 1, trim = 0 ) :

        other = otherToSelfsUnits( self, other, checkXOnly = True )
        t = pointwiseXY.union( self, other, fillWithSelf = fillWithSelf, trim = trim  )
        return( self.returnAsClass( self, t ) )

    @property
    def domainMin( self ) :

        return( pointwiseXY.domainMin( self ) )

    @property
    def domainMax( self ) :

        return( pointwiseXY.domainMax( self ) )

    @property
    def domainUnit( self ) :

        return( self.getAxisUnitSafely( xAxisIndex ) )

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQU.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    @property
    def domainGrid( self ) :

        return( pointwiseXY.domainGrid( self, 1 ) )

    @property
    def rangeMin( self ) :

        return( pointwiseXY.rangeMin( self ) )

    @property
    def rangeMax( self ) :

        return( pointwiseXY.rangeMax( self ) )

    @property
    def rangeUnit( self ) :

        return( self.getAxisUnitSafely( yAxisIndex ) )

    def rangeUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQU.PQU( '1 ' + self.rangeUnit ).getValueAs( unitTo ) )

    def domainSlice( self, domainMin = None, domainMax = None, fill = 1, dullEps = 0. ) :
        '''
        Returns a new instance with self sliced between ``domainMin`` and ``domainMax``.

        :param domainMin:   [optional] the lower x-value of the slice, default is domain minimum of self,
        :param domainMax:   [optional] the upper x-value of the slice, default is domain maximum of self,
        :param fill:        [optional] if True, points are added at domainMin and domainMax if they are not in self, 
                                       else only existing points in the range [domainMin, domainMax] are included.
        :param dullEps:     [optional] (Currently not implemented) the lower and upper points are dulled, default is 0.
        '''

        if( domainMin is None ) : domainMin = self.domainMin
        if( domainMax is None ) : domainMax = self.domainMax
        s = pointwiseXY.domainSlice( self, domainMin = domainMin, domainMax = domainMax, fill = fill, dullEps = dullEps )
        return( self.returnAsClass( self, s ) )

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
        The arguments f2 (and f3) must be None or an XYs1d instance.

        If ``asXYs`` is ``False``, then ``len( xs ) - 1`` integrals are returned.
        If ``asXYs`` is ``True``, the last integral's values is appended to the end to make a list of length ``len( xs )``, and 
        an instance of class ``XYs1d`` is returned with the x-values from ``xs``, the y-values from the integrals and the 
        interpolation is 'flat'.

        Historical note: the word group comes from deterministic neutron transport (e.g., transport used to simulate nuclear reactors).
        """

        yUnit = PQU.PQU( 1, self.getAxisUnitSafely( yAxisIndex ) )
        if( f2 is None ) :
            if( f3 is None ) : 
                grouped = self.groupOneFunction( xs, norm = norm )
            else :
                grouped = self.groupTwoFunctions( xs, f3, norm = norm )
                yUnit = yUnit * PQU.PQU( 1,  f3.getAxisUnitSafely( yAxisIndex ) )
        else :
            yUnit = yUnit * PQU.PQU( 1, f2.getAxisUnitSafely( yAxisIndex ) )
            if( f3 is None ) :
                grouped = self.groupTwoFunctions( xs, f2, norm = norm )
            else :
                grouped = self.groupThreeFunctions( xs, f2, f3, norm = norm )
                yUnit = yUnit * PQU.PQU( 1, f3.getAxisUnitSafely( yAxisIndex ) )
        if( norm is None ) :
            yUnit = PQU.PQU( 1, self.getAxisUnitSafely( xAxisIndex ) ) * yUnit
        elif( norm != 'dx' ) :
            pass                    # Need to add units to norm. That is, norm, grouped and xs should be an instance of Ys.
        if( asXYs ) :
            grouped.append( grouped[-1] )
# BRB: FIXME, the next line probably had indicies reversed.
            axes = axesModule.axes( labelsUnits = { 0 : [ self.axes[xAxisIndex].label, self.axes[xAxisIndex].unit ], 1 : [ "", yUnit.getUnitSymbol( ) ] } )
            grouped = XYs1d( [ xs, grouped ], dataForm = 'xsandys', 
                interpolation = standardsModule.interpolation.flatToken, axes = axes )
        return( grouped )

    def groupOneFunction( self, xs, norm = None ) :
        '''.. note:: Need unit of xs.'''

        if type(xs)==list:          boundaries = xs
        elif type(xs.values)==list: boundaries = xs.values
        else:                       boundaries = xs.values.values
        return( pointwiseXY.groupOneFunction( self, boundaries, norm = norm ) )

    def groupTwoFunctions( self, xs, f2, norm = None ) :
        '''.. note:: Need unit of xs.'''

        if type(xs)==list:          boundaries = xs
        elif type(xs.values)==list: boundaries = xs.values
        else:                       boundaries = xs.values.values
        return( pointwiseXY.groupTwoFunctions( self, boundaries, f2, norm = norm ) )

    def groupThreeFunctions( self, xs, f2, f3, norm = None ) :
        '''.. note:: Need unit of xs.'''

        if type(xs)==list:          boundaries = xs
        elif type(xs.values)==list: boundaries = xs.values
        else:                       boundaries = xs.values.values
        return( pointwiseXY.groupThreeFunctions( self, boundaries, f2, f3, norm = norm ) )

    def integrate( self, domainMin = None, domainMax = None ) :
        """
        Definite integral of current ``XYs1d`` instance from ``domainMin`` to ``domainMax``:
        
        .. math::
            \int_{domainMin}^{domainMax} dx \; XYs(x)

        If ``domainMin`` or ``domainMax`` is unspecified, it is taken from the domain of the self.
        """

        unit = self.getAxisUnitSafely( xAxisIndex )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin )
        domainMax = min( domainMax, self.domainMax )
        return( PQU.PQU( pointwiseXY.integrate( self, domainMin = domainMin, domainMax = domainMax ), 
                baseModule.processUnits( unit, self.getAxisUnitSafely( yAxisIndex ), '*' ), checkOrder = False ) )

    def indefiniteIntegral( self, domainMin = None, domainMax = None ) :
        '''
        Indefinite integral of self:
        
        .. math::
            \int_0^x dx \; XYs(x)
            
        The new ``XYs1d`` instance is defined on the range of the old one and the units are wrong.
        '''

# BRB: I think this is just a running sum and may be implemented already. Needs units.
        myAxes = self.axes
        myData = [ [ self.domainMin, 0.0 ] ]
        for i in range( len( self ) - 1 ):
            domainMin = self[i][0]
            domainMax = self[i+1][0]
            myData.append( [ domainMax, myData[-1][1]+pointwiseXY.integrate( self, domainMin = domainMin, domainMax = domainMax ) ] )
        return XYs1d( myData, axes = myAxes )

    def integrateTwoFunctions( self, f2, domainMin = None, domainMax = None ) :
        """

        :param f2:
        :param domainMin:
        :param domainMax:
        :return:
        """

        if( not isinstance( f2, XYs1d ) ) : raise TypeError( "f2 must be an instance of an XYs1d" )
        unit = self.axes[xAxisIndex].unit
        if( f2.axes[xAxisIndex].unit != unit ) :
            f2 = f2.copy( )
            f2.convertAxisToUnit( xAxisIndex, unit )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin, f2.domainMin )
        domainMax = min( domainMax, self.domainMax, f2.domainMax )
        return( PQU.PQU( pointwiseXY.groupTwoFunctions( self, [ domainMin, domainMax ], f2 )[0],
                baseModule.processUnits( baseModule.processUnits( unit, self.axes[yAxisIndex].unit, '*' ), f2.axes[yAxisIndex].unit, '*' ), checkOrder = False ) )

    def integrateThreeFunctions( self, f2, f3, domainMin = None, domainMax = None ) :
        """

        :param f2:
        :param f3:
        :param domainMin:
        :param domainMax:
        :return:
        """
        if( not isinstance( f2, XYs1d ) ) : raise TypeError( "f2 must be an instance of an XYs1d" )
        if( not isinstance( f3, XYs1d ) ) : raise TypeError( "f3 must be an instance of an XYs1d" )
        unit = self.axes[xAxisIndex].unit
        if( f2.axes[xAxisIndex].unit != unit ) :
            f2 = f2.copy( )
            f2.convertAxisToUnit( xAxisIndex, unit )
        if( f3.axes[xAxisIndex].unit != unit ) :
            f3 = f3.copy( )
            f3.convertAxisToUnit( xAxisIndex, unit )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin, f2.domainMin, f3.domainMin )
        domainMax = min( domainMax, self.domainMax, f2.domainMax, f3.domainMax )
        return( pointwiseXY.groupThreeFunctions( self, [ domainMin, domainMax ], f2, f3 )[0] )

    def scaleDependent( self, value, insitu = False ) :

        xys = self
        if( not( insitu ) ) : xys = self.copy( )
        xys.scaleOffsetXAndY( yScale = value, insitu = True )
        return( xys )

    def splitInTwo( self, domainValue, epsilon = domainEpsilon ) : 
        """
        Splits self at domainValue into two XYs1d instances. Returns either None if domainValue not in
        domain or two XYs1d instances.
        """

        domainMin, domainMax = self.domainMin, self.domainMax
        if( domainValue <= ( 1 + epsilon ) * domainMin ) : return( None )
        if( domainValue >= ( 1 - epsilon ) * domainMax ) : return( None )
        return( self.domainSlice( domainMax = domainValue ), self.domainSlice( domainMin = domainValue ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Returns a new instance, converted to lin-lin interpolation with added points to maintain desired accuracy.

        Optional (key-word) arguments:
        :param accuracy: desired accuracy. Controls how many points are added when switching interpolation
        :param lowerEps: has no effect, kept for compatibility with regions1d.toPointwise_withLinearXYs
        :param upperEps: has no effect, kept for compatibility with regions1d.toPointwise_withLinearXYs
        :param cls: class to return. Defaults to xData.XYs.XYs1d
        :return:
        """

        arguments = self.getArguments( kwargs, { 'accuracy' : defaultAccuracy, 'lowerEps' : 0, 'upperEps' : 0, 'cls' : None } )
        accuracy = arguments['accuracy']
        lowerEps = arguments['lowerEps']
        upperEps = arguments['upperEps']
        cls = arguments['cls']
        return( self.changeInterpolation( standardsModule.interpolation.linlinToken, accuracy, lowerEps = lowerEps, 
                upperEps = upperEps, cls = cls ) )

    def toXMLList( self, indent = '', **kwargs ) :

        incrementalIndent = kwargs.get( 'incrementalIndent', '  ' )
        oneLine = kwargs.get( 'oneLine', False )

        indent2 = indent + incrementalIndent
        if( oneLine ) : indent2 = ''

        attributeStr = baseModule.xDataFunctional.attributesToXMLAttributeStr( self )
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

        domainMin = getUnitlessNumber( domainMin, xUnit, self.domainMin )
        domainMax = getUnitlessNumber( domainMax, xUnit, self.domainMax )
        rangeMin = getUnitlessNumber( rangeMin, yUnit, self.rangeMin )
        rangeMax = getUnitlessNumber( rangeMax, yUnit, self.rangeMax )

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
        Returns other as a class of cls. Other must be a sub-class pointwiseXY. cls must be a sub-class of XYs1d. 
        If index and axes are not specified, they are taken from self. The main use of this classmethod is for 
        methods like __add__ where the addends may be a class derived from XYs1d. For example, the crossSection.XYs1d
        class is derived from The XYs1d class.  If two instances xSec1 and xSec2 of the crossSection.pointwise class 
        are added together (i.e., xSec1 + xSec2) then, since the __add__ method used this classmethod, the returned 
        instance will also be an instance of crossSection.pointwise class.
        """

        if( index is None ) : index = self.index
        if( value is None ) : value = self.value
        if( axes is None ) : axes = self.axes
        if( interpolation is None ) : interpolation = self.interpolation
        return( cls( data = other, interpolation = interpolation, axes = axes, 
                overflowSize = 10, infill = other.getInfill( ), 
                safeDivide = other.getSafeDivide( ), index = index, value = value ) )

    @classmethod
    def parseXMLNode( cls, xDataElement, xPath, linkData, axes = None, **kwargs ) :
        """
        Translate an XYs1d XML element into the python XYs1d xData class.
        """

        xPath.append( xDataElement.tag )

        attrs = { 'interpolation' : standardsModule.interpolation.linlinToken, 'label' : None, 
                'index' : None, 'value' : None }
        attributes = { 'interpolation' : str, 'label' : str, 'index' : int, 'value' : float }
        for key, item in xDataElement.items( ) :
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key]( item )

        values, uncertainties = None, None
        for subElement in xDataElement :
            if( subElement.tag == 'axes' ) :
                axes = axesModule.axes.parseXMLNode( subElement, xPath, linkData )
            elif( subElement.tag == 'values' ) :
                values = valuesModule.values.parseXMLNode( subElement, xPath, linkData )
            elif( subElement.tag == 'uncertainties' ) :
                uncertainties = uncertaintiesModule.uncertainties.parseXMLNode( subElement, xPath, linkData )
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

        return( cls.parseXMLNode( cElementTree.fromstring( XMLString ), xPath=[], linkData={} ) )

    @classmethod
    def createFromFunction( cls, axes, Xs, func, parameters, accuracy, biSectionMax, checkForRoots = False, infill = 1, safeDivide = 1 ) :
        """
        Given an ascending list of x-values (Xs) and a function (func), create a linear pointwise 
        representation of the function over the domain of Xs. There must be at least 2 values in Xs.
        See pointwiseXY_C.createFromFunction for all other arguments.
        """

        xys = pointwiseXY_C.createFromFunction( Xs, func, parameters, accuracy, biSectionMax, checkForRoots = checkForRoots, infill = infill, safeDivide = safeDivide )
        return( cls( data = xys, axes = axes, infill = infill, safeDivide = safeDivide ) )

    @staticmethod
    def defaultAxes( labelsUnits = None ) :

        return( axesModule.axes( rank = 2, labelsUnits = labelsUnits or {} ) )
