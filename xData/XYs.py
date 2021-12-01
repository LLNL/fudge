# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
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

from pqu import PQU as PQUModule

from numericalFunctions import pointwiseXY_C, pointwiseXY

from . import base as baseModule
from . import link as linkModule
from . import axes as axesModule
from . import values as valuesModule
from . import standards as standardsModule
from . import uncertainties as uncertaintiesModule

defaultAccuracy = 1e-3
xAxisIndex = 1
yAxisIndex = 0
domainEpsilon = 1e-15

def return_pointwiseXY_AsXYs( self, data, units = None, index = None, outerDomainValue = None, axes = None ) :
    """For internal use only."""

    if( not( isinstance( self, XYs1d ) ) ) : raise TypeError( 'Self must be an instance of XYs1d.' )

    if( index is None ) : index = self.index
    if( axes is None ) : axes = self.axes
    if( outerDomainValue is None ) : outerDomainValue = self.outerDomainValue
    xys1d = self.__class__( data = data, axes = axes, infill = True, safeDivide = False, index = index, outerDomainValue = outerDomainValue )
    if( ( xys1d.axes is not None ) and ( units is not None ) ) :
        for k in units : xys1d.axes[k].unit = units[k]
    return( xys1d )

def otherToSelfsUnits( self, other, checkXOnly = False ) :
    """For internal use only."""

    if( not( isinstance( other, XYs1d ) ) ) : raise TypeError( 'other instance not XYs1d instance' )

    if( ( self.axes is None ) and ( other.axes is None ) ) : return( other.nf_pointwiseXY )
    yScale = 1
    xUnitSelf = PQUModule._getPhysicalUnit( self.axes[xAxisIndex].unit )
    xScale = xUnitSelf.conversionFactorTo( other.axes[xAxisIndex].unit )
    if( not( checkXOnly ) ) :
        yUnitSelf = PQUModule._getPhysicalUnit( self.axes[yAxisIndex].unit )
        yScale = yUnitSelf.conversionFactorTo( other.axes[yAxisIndex].unit )

    other = other.nf_pointwiseXY
    if( ( xScale != 1 ) or ( yScale != 1 ) ) : other = other.scaleOffsetXAndY( xScale = xScale, yScale = yScale )

    return( other )

def getValueAsUnit( unit, quantity ) :
    """For internal use only."""

    if( not( isinstance( quantity, PQUModule.PQU ) ) ) : raise TypeError( 'Quantity is not an instance of PQUModule.PQU' )
    return( quantity.getValueAs( unit ) )

def getValueAsAxis( axis, quantity ) :

    return( getValueAsUnit( axis.unit, quantity ) )

def getOtherAndUnit( self, other ) :
    """
    For internal use only. This function is used by the multiply and divide methods. Other must be a XYs1d instance,
    a pointwiseXY instance or an object convertible to a PQU object. The first item returned ia a pointwiseXY instance 
    or number is returned.
    """

    if( isinstance( other, XYs1d ) ) :
        yUnit = other.getAxisUnitSafely( yAxisIndex )
        other = otherToSelfsUnits( self, other, checkXOnly = True )                     # Returned other is a pointwiseXY instance
    elif( isinstance( other, pointwiseXY ) ) :
        yUnit = ''
    else :                  # Must be an object convertible to a PQU instance.
        if( not( isinstance( other, PQUModule.PQU ) ) ) : other = PQUModule.PQU( other )
        yUnit = other.getUnitSymbol( )
        other = other.getValue( )
    return( other, yUnit )

def allow_XYsWithSameUnits_orNumberWithSameUnit( self, other ) :
    """
    For internal use only. Returns a pointwiseXY instance or a number.
    """

    if( isinstance( other, XYs1d ) ) :
         other = otherToSelfsUnits( self, other )                                       # Returned other is a pointwiseXY instance
    else :
        yUnit = self.getAxisUnitSafely( yAxisIndex )
        other = PQUModule.PQU( other, checkOrder = False ).getValueAs( yUnit )
    return( other )

class XYs1d( baseModule.xDataFunctional ) :

    moniker = 'XYs1d'
    dimension = 1
    mutableYUnit = True     # For __imul__ and __idiv__.

    def __init__( self, data = None, dataForm = "xys", interpolation = standardsModule.interpolation.linlinToken, axes = None,
            index = None, valueType = standardsModule.types.float64Token, outerDomainValue = None, label = None, 
            sep = ' ', initialSize = 10, overflowSize = 10, infill = True, safeDivide = False ) :
        """
        Constructor for XYs1d class. dataForm can be 'xys', 'xsandys' or 'list'.
        """

        baseModule.xDataFunctional.__init__( self, self.moniker, axes, index = index, valueType = valueType,
                outerDomainValue = outerDomainValue, label = label )

        if( not( isinstance( interpolation, str ) ) ) : raise TypeError( 'interpolation must be a string' )

        if( not( isinstance( sep, str ) ) ) : raise TypeError( 'sep must be of type str' )
        if( len( sep ) != 1 ) : raise TypeError( 'sep length must be 1 not %d' % len( sep ) )
        self.__sep = sep

        if data is None: data = []
        if( isinstance( data, XYs1d ) ) : data = data.nf_pointwiseXY

        initialSize = max( initialSize, len( data ) )
        self.__nf_pointwiseXY = pointwiseXY( data = data, dataForm = dataForm, initialSize = initialSize, overflowSize = overflowSize,
            accuracy = defaultAccuracy, biSectionMax = 16, interpolation = interpolation, infill = infill, 
            safeDivide = safeDivide )

    def __getstate__( self ) :

        state = self.__dict__.copy( )

        for key in state :
            item = getattr( self, key )
            if( isinstance( item, pointwiseXY ) ) : break

        nf_pointwiseXY_pickled = state.pop( key )
        state['nf_pointwiseXY_pickled'] = { 'name' : key, 'interpolation' : self.interpolation, 'data' : [ xy for xy in self.__nf_pointwiseXY ] }

        return( state )

    def __setstate__( self, state ) :

        nf_pointwiseXY_pickled = state.pop( 'nf_pointwiseXY_pickled' )
        self.__dict__ = state

        name = nf_pointwiseXY_pickled['name']
        interpolation = nf_pointwiseXY_pickled['interpolation']
        data = nf_pointwiseXY_pickled['data']

        nf_pointwiseXY_pickled = pointwiseXY( data = data, initialSize = len( data ), accuracy = defaultAccuracy,
                biSectionMax = 16, interpolation = interpolation, infill = True, safeDivide = False )

        setattr( self, name, nf_pointwiseXY_pickled )

    def __len__( self ) :

        return( len( self.nf_pointwiseXY ) )

    def __abs__( self ) :

        return( self.returnAsClass( self, self.nf_pointwiseXY.__abs__( ) ) )

    def __neg__( self ) :

        return( self.returnAsClass( self, self.nf_pointwiseXY.__neg__( ) ) )

    def __add__( self, other ) :

        other = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )      # Returned other is a pointwiseXY instance or number.
        return( self.returnAsClass( self, self.nf_pointwiseXY.__add__( other ) ) )

    __radd__ = __add__

    def __iadd__( self, other ) :

        other = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )      # Returned other is a pointwiseXY instance or number.
        self.nf_pointwiseXY.__iadd__( other )
        return( self )

    def __sub__( self, other ) :

        other = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )      # Returned other is a pointwiseXY instance or number.
        return( self.returnAsClass( self, self.nf_pointwiseXY.__sub__( other ) ) )

    def __rsub__( self, other ) :

        sub = self.__sub__( other )
        return( sub.__neg__( ) )

    def __isub__( self, other ) :

        other = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )      # Returned other is a pointwiseXY instance or number.
        self.nf_pointwiseXY.__isub__( other )
        return( self )

    def __mul__( self, other ) :

        unit1 = self.getAxisUnitSafely( yAxisIndex )
        other, unit2 = getOtherAndUnit( self, other )                           # Returned other is a pointwiseXY instance or number.
        unit = baseModule.processUnits( unit1, unit2, '*' )
        points = self.nf_pointwiseXY.__mul__( other )
        return( return_pointwiseXY_AsXYs( self, points, units = { yAxisIndex : unit } ) )

    __rmul__ = __mul__

    def __imul__( self, other ) :

        other, otherUnit1 = getOtherAndUnit( self, other )                      # Returned other is a pointwiseXY instance or number.
        if( not( self.mutableYUnit ) ) :
            if( otherUnit1 != '' ) : raise Exception( "Self's y-unit is immutable and other has unit of '%s'" % otherUnit1 )
        self.nf_pointwiseXY.__imul__( other )
        if( otherUnit1 != '' ) : self.axes[yAxisIndex].unit = baseModule.processUnits( self.axes[yAxisIndex].unit, otherUnit1, '*' )
        return( self )

    # division operators for Python 2.7:
    def __div__( self, other ) :

        unit1 = self.getAxisUnitSafely( yAxisIndex )
        other, unit2 = getOtherAndUnit( self, other )                           # Returned other is a pointwiseXY instance or number.
        unit = baseModule.processUnits( unit1, unit2, '/' )
        points = self.nf_pointwiseXY.__div__( other )
        return( return_pointwiseXY_AsXYs( self, points, units = { yAxisIndex : unit } ) )

    def __rdiv__( self, other ) :

        unit2 = self.getAxisUnitSafely( yAxisIndex )
        other, unit1 = getOtherAndUnit( self, other )                           # Returned other is a pointwiseXY instance or number.
        points = self.nf_pointwiseXY.__rdiv__( other )
        unit = baseModule.processUnits( unit1, unit2, '/' )
        return( return_pointwiseXY_AsXYs( self, points, units = { yAxisIndex : unit } ) )

    def __idiv__( self, other ) :

        other, otherUnit1 = getOtherAndUnit( self, other )                      # Returned other is a pointwiseXY instance or number.
        if( not( self.mutableYUnit ) ) :
            if( otherUnit1 != '' ) : raise Exception( "Self's y-unit is immutable and other has unit of '%s'" % otherUnit1 )
        self.nf_pointwiseXY.__idiv__( other )
        if( otherUnit1 != '' ) : self.axes[yAxisIndex].unit = baseModule.processUnits( self.axes[yAxisIndex].unit, otherUnit1, '/' )
        return( self )

    # division operators for Python 3.x:
    def __truediv__( self, other ) :

        unit1 = self.getAxisUnitSafely( yAxisIndex )
        other, unit2 = getOtherAndUnit( self, other )                           # Returned other is a pointwiseXY instance or number.
        unit = baseModule.processUnits( unit1, unit2, '/' )
        points = self.nf_pointwiseXY.__truediv__( other )
        return( return_pointwiseXY_AsXYs( self, points, units = { yAxisIndex : unit } ) )

    def __rtruediv__( self, other ) :

        unit2 = self.getAxisUnitSafely( yAxisIndex )
        other, unit1 = getOtherAndUnit( self, other )                           # Returned other is a pointwiseXY instance or number.
        points = self.nf_pointwiseXY.__rtruediv__( other )
        unit = baseModule.processUnits( unit1, unit2, '/' )
        return( return_pointwiseXY_AsXYs( self, points, units = { yAxisIndex : unit } ) )

    def __itruediv__( self, other ) :

        other, otherUnit1 = getOtherAndUnit( self, other )                      # Returned other is a pointwiseXY instance or number.
        if( not( self.mutableYUnit ) ) :
            if( otherUnit1 != '' ) : raise Exception( "Self's y-unit is immutable and other has unit of '%s'" % otherUnit1 )
        self.nf_pointwiseXY.__itruediv__( other )
        if( otherUnit1 != '' ) : self.axes[yAxisIndex].unit = baseModule.processUnits( self.axes[yAxisIndex].unit, otherUnit1, '/' )
        return( self )

    def __getitem__( self, indexOrSlice ) :

        if( isinstance( indexOrSlice, slice ) ) :
            start, stop, step = indexOrSlice.indices( len( self ) )
            if( step != 1 ) :  raise ValueError( "For slicing, only a step of 1 (or None) is allowed. Entered step = %d." % step )
            return( self.__getslice__( start, stop ) )
        else :
            return( self.nf_pointwiseXY.__getitem__( indexOrSlice ) )

    def __setitem__( self, indexOrSlice, xy ) :

        if( isinstance( indexOrSlice, slice ) ) :
            start, stop, step = indexOrSlice.indices( len( self ) )
            if( step != 1 ) :  raise ValueError( "For slicing, only a step of 1 (or None) is allowed. Entered step = %d." % step )
            return( self.__setslice__( start, stop, xy ) )
        else :
            if( len( xy ) != 2 ) : raise ValueError( 'right-hand-side must be list of length 2 not %s' % len( xy ) )
            self.nf_pointwiseXY.__setitem__( indexOrSlice, xy )

    def __getslice__( self, index1, index2 ) :

        return( self.returnAsClass( self, self.nf_pointwiseXY.getslice( index1, index2 ) ) )

    def __setslice__( self, index1, index2, slice1 ) :

        slice1 = otherToSelfsUnits( self, slice1 )                                  # Returned other is a pointwiseXY instanc
        self.nf_pointwiseXY.__setslice__( index1, index2, slice1 )

    @property
    def nf_pointwiseXY( self ) :
        """Returns a reference to self's hidden nf_pointwiseXY member."""

        return( self.__nf_pointwiseXY )

    @property
    def interpolation( self ) :

        return( self.nf_pointwiseXY.getInterpolation( ) )

    @interpolation.setter
    def interpolation( self, interpolation ) :

        self.nf_pointwiseXY.setInterpolation( interpolation )

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

        return( self.returnAsClass( self, self.nf_pointwiseXY.applyFunction( f, parameters, accuracy = accuracy, biSectionMax = biSectionMax, checkForRoots = checkForRoots ) ) )

    def changeInterpolation( self, interpolation, accuracy, lowerEps = 0, upperEps = 0, cls = None ) :

        if( interpolation != standardsModule.interpolation.linlinToken ) : raise ValueError( 'Only "%s" interpolation currently supported: not %s' %
                ( standardsModule.interpolation.linlinToken, interpolation ) )
        c1 = self.nf_pointwiseXY.changeInterpolation( interpolation = interpolation, accuracy = accuracy, lowerEps = lowerEps, 
                upperEps = upperEps )

        axes = self.axes
        c1 = return_pointwiseXY_AsXYs( self, c1, axes = axes, outerDomainValue = self.outerDomainValue )
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

        c1 = self.nf_pointwiseXY.cloneToInterpolation( interpolation )
        return( self.__class__.returnAsClass( self, c1, axes = self.axes, interpolation = interpolation ) )

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the form { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        if( self.axes is None ) : return
        factors = self.axes.convertUnits( unitMap )
        if( factors[:2] != [ 1., 1. ] ) : self.nf_pointwiseXY.scaleOffsetXAndY( xScale = factors[1], yScale = factors[0], insitu = True )
        self.fixValuePerUnitChange( factors )

    def clip( self, rangeMin = None, rangeMax = None ) :

        if( rangeMin is None ) : rangeMin = self.rangeMin
        if( rangeMax is None ) : rangeMax = self.rangeMax
        return( self.returnAsClass( self, self.nf_pointwiseXY.clip( rangeMin, rangeMax ) ) )

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
        factor = PQUModule.PQU( '1 ' + axis.unit ).getValueAs( newUnit )
        data = []
        xyIndex = 1 - index
        for xy in self :
            xy[xyIndex] *= factor
            data.append( xy )
        n = return_pointwiseXY_AsXYs( self, data, units = { index : newUnit } )
        return( self.returnAsClass( self, n, axes = n.axes ) )

    def copy( self ) :

        xys = self.nf_pointwiseXY.copy( )
        axes = self.axes
        if( axes is not None ) : axes = axes.copy( )
        return( self.returnAsClass( self, xys, index = self.index, outerDomainValue = self.outerDomainValue, axes = axes ) )

    __copy__ = copy

    def __deepcopy__( self, memodict = {} ) :

        copy_ = self.copy( )
        memodict[ id(self.axes) ] = copy_.axes
        return copy_

    def copyDataToNestedLists( self ) :

        return( self.nf_pointwiseXY.copyDataToXYs( ) )

    def areDomainsMutual( self, other ) :

        return( self.nf_pointwiseXY.areDomainsMutual( other.nf_pointwiseXY ) )

    def copyDataToXYs( self ) :

        return( self.__nf_pointwiseXY.copyDataToXYs( ) )

    def copyDataToXsAndYs( self ) :

        return( self.__nf_pointwiseXY.copyDataToXsAndYs( ) )

    def domain( self ) :
        """This should be deprecated."""

        return( self.nf_pointwiseXY.domain( ) )

    def dullEdges( self, lowerEps = 0., upperEps = 0., positiveXOnly = 0 ) :

        dulled = self.nf_pointwiseXY.dullEdges( lowerEps = lowerEps, upperEps = upperEps, positiveXOnly = positiveXOnly )
        return( self.returnAsClass( self, dulled ) )

    def equalProbableBins(self, numberOfBins):
        """Returns a list that that are equal probable bins of self. Currently, only supports lin-lin interpolation."""

        if self.rangeMin < 0.0: raise Exception('Cannot calculate equal probable bins for a function with negative values.')

        trimmed = self.trim()
        if len(trimmed) == 0: raise Exception('Self does not have any area - 1.')
        runningIntegrals = trimmed.runningIntegral()
        runningIntegralMax = runningIntegrals[-1]
        if runningIntegralMax == 0.0: raise Exception('Self does not have any area - 2.')

        epbs = [ trimmed[0][0] ]
        indexOfBin = 1
        integral = runningIntegralMax / numberOfBins
        absDomainMax = max(abs(trimmed[0][0]), abs(trimmed[-1][0]))
        for index, runningIntegral in enumerate(runningIntegrals):
            if indexOfBin >= numberOfBins: break
            x2, y2 = trimmed[index]
            if index > 0:
                while integral <= runningIntegral:
                    deltaArea = integral - priorRunningIntegral
                    if deltaArea == 0.0:
                        nextX = x2
                    else:
                        if y1 == y2:
                            if y1 == 0.0:
                                nextX = None
                            else:
                                nextX = x1 + deltaArea / y1
                        else:
                            slope = ( y2 - y1 ) / ( x2 - x1 )
                            sqrtArgument = y1 * y1 + 2.0 * slope * deltaArea
                            if sqrtArgument <= 0:
                                nextX = x2
                            else:
                                nextX = x1 + 2.0 * deltaArea / ( y1 + math.sqrt( sqrtArgument ) )
                    if nextX is not None: epbs.append(nextX)
                    indexOfBin += 1
                    if indexOfBin >= numberOfBins: break
                    integral = runningIntegralMax * indexOfBin / numberOfBins
            x1 = x2
            y1 = y2
            priorRunningIntegral = runningIntegral

        epbs.append(trimmed[-1][0])

# FIXME: Should look for case when one of the epbs is finite and very close to 0.0 (i.e., probably should be 0.0). Logic would look something like:
# if epbs[index] * epbs[index-1] < 0.0 # going from a negative value to a positive value.
# Now compare epbs[index] to 1e-12 * epbs[index-1] and to 1e-12 * epbs[index+1]. If smaller, set epbs[index] 0.0.
# Note, in the above epbs[index-1] may be the smaller value (in magnitude) and should be the one testing.

        return epbs

    def evaluate( self, x ) :

        return( self.__nf_pointwiseXY.evaluate( x ) )

    def lowerIndexBoundingX( self, x ) :

        return( self.__nf_pointwiseXY.lowerIndexBoundingX( x ) )

    def mutualify( self, lowerEps1, upperEps1, positiveXOnly1, other, lowerEps2, upperEps2, positiveXOnly2 ) :
        '''
        .. note:: Need to check that x units are the same.
        '''

        m1, m2 = self.nf_pointwiseXY.mutualify( lowerEps1, upperEps1, positiveXOnly1, other.nf_pointwiseXY, lowerEps2, upperEps2, positiveXOnly2 )
        return( self.returnAsClass( self, m1 ), other.returnAsClass( other, m2 ) )

    def normalize( self, insitu = False, dimension = 1 ) :
        """
        The dimension argument is ignored. Only here to be compatable with calling from XYsnd.normalize.
        """

        xys = self.nf_pointwiseXY.normalize( insitu = insitu )
        if( insitu ) : return( self )
        return( self.returnAsClass( self, xys ) )

    def setData( self, point ) :

        self.nf_pointwiseXY.setData( point )

    def setValue( self, xValue, yValue ) :

        self.nf_pointwiseXY.setValue( xValue, yValue )

    def thicken( self, sectionSubdivideMax = 1, dDomainMax = 0., fDomainMax = 1. ) :
        """
        .. note:: Need unit for dDomainMax.
        """

        t = self.nf_pointwiseXY.thicken( sectionSubdivideMax = sectionSubdivideMax, dDomainMax = dDomainMax, fDomainMax = fDomainMax )
        return( self.returnAsClass( self, t ) )

    def thin( self, accuracy ) :

        return( self.returnAsClass( self, self.nf_pointwiseXY.thin( accuracy ) ) )

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

        return( self.returnAsClass( self, self.nf_pointwiseXY.trim( ) ) )

    def union( self, other, fillWithSelf = 1, trim = 0 ) :

        other = otherToSelfsUnits( self, other, checkXOnly = True )                         # Returned other is a pointwiseXY instanc
        t = self.nf_pointwiseXY.union( other, fillWithSelf = fillWithSelf, trim = trim  )
        return( self.returnAsClass( self, t ) )

    @property
    def domainMin( self ) :

        return( self.nf_pointwiseXY.domainMin( ) )

    @property
    def domainMax( self ) :

        return( self.nf_pointwiseXY.domainMax( ) )

    @property
    def domainUnit( self ) :

        return( self.getAxisUnitSafely( xAxisIndex ) )

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    @property
    def domainGrid( self ) :

        return( self.nf_pointwiseXY.domainGrid( 1.0 ) )

    @property
    def rangeMin( self ) :

        return( self.nf_pointwiseXY.rangeMin( ) )

    @property
    def rangeMax( self ) :

        return( self.nf_pointwiseXY.rangeMax( ) )

    @property
    def rangeUnit( self ) :

        return( self.getAxisUnitSafely( yAxisIndex ) )

    def rangeUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.rangeUnit ).getValueAs( unitTo ) )

    def domainSlice( self, domainMin = None, domainMax = None, fill = 1, dullEps = 0. ) :
        """
        Returns a new instance with self sliced between ``domainMin`` and ``domainMax``.

        :param domainMin:   [optional] the lower x-value of the slice, default is domain minimum of self,
        :param domainMax:   [optional] the upper x-value of the slice, default is domain maximum of self,
        :param fill:        [optional] if True, points are added at domainMin and domainMax if they are not in self, 
                                       else only existing points in the range [domainMin, domainMax] are included.
        :param dullEps:     [optional] (Currently not implemented) the lower and upper points are dulled, default is 0.
        """

        if( domainMin is None ) : domainMin = self.domainMin
        if( domainMax is None ) : domainMax = self.domainMax
        s = self.nf_pointwiseXY.domainSlice( domainMin = domainMin, domainMax = domainMax, fill = fill, dullEps = dullEps )
        return( self.returnAsClass( self, s ) )

    def __mod__( self, other ) : raise NotImplementedError( 'Currently, mod is not implemented' )

    def __pow__( self, other ) : raise NotImplementedError( 'Currently, pow is not implemented' )

    def __exp__( self, other ) : raise NotImplementedError( 'Currently, __exp__ is not implemented' )

    def convolute( self, func ) : 
        """
        Uses pointwiseXY's convolute function to do the grunt work.  
        """

        if( isinstance( func, XYs1d ) ) : func = func.nf_pointwiseXY
        if( not( isinstance( func, pointwiseXY ) ) ) : raise TypeError( 'func argument of convolute() must be an instance of pointwiseXY' )
        return( self.returnAsClass( self, self.nf_pointwiseXY.convolute( func, 0 ) ) )

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

        yUnit = PQUModule.PQU( 1, self.getAxisUnitSafely( yAxisIndex ) )
        if( f2 is None ) :
            if( f3 is None ) : 
                grouped = self.groupOneFunction( xs, norm = norm )
            else :
                grouped = self.groupTwoFunctions( xs, f3, norm = norm )
                yUnit = yUnit * PQUModule.PQU( 1,  f3.getAxisUnitSafely( yAxisIndex ) )
        else :
            yUnit = yUnit * PQUModule.PQU( 1, f2.getAxisUnitSafely( yAxisIndex ) )
            if( f3 is None ) :
                grouped = self.groupTwoFunctions( xs, f2, norm = norm )
            else :
                grouped = self.groupThreeFunctions( xs, f2, f3, norm = norm )
                yUnit = yUnit * PQUModule.PQU( 1, f3.getAxisUnitSafely( yAxisIndex ) )
        if( norm is None ) :
            yUnit = PQUModule.PQU( 1, self.getAxisUnitSafely( xAxisIndex ) ) * yUnit
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

        if( type( xs ) == list ) :
            boundaries = xs
        elif( type( xs.values ) == list ) :
            boundaries = xs.values
        else :
            boundaries = xs.values.values

        if( isinstance( norm, XYs1d ) ) : norm = norm.nf_pointwiseXY

        return( self.nf_pointwiseXY.groupOneFunction( boundaries, norm = norm ) )

    def groupTwoFunctions( self, xs, f2, norm = None ) :
        '''.. note:: Need unit of xs.'''

        if( type( xs ) == list ) :
            boundaries = xs
        elif( type( xs.values ) == list ) :
            boundaries = xs.values
        else :
            boundaries = xs.values.values

        if( isinstance( f2, XYs1d ) ) : f2 = f2.nf_pointwiseXY
        if( isinstance( norm, XYs1d ) ) : norm = norm.nf_pointwiseXY

        return( self.nf_pointwiseXY.groupTwoFunctions( boundaries, f2, norm = norm ) )

    def groupThreeFunctions( self, xs, f2, f3, norm = None ) :
        '''.. note:: Need unit of xs.'''

        if( type( xs ) == list ) :
            boundaries = xs
        elif( type( xs.values ) == list ) :
            boundaries = xs.values
        else :
            boundaries = xs.values.values

        if( isinstance( f2, XYs1d ) ) : f2 = f2.nf_pointwiseXY
        if( isinstance( f3, XYs1d ) ) : f3 = f3.nf_pointwiseXY
        if( isinstance( norm, XYs1d ) ) : norm = norm.nf_pointwiseXY

        return( self.nf_pointwiseXY.groupThreeFunctions(  boundaries, f2, f3, norm = norm ) )

    def integrate( self, domainMin = None, domainMax = None ) :
        """
        Definite integral of current ``XYs1d`` instance from ``domainMin`` to ``domainMax``:
        
        .. math::
            \int_{domainMin}^{domainMax} dx \; XYs(x)

        If ``domainMin`` or ``domainMax`` is unspecified, it is taken from the domain of the self.
        """

        if( len( self ) == 0 ) : return( 0.0 )
        unit = self.getAxisUnitSafely( xAxisIndex )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin )
        domainMax = min( domainMax, self.domainMax )
        return( PQUModule.PQU( self.nf_pointwiseXY.integrate( domainMin = domainMin, domainMax = domainMax ), 
                baseModule.processUnits( unit, self.getAxisUnitSafely( yAxisIndex ), '*' ), checkOrder = False ) )

    def integrateWithWeight_x( self, domainMin = None, domainMax = None ) :

        if( len( self ) == 0 ) : return( 0.0 )
        unit = self.getAxisUnitSafely( xAxisIndex )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin )
        domainMax = min( domainMax, self.domainMax )
        return( self.__nf_pointwiseXY.integrateWithWeight_x( domainMin, domainMax ) )

    def integrateWithWeight_sqrt_x( self, domainMin = None, domainMax = None ) :

        if( len( self ) == 0 ) : return( 0.0 )
        unit = self.getAxisUnitSafely( xAxisIndex )
        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin )
        domainMax = min( domainMax, self.domainMax )
        return( self.__nf_pointwiseXY.integrateWithWeight_sqrt_x( domainMin, domainMax ) )

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
            myData.append( [ domainMax, myData[-1][1] + self.nf_pointwiseXY.integrate( domainMin = domainMin, domainMax = domainMax ) ] )
        return XYs1d( myData, axes = myAxes )

    def integrateTwoFunctions( self, f2, domainMin = None, domainMax = None ) :
        """

        :param f2:
        :param domainMin:
        :param domainMax:
        :return:
        """

        if( not isinstance( f2, XYs1d ) ) : raise TypeError( "f2 must be an instance of an XYs1d" )

        unit = ''
        integrationUnit = ''
        if( ( self.axes is not None ) and ( f2.axes is not None ) ) :
            unit = self.axes[xAxisIndex].unit
            if( f2.axes[xAxisIndex].unit != unit ) :
                f2 = f2.copy( )
                f2.convertAxisToUnit( xAxisIndex, unit )
            integrationUnit = baseModule.processUnits( baseModule.processUnits( unit, self.axes[yAxisIndex].unit, '*' ), f2.axes[yAxisIndex].unit, '*' )

        domainMin, domainMax = baseModule.getDomainLimits( self, domainMin, domainMax, unit )
        domainMin = max( domainMin, self.domainMin, f2.domainMin )
        domainMax = min( domainMax, self.domainMax, f2.domainMax )

        return( PQUModule.PQU( self.groupTwoFunctions( [ domainMin, domainMax ], f2.nf_pointwiseXY )[0], integrationUnit, checkOrder = False ) )

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

        unit = ''
        if( ( self.axes is not None ) and ( f2.axes is not None ) and ( f3.axes is not None ) ) :
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

        return( self.groupThreeFunctions( [ domainMin, domainMax ], f2, f3 )[0] )

    def inverse( self ) :

        return( self.returnAsClass( self, self.nf_pointwiseXY.inverse( ) ) )

    def pdfOfY( self, epsilon, domainMin = None, domainMax = None ) :
        """     
        This method calculates the pdf of the y-values of self (i.e., pdf(y)).
        """

        def addData( pdf, data, epsilon ) :

            if( pdf is None ) : return( XYs1d( data = data ) )
            data = XYs1d( data = data )
            pdf, data = pdf.mutualify( epsilon, epsilon, False, data, epsilon, epsilon, False )
            pdf += data
            return( pdf )
                
        region = self
        
        if( domainMin is None ) : domainMin = region.domainMin
        if( domainMax is None ) : domainMax = region.domainMax
            
        if( len( region ) < 2 ) : raise ValueError( 'self must have at least 2 points' )
        if( domainMin < region.domainMin ) :
            if( region[0][1] != 0.0 ) : raise ValueError( "Specified domainMin is less than self's domainMin." )
            region[0] = [ domainMin, 0.0 ]
        if( domainMax > region.domainMax ) : 
            if( region[-1][1] != 0.0 ) : raise ValueError( "Specified domainMax is greater than self's domainMax." )
            region[len(self)] = [ domainMax, 0.0 ]

        region = region.domainSlice( domainMin, domainMax, fill = True )

        domainWidth = domainMax - domainMin
        deltas = {} 
        pdf = None

        if( region.interpolation == standardsModule.interpolation.flatToken ) :
            x1 = None
            for x2, y2 in region :
                if( x1 is not None ) :
                    if( y1 not in deltas ) : deltas[y1] = 0.0
                    deltas[y1] += ( x2 - x1 ) / domainWidth
                x1 = x2
                y1 = y2
        elif( region.interpolation == standardsModule.interpolation.linlinToken ) :
            slope = 1
            data = []
            for i1, ( x2, y2 ) in enumerate( region ) :
                if( i1 > 0 ) :
                    if( y1 == y2 ) :
                        if( y1 not in deltas ) : deltas[y1] = 0.0
                        deltas[y1] += ( x2 - x1 ) / domainWidth
                    else :
                        dy = y2 - y1
                        dXdY = ( x2 - x1 ) / ( abs( dy ) * domainWidth )
                        if( i1 == 1 ) :
                            data = [ [ y1, dXdY ] ]
                            slope = dy
                        if( dy * slope < 0 ) :
                            pdf = addData( pdf, data, epsilon )
                            data = [ [ y1, dXdY ] ]
                        if( dy < 0 ) :
                            data.insert( 0, [ y2, dXdY ] )
                        else :
                            data.append( [ y2, dXdY ] )
                        slope = dy
                x1 = x2
                y1 = y2
            pdf = addData( pdf, data, epsilon )
        else :      
            raise Exception( 'Unsupported interpolation = "%s"' % self.interpolation )

        for yValue in sorted( deltas ) :
            if( yValue == 0.0 ) :
                data = [ [ 0.0, 1.0 ], [ epsilon, 0.0 ] ]
            else :
                delta = [ [ yValue * ( 1.0 - epsilon ), 0.0 ], [ yValue, 1.0 ], [ yValue * ( 1.0 + epsilon ), 0.0 ] ]
            delta = XYs1d( data ).normalize( )
            delta *= deltas[yValue]
            pdf = addData( pdf, delta, epsilon )

        pdf.normalize( insitu = True )

        return( pdf )

    def runningIntegral( self ) :

        return( self.__nf_pointwiseXY.runningIntegral( ) )

    def scaleDependent( self, value, insitu = False ) :

        xys = self
        if( not( insitu ) ) : xys = self.copy( )
        xys.nf_pointwiseXY.scaleOffsetXAndY( yScale = value, insitu = True )
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

    def toLinearXYsClass( self ) :

        return( XYs1d )

    def toPointwiseLinear( self, **kwargs ) :
        """
        Returns a new instance, converted to lin-lin interpolation with added points to maintain desired accuracy.
        For more details see toPointwise_withLinearXYs.
        """

        return( self.toPointwise_withLinearXYs( **kwargs ) )

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

    def toString( self, pairsPerLine = 1, format = " %16.8e %16.8e", pairSeparator = "" ) :

        return( self.nf_pointwiseXY.toString( pairsPerLine = pairsPerLine, format = format, pairSeparator = pairSeparator ) )

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
        for x, y in self.nf_pointwiseXY.copyDataToXYs( ) :
            xys.append( x )
            xys.append( y )
        XMLList += valuesModule.values( xys, valueType = self.valueType, sep = self.__sep ).toXMLList( indent2, **kwargs )
        if( self.uncertainty ) : XMLList += self.uncertainty.toXMLList( indent = indent2, **kwargs )
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

    def plot( self, xylog = 0, domainMin = None, domainMax = None, rangeMin = None , rangeMax = None, title = '',
              multiPlot = False, dataKey = None ) :

        from .interactivePlot import plotbase as plotbaseModule
        from .interactivePlot import multiplot as interactivePlotModule

        def getUnitlessNumber( value, unit, default ) :

            if( value is None ) : return( default )
            if( not( isinstance( value, PQUModule.PQU ) ) ) : value = PQUModule.PQU( value )
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

        plotOptions = plotbaseModule.parsePlotOptions( domainMin, domainMax, rangeMin, rangeMax, xLabel, yLabel, title )
        if dataKey is None:
            if self.label:
                dataKey = self.label
            elif not multiPlot:
                dataKey = 'noLabel' if title == '' else title

        plotData = {dataKey: self.copyDataToXsAndYs()}

        if multiPlot:
            return plotOptions, plotData

        else:
            interactivePlotModule.MultiPlotWithPyQt5(plotOptions, plotData)

    def ysMappedToXs( self, cls, grid, label = None, extendToEnd = False ) :

        offset, Ys = self.nf_pointwiseXY.ysMappedToXs( grid.values )
        if( extendToEnd ) : Ys += ( len( grid.values ) - ( offset + len( Ys ) ) ) * [ 0.0 ]
        Ys1d = cls( valuesModule.values( Ys, start = offset ), axes = self.axes, label = label )
        Ys1d.axes[1] = linkModule.link2( grid.moniker, instance = grid, keyName = 'index', keyValue = 1 )
        return( Ys1d )

    @classmethod
    def returnAsClass( cls, self, other, index = None, outerDomainValue = None, axes = None, interpolation = None ) :
        """
        Returns other as a class of cls. Other must be a sub-class pointwiseXY. cls must be a sub-class of XYs1d. 
        If index and axes are not specified, they are taken from self. The main use of this classmethod is for 
        methods like __add__ where the addends may be a class derived from XYs1d. For example, the crossSection.XYs1d
        class is derived from The XYs1d class.  If two instances xSec1 and xSec2 of the crossSection.pointwise class 
        are added together (i.e., xSec1 + xSec2) then, since the __add__ method used this classmethod, the returned 
        instance will also be an instance of crossSection.pointwise class.
        """

        if( not( isinstance( self, XYs1d ) ) ) : raise TypeError( 'Self is not a XYs1d instance.' )
        if( isinstance( other, XYs1d ) ) : other = other.nf_pointwiseXY

        if( index is None ) : index = self.index
        if( outerDomainValue is None ) : outerDomainValue = self.outerDomainValue
        if( axes is None ) : axes = self.axes
        if( interpolation is None ) :
            interpolation = self.interpolation
            if( len( self ) == 0 ) : interpolation = other.getInterpolation( )

        xys = cls( data = other, interpolation = interpolation, axes = axes, overflowSize = 10, infill = other.getInfill( ), 
                safeDivide = other.getSafeDivide( ), index = index, outerDomainValue = outerDomainValue )

        return( xys )

    @classmethod
    def parseXMLNode( cls, xDataElement, xPath, linkData, axes = None, **kwargs ) :
        """
        Translate an XYs1d XML element into the python XYs1d xData class.
        """

        xmlAttr = False
        for attrName in ( 'label', 'outerDomainValue' ) :
            if xDataElement.get(attrName) is not None:
                xmlAttr = True
                xPath.append( '%s[@%s="%s"]' % (xDataElement.tag, attrName, xDataElement.get(attrName) ) )
                break
        if( not xmlAttr ) : xPath.append( xDataElement.tag )

        attrs = { 'interpolation' : standardsModule.interpolation.linlinToken, 'label' : None, 
                'index' : None, 'outerDomainValue' : None }
        attributes = { 'interpolation' : str, 'label' : str, 'index' : int, 'outerDomainValue' : float }
        for key, item in list( xDataElement.items( ) ) :
            if( key not in attributes ) : raise TypeError( 'Invalid attribute "%s"' % key )
            attrs[key] = attributes[key]( item )

        values, uncertainty = None, None
        for subElement in xDataElement :
            if( subElement.tag == 'axes' ) :
                axes = axesModule.axes.parseXMLNode( subElement, xPath, linkData )
            elif( subElement.tag == 'values' ) :
                values = valuesModule.values.parseXMLNode( subElement, xPath, linkData )
            elif( subElement.tag == 'uncertainty' ) :
                uncertainty = uncertaintiesModule.uncertainty.parseXMLNode( subElement, xPath, linkData )
            else :
                raise TypeError( 'sub-element "%s" not valid' % subElement.tag )

        if( values is None ) : raise Exception( 'values element missing' )
        attrs['sep'] = values.sep

        xys = cls( data = values, dataForm = "list", axes = axes, **attrs )
        if uncertainty is not None: xys.uncertainty = uncertainty
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

    @staticmethod
    def multiPlot( curve1ds, **kwargs ) :
        """
        Plots a list of 1d curves on the same plot. Uses each curve's 'plotLegendKey' as the legend key. Each curve must be
        and XYs1d instance.
        """

        def argumentValue( name, default ) :

            kwargsValue = kwargs.get( name, None )
            if( kwargsValue is not None ) : return( kwargsValue )
            return( default )

        from .interactivePlot import plotbase as plotbaseModule
        from .interactivePlot import multiplot as multiplotModule

        if( len( curve1ds ) == 0 ) : return

        options = { 'xylog' : 0, 'title' : None, 'title' : None, 'xLabel' : None, 'yLabel' : None, 
                    'domainMin' : None, 'domainMax' : None, 'rangeMin' : None, 'rangeMax' : None }

        plotData = {}
        domainMin = []
        domainMax = []
        rangeMin = []
        rangeMax = []
        for index, curve1d in enumerate( curve1ds ) :
            if( hasattr( curve1d, 'plotLegendKey' ) ) :
                plotLegendKey = curve1d.plotLegendKey
            else :
                plotLegendKey = 'curve %d' % index
            if( not( isinstance( curve1d, XYs1d ) ) ) : raise TypeError( 'Curve is not an XYs1d instance.' )
            plotData[plotLegendKey] = curve1d.copyDataToXsAndYs( )
            domainMin.append( curve1d.domainMin )
            domainMax.append( curve1d.domainMax )
            rangeMin.append( curve1d.rangeMin )
            rangeMax.append( curve1d.rangeMax )
        domainMin = argumentValue( 'domainMin', min( domainMin ) )
        domainMax = argumentValue( 'domainMax', max( domainMax ) )
        rangeMin =  argumentValue( 'rangeMin',  min( rangeMin ) )
        rangeMax =  argumentValue( 'rangeMax',  max( rangeMax ) )

        title = kwargs.get( 'title', 'No title' )
        xLabel = kwargs.get( 'xLabel', 'x' )
        yLabel = kwargs.get( 'yLabel', 'y' )

        plotOptions = plotbaseModule.parsePlotOptions( domainMin, domainMax, rangeMin, rangeMax, xLabel, yLabel, title )
        plotOptions['logs'] = kwargs.get( 'logs', 0 )

        multiplotModule.MultiPlotWithPyQt5( plotOptions, plotData )

    def checkEqualProbableBinsResult(self, epbs, epsilon=1e-12, printResults=False):
        """This method checks the results of XYs1d.equalProbableBins and executes a raise if it finds an issue.

            :param self:            The XYs1d instance that equalProbableBins was call for.
            :param epbs:            The results returned by equalProbableBins.
            :param epsilon:         The relative accuracy to compare the values in epbs to.
            :param printResults:    If true, prints information on each value in epbs.

            :return:                The number of errors found.
        """

        errCounter = 0

        trimmed = self.trim()
        runningIntegrals = trimmed.runningIntegral()
        try:
            inverse = XYs1d( data = [ [ x, runningIntegrals[index] ] for index, (x, y) in enumerate(trimmed) ] ).inverse()
        except:
            print( runningIntegrals )
            print( self.toString())
            raise

        numberOfBins = len(epbs) - 1
        for index in range(numberOfBins+1):
            runningIntegral = index * runningIntegrals[-1] / numberOfBins
            integral = float(self.integrate( domainMax = epbs[index] ))
            delta = runningIntegral - integral
            ratio = 0.0
            if( delta != 0 ):
                ratio = delta / runningIntegral 
                if abs(ratio) > epsilon:
                    errCounter += 1
            if printResults:
                s = ''
                if abs(ratio) > epsilon: s = ' *********'
                print( '%4d %20.13e %10.10e %18.10e %10.2e, %10.2e' % (index, epbs[index], runningIntegral, integral, delta, ratio ), s )

        return errCounter
