# <<BEGIN-copyright>>
# <<END-copyright>>

"""
    This module contains the class ``XYs``. This class treats a list of :math:`(x_i, y_i)` pairs as if they were the function :math:`y(x)`.
    That is, it is a numerical representation of :math:`f(x)`. As an example, let 
    
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
        :axes:            Description of the x and y data attributes (e.g., label, interpolation, units).

    A ``XYs`` object can be the ``XY`` instances of higher dimensional math objects like ``W_XYs``. In this
    case, the ``XYs`` instances are secondary instances.
"""

import axes
from fudge.core.ancestry import ancestry
from pqu import physicalQuantityWithUncertainty
from fudge.core.utilities import subprocessing
from fudge.core.math import fudgemath
from pqu.physicalQuantityWithUncertainty import PhysicalQuantityWithUncertainty

useExtension = True
try:
    from numericalFunctions import pointwiseXY_C, pointwiseXY
except ImportError:
    # c extensions are missing. For better performance, please build extensions
    useExtension = False
    from fudge.core.math.pointwiseXY import pointwiseXY
from fudge.core.utilities import brb

monikerXYs = 'XYs'

__metaclass__ = type

def return_pointwiseXY_AsXYs( self, C, units = {}, parent = None, index = None, value = None, axes = None, moniker = monikerXYs, isPrimaryXData = None,
        template = None, changeInterpolationSubFunction = None ) :

    if( index is None ) : index = self.index
    if( axes is None ) : axes = self.axes
    if( useExtension ) :
        t = C
        if( template is not None ) : t = template
        c = XYs( axes, C, t.getAccuracy( ), biSectionMax = t.getBiSectionMax( ), infill = t.getInfill( ), safeDivide = t.getSafeDivide( ), \
            parent = parent, index = index, value = value, moniker = moniker, isPrimaryXData = isPrimaryXData, changeInterpolationSubFunction = changeInterpolationSubFunction )
    else :
        c = XYs( axes, C, 0.001, biSectionMax = 4, infill = 1, safeDivide = 0, parent = parent, index = index, value = value,
                moniker = moniker, isPrimaryXData = isPrimaryXData, changeInterpolationSubFunction = changeInterpolationSubFunction )
    for k in units : c.axes[k].unit = units[k]
    return( c )

def raiseNotSameUnit( v1, v2 ) :

    if( v1.getUnit( ) != v2.getUnit( ) ) : raise ValueError( 'unit = %s != unit = %s' % ( v1.getUnit( ), v2.getUnit( ) ) )

def raiseNotSameUnits( self, other, checkXOnly = False ) :

    if( not( isinstance( other, XYs ) ) ) : raise TypeError( 'other instance of %s and not XYs' % brb.getType( other ) )
    if( len( self.axes ) != len( other.axes ) ) : raise Exception( 'axes not the same' )
    raiseNotSameUnit( self.axes[0], other.axes[0] )
    if( not( checkXOnly ) ) : raiseNotSameUnit( self.axes[1], other.axes[1] )

def getValueAsUnit( unit, quantity ) :

    if( unit == '' ) :         # Must handle special case of dimensionless unit until PhysicalQuantityWithUncertainty can handle it.
        if( fudgemath.isNumber( quantity ) ) : return( quantity )
    if( not( isinstance( quantity, physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty ) ) ) : 
        raise TypeError( 'quantity is not an instance of physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty' )
    return( quantity.getValueAs( unit ) )

def getValueAsAxis( axis, quantity ) :

    return( getValueAsUnit( axis.getUnit( ), quantity ) )

def processUnits( unit1, unit2, operator ) :

    if( operator not in [ '*', '/' ] ) : raise ArithmeticError( 'unsupported unit operation "%s"' % operator )
    if( unit1.strip( ) == '' ) :
        if( unit2.strip( ) == '' ) : return( '' )
        if( operator == '*' ) : return( unit2 )
        if( operator == '/' ) : 
            pq = 1 / physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( 1, unit2 )
            return( pq.getUnitSymbol( ) )
    elif( unit2.strip( ) == '' ) :
        return( unit1 )
    pq1 = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( 1, unit1 )
    pq2 = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( 1, unit2 )
    v = eval( 'pq1 %s pq2' % operator )
    if( fudgemath.isNumber( v ) ) : return( '' )
    return( v.getUnitSymbol( ) )

def getOtherAndUnit( self, other ) :

    if( fudgemath.isNumber( other ) ) :
        unit2 = ''
    elif( isinstance( other, physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty ) ) :
        unit2 = other.getUnitName( )
        other = other.getValue( )
    elif( isinstance( other, XYs ) ) :
        raiseNotSameUnit( self.axes[0], other.axes[0] )
        unit2 = other.axes[1].getUnit( )
    else :
        raise TypeError( 'operation of XYs with object of type %s is not supported' % brb.getType( other ) )
    return( other, unit2 )

def allow_XYsWithSameUnits_orNumberWithSameUnit( self, other ) :

    unit = self.axes[1].getUnit( )
    if( type( other ) == type( '' ) ) :     # If its is a string try to convert to a number or PhysicalQuantityWithUncertainty.
        try :
            other = float( other )
        except :
            other = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( other )
    if( fudgemath.isNumber( other ) ) :
        if( unit != '' ) : raise ValueError( 'unitless number found where unit = "%s" needed' % unit )
    elif( isinstance( other, physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty ) ) :
        if( unit == '' ) : raise ValueError( 'unit number = "%s" found where unitless needed' % other )
        other = other.getValueAs( unit )
    else :
         raiseNotSameUnits( self, other )
    return( other )

def evaluateValueAsUnit( unit, value ) :

    if( unit == '' ) :
        if( not( fudgemath.isNumber( value ) ) ) : raise Exception( 'instance is type %s and not a number: "%s"' % ( brb.getType( value ), value ) )
        return( value )
    if( fudgemath.isNumber( value ) ) : raise ValueError( 'instance is a number and not of unit "%s"' % unit )
    if( type( value ) == type( '' ) ) : value = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( value )
    if( not( isinstance( value, physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty ) ) ) : 
        raise TypeError( 'instance is type %s and not a physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty' % brb.getType( value ) )
    return( value.getValueAs( unit ) )

def initializeFromTemplate( self, template ) :
    """This is a routine to copy all the data from template to self. This should only be called once for self."""

    XYs.__init__( self, template.axes, template, template.getAccuracy( ), biSectionMax = template.getBiSectionMax( ), \
        infill = template.getInfill( ), safeDivide = template.getSafeDivide( ), index = template.index, value = template.value, \
        parent = template.parent )

def compareAxisToParent( axis, parent_, parentAxisIndex ) :

    parent = parent_
    while( parent.parent != None ) :
        if( not( hasattr( parent.parent, 'axes' ) ) ) : break
        parent = parent.parent
    parentAxis = parent.axes[parentAxisIndex]
    if( isinstance( axis, axes.interpolationAxis ) and ( isinstance( parentAxis, axes.axis ) ) ) : return
    if( axis.getLabel( ) != parentAxis.getLabel( ) ) : 
        raise ValueError( "axis' label not the same: '%s' vs. '%s'" % ( axis.getLabel( ), parentAxis.getLabel( ) ) )
    if( axis.getUnit( ) != parentAxis.getUnit( ) ) : raise ValueError( "axis' units not the same: '%s' vs. '%s'" % ( axis.getUnit( ), parentAxis.getUnit( ) ) )
    if( axis.frame != parentAxis.frame ) : raise ValueError( "axis' frames not the same: '%s' vs. '%s'" % ( axis.frame, parentAxis.frame ) )

class XYs( pointwiseXY, ancestry ) :

    xData = monikerXYs

    def __init__( self, axes_, data, accuracy, dataForm = "xys", initialSize = None, overflowSize = 10, biSectionMax = 6, 
            infill = True, safeDivide = False, index = None, value = None, parent = None, moniker = xData, isPrimaryXData = None, 
            changeInterpolationSubFunction = None ) :
        """
        Constructor for XYs class. dataForm can be 'xys', 'xsandys' or list.
        """

        attribute = None
        if( parent is not None ) :
            attribute = 'index'
#            compareAxisToParent( axes_[0], parent, -2 )        # This needs work?????????
#            compareAxisToParent( axes_[1], parent, -1 )
            if( index is None ) : raise ValueError( 'index argument cannot be None when parent is not None' )
        ancestry.__init__( self, moniker, parent, attribute = attribute )

        axes.isValidAxes( axes_ )
        axes_.checkDimension( 2 )
        self.changeInterpolationSubFunction = changeInterpolationSubFunction
        interpolation = str( axes_[0].interpolation )
        if( axes.chargedParticleToken in interpolation ) :
            interpolation = '%s,%s' % ( axes.linearToken, axes.otherToken )
            if( changeInterpolationSubFunction is None ) : raise ValueError( '"changeInterpolationSubFunction" argument needed for "%s" interpolation' % axes.chargedParticleToken )
        if( initialSize is None ) : initialSize = len( data )

        pointwiseXY.__init__( self, initialSize = initialSize, overflowSize = overflowSize, data = data, dataForm = dataForm,
            accuracy = accuracy, biSectionMax = biSectionMax, interpolation = interpolation, infill = infill, safeDivide = safeDivide )

        self.axes = axes_.copy( )           # Even secondary must have axes data
        if( isinstance( self.axes, axes.referenceAxes ) ) :
            self.axes.setParent( axes_.parent )
        else :
            self.axes.setParent( self )
        self.index = index
        self.value = value
        if( isPrimaryXData is not None ) : self.isPrimaryXData = isPrimaryXData

    def __abs__( self ) :

        return( self.returnAsClass( self, pointwiseXY.__abs__( self ) ) )

    def __neg__( self ) :

        return( self.returnAsClass( self, pointwiseXY.__neg__( self ) ) )

    def __add__( self, other ) :

        other_ = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )
        return( self.returnAsClass( self, pointwiseXY.__add__( self, other_ ) ) )

    __radd__ = __add__

    def __sub__( self, other ) :

        other_ = allow_XYsWithSameUnits_orNumberWithSameUnit( self, other )
        return( self.returnAsClass( self, pointwiseXY.__sub__( self, other_ ) ) )

    def __rsub__( self, other ) :

        sub = self.__sub__( other )
        return( sub.__neg__( ) )

    def __mul__( self, other ) :

        unit1 = self.axes[1].getUnit( )
        other, unit2 = getOtherAndUnit( self, other )
        unit = processUnits( unit1, unit2, '*' )
        m = pointwiseXY.__mul__( self, other )
        return( return_pointwiseXY_AsXYs( self, m, units = { 1 : unit } ) )

    __rmul__ = __mul__

    def __div__( self, other ) :

        unit1 = self.axes[1].getUnit( )
        other, unit2 = getOtherAndUnit( self, other )
        unit = processUnits( unit1, unit2, '/' )
        d = pointwiseXY.__div__( self, other )
        return( return_pointwiseXY_AsXYs( self, d, units = { 1 : unit } ) )

    def __rdiv__( self, other ) :

        unit2 = self.axes[1].getUnit( )
        other, unit1 = getOtherAndUnit( self, other )
        d = pointwiseXY.__rdiv__( self, other )
        unit = processUnits( unit1, unit2, '/' )
        return( return_pointwiseXY_AsXYs( self, d, units = { 1 : unit } ) )

    def getitem_units( self, index ) :

        x, y = pointwiseXY.__getitem__( self, index )
        if( self.axes[0].getUnit( ) != '' ) : x = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( x, self.axes[0].getUnit( ) )
        if( self.axes[1].getUnit( ) != '' ) : y = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( y, self.axes[1].getUnit( ) )
        return( [ x, y ] ) 

    def __setitem__( self, index, xy ) :

        if( len( xy ) != 2 ) : raise ValueError( 'right-hand-side must be list of length 2 not %s' % len( xy ) )
        pointwiseXY.__setitem__( self, index, xy )

    def setitem_units( self, index, xy ) :

        if( len( xy ) != 2 ) : raise ValueError( 'right-hand-side must be list of length 2 not %s' % len( xy ) )
        xy = [ getValueAsAxis( self.axes[0], xy[0] ), getValueAsAxis( self.axes[1], xy[1] ) ]
        pointwiseXY.__setitem__( self, index, xy )

    def __getslice__( self, index1, index2 ) :

        return( self.returnAsClass( self, pointwiseXY.__getslice__( self, index1, index2 ) ) )

    def __setslice__( self, index1, index2, slice ) :

        raiseNotSameUnits( self, slice )
        pointwiseXY.__setslice__( self, index1, index2, slice )

    def applyFunction( self, f, parameters, accuracy = -1, biSectionMax = -1, checkForRoots = False ) :
        """applyFunction maps the function 'f' onto each y-value in an XYs object, returning a new XYs object.
        Additional points may be added to preserve the desired accuracy.

        For example, the following will transform all negative y-values in an XYs object to zero:
        >>> newXYs = XYs.applyFunction( lambda y,tmp: 0 if y<0 else y, None )

        Extra parameters to the applied function can be passed in the 'parameters' argument.
        See XYs.XYs documentation for details on 'accuracy' and 'biSectionMax'"""

        return( self.returnAsClass( self, pointwiseXY.applyFunction( self, f, parameters, accuracy = accuracy, biSectionMax = biSectionMax, checkForRoots = checkForRoots ) ) )

    def changeInterpolation( self, xInterpolation, yInterpolation, accuracy = None, lowerEps = 0, upperEps = 0, cls = None ) :

        if( xInterpolation != axes.linearToken ) : raise ValueError( 'Only linear interpolation currently supported for x-axis, not %s' % xInterpolation )
        if( yInterpolation != axes.linearToken ) : raise ValueError( 'Only linear interpolation currently supported for y-axis, not %s' % yInterpolation )
        if( accuracy is None ) : accuracy = self.getAccuracy( )
        if( self.changeInterpolationSubFunction is None ) :
            c = pointwiseXY.changeInterpolation( self, interpolation = '%s,%s' % ( xInterpolation, yInterpolation ), 
                accuracy = accuracy, lowerEps = lowerEps, upperEps = upperEps )
        else :
            c = self.changeInterpolationSubFunction( self, xInterpolation, yInterpolation, accuracy, lowerEps = lowerEps, upperEps = upperEps )
        axes_ = self.axes.copy( standAlone = True )
        axes_[0].setInterpolation( axes.interpolationXY( xInterpolation, yInterpolation ) )
        c = return_pointwiseXY_AsXYs( self, c, axes = axes_, template = self, value = self.value )
        if( cls is None  ) : cls = self
        return( cls.returnAsClass( self, c, axes_ = axes_, moniker = self.moniker ) )

    def changeInterpolationIfNeeded( self, allowedInterpolations, accuracy = None, lowerEps = 0, upperEps = 0, cls = None ) :

        independent, dependent, qualifier = self.axes[0].interpolation.getInterpolationTokens( )
        for interpolation in allowedInterpolations :
            if( [ independent, dependent ] == interpolation ) : return( self )
        return( self.changeInterpolation( allowedInterpolations[0][0], allowedInterpolations[0][1], accuracy = accuracy, lowerEps = lowerEps, 
            upperEps = upperEps, cls = cls ) )

    def clip( self, yMin = None, yMax = None ) :

        unit = self.axes[1].getUnit( )
        if( yMin is None ) : yMin = self.yMin( )
        if( yMax is None ) : yMax = self.yMax( )
        return( self.returnAsClass( self, pointwiseXY.clip( self, yMin, yMax ) ) )

    def clip_units( self, yMin = None, yMax = None ) :

        unit = self.axes[1].getUnit( )
        if( yMin is None ) : yMin = self.yMin( )
        yMin = evaluateValueAsUnit( unit, yMin )
        if( yMax is None ) : yMax = self.yMax( )
        yMax = evaluateValueAsUnit( unit, yMax )
        c = pointwiseXY.clip( self, yMin, yMax )
        return( self.returnAsClass( self, pointwiseXY.clip( self, yMin, yMax ) ) )

    def convertAxisToUnit( self, indexOrName, newUnit ) :

        index = self.getAxisIndexByIndexOrName( indexOrName )
        axis = self.axes[index]
        factor = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( '1 ' + axis.getUnit( ) ).getValueAs( newUnit )
        data = []
        for xy in self :
            xy[index] *= factor
            data.append( xy )
        n = return_pointwiseXY_AsXYs( self, data, units = { index : newUnit }, template = self )
        return( self.returnAsClass( self, n, axes_ = n.axes ) )

    def getAxisIndexByIndexOrName( self, indexOrName ) :

        if( type( indexOrName ) == type( 1 ) ) :
            return( indexOrName )
        elif( type( indexOrName ) == type( '' ) ) :
            for index, axis in enumerate( self.axes ) :
                if( axis.getLabel( ) == indexOrName ) : return( index )
        raise TypeError( 'argument must be an integer of a string, not type %s' % brb.getType( indexOrName ) )

    def copy( self, parent = None, index = None, value = None, moniker = xData, axes_ = None, isPrimaryXData = None ) :

        c = pointwiseXY.copy( self )
        if( axes_ is None ) : axes_ = self.axes
        if( parent is not None ) :
            import regions
            if( isinstance( parent, regions.regions ) ) :
                n = len( parent.axes.getRootParent( ).axes )
                axes_ = axes.interpolationAxes( n - 2, axes.interpolationXY( self.axes[0].interpolation.independent, self.axes[0].interpolation.dependent ), self )
        return( self.returnAsClass( self, c, parent = parent, index = index, value = value, axes_ = axes_, moniker = moniker, isPrimaryXData = isPrimaryXData ) )

    def copyDataToXYs( self, xUnit = None, yUnit = None ) :

        xScale, yScale = 1.0, 1.0
        if( xUnit is not None ) : xScale = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( '1 ' + self.axes[0].getUnit( ) ).getValueAs( xUnit )
        if( yUnit is not None ) : yScale = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( '1 ' + self.axes[1].getUnit( ) ).getValueAs( yUnit )
        return( pointwiseXY.copyDataToXYs( self, xScale = xScale, yScale = yScale ) )

    def dullEdges( self, lowerEps = 0., upperEps = 0., positiveXOnly = 0 ) :

        d = pointwiseXY.dullEdges( self, lowerEps = lowerEps, upperEps = upperEps, positiveXOnly = positiveXOnly );
        return( self.returnAsClass( self, d ) )

    def getValue_units( self, x ) :

        x = evaluateValueAsUnit( self.axes[0].getUnit( ), x )
        y = pointwiseXY.getValue( self, x )
        if( self.axes[1].getUnit( ) == '' ) : return( y )
        return( physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( y, self.axes[1].getUnit( ) ) )

    def setValue_units( self, x, y ) :

        x = evaluateValueAsUnit( self.axes[0].getUnit( ), x )
        y = evaluateValueAsUnit( self.axes[1].getUnit( ), y )
        pointwiseXY.setValue( self, x, y )

    def mutualify( self, lowerEps1, upperEps1, positiveXOnly1, other, lowerEps2, upperEps2, positiveXOnly2 ) :
        '''
        .. note:: Need to check that x units are the same.
        '''
        m1, m2 = pointwiseXY.mutualify( self, lowerEps1, upperEps1, positiveXOnly1, other, lowerEps2, upperEps2, positiveXOnly2 )
        return( self.returnAsClass( self, m1 ), other.returnAsClass( other, m2 ) )

    def normalize( self, insitu = False ) :

        n = pointwiseXY.normalize( self )
        if( insitu ) :
            pointwiseXY.setData( self, n )
            return( self )
        return( self.returnAsClass( self, n ) )

    def thicken( self, sectionSubdivideMax = 1, dxMax = 0., fxMax = 1. ) :
        '''
        .. note:: Need unit for dxMax.
        '''
        t = pointwiseXY.thicken( self, sectionSubdivideMax = sectionSubdivideMax, dxMax = dxMax, fxMax = fxMax )
        return( self.returnAsClass( self, t ) )

    def thin( self, accuracy ) :

        return( self.returnAsClass( self, pointwiseXY.thin( self, accuracy ) ) )

    def trim( self ) :

        return( self.returnAsClass( self, pointwiseXY.trim( self ) ) )

    def union( self, other, fillWithSelf = 1, trim = 0 ) :

        raiseNotSameUnits( self, other, checkXOnly = True )
        t = pointwiseXY.union( self, other, fillWithSelf = fillWithSelf, trim = trim  )
        return( self.returnAsClass( self, t ) )

    def xMin( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ) )

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( '1 ' + self.getDomainUnit( ) ).getValueAs( unitTo ) )

    def domainMin( self, unitTo = None, asPQU = False ) :

        return( physicalQuantityWithUncertainty.valueOrPQ( pointwiseXY.xMin( self ), unitFrom = self.axes[0].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def xMax( self, unitTo = None, asPQU = False ) :

        return( self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def domainMax( self, unitTo = None, asPQU = False ) :

        return( physicalQuantityWithUncertainty.valueOrPQ( pointwiseXY.xMax( self ), unitFrom = self.axes[0].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def getDomain( self, unitTo = None, asPQU = False ) :

        return( self.domainMin( unitTo = unitTo, asPQU = asPQU ), self.domainMax( unitTo = unitTo, asPQU = asPQU ) )

    def getDomainGrid( self, unitTo = None ) :

        scale = self.domainUnitConversionFactor( unitTo )
        return( pointwiseXY.getDomainGrid( self, scale ) )

    def getDomainUnit( self ) :

        return( self.axes[0].getUnit( ) )

    def yMin( self, unitTo = None, asPQU = False ) :

        return( physicalQuantityWithUncertainty.valueOrPQ( pointwiseXY.yMin( self ), unitFrom = self.axes[1].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def yMax( self, unitTo = None, asPQU = False ) :

        return( physicalQuantityWithUncertainty.valueOrPQ( pointwiseXY.yMax( self ), unitFrom = self.axes[1].getUnit( ), unitTo = unitTo, asPQU = asPQU ) )

    def xSlice( self, xMin = None, xMax = None, fill = 1, dullEps = 0. ) :
        '''
        Returns a new instance with self sliced between ``xMin`` and ``xMax``.

        :param xMin:    [optional] the lower x-value of the slice, default is xMin of self,
        :param xMax:    [optional] the upper x-value of the slice, default is xMax of self,
        :param fill:    [optional] if True, points are added at xMin and xMax if they are not in self, else only existing points in the range [xMin, xMax] are included.
        :param dullEps: [optional] (Currently not implemented) the lower and upper points are dulled, default is 0.
        '''

        if( xMin is None ) : xMin = self.xMin( )
        if( xMax is None ) : xMax = self.xMax( )
        s = pointwiseXY.xSlice( self, xMin = xMin, xMax = xMax, fill = fill, dullEps = dullEps )
        return( self.returnAsClass( self, s ) )

    def xSlice_units( self, xMin = None, xMax = None, fill = 1, dullEps = 0. ) :
        '''
        Same as xSlice, only xMin and xMax are now PhysicalQuantitiesWithUncertainty and have units.
        '''

        if( xMin is None ) : xMin = self.xMin( asPQU=True )
        if( xMax is None ) : xMax = self.xMax( asPQU=True )
        unit = self.axes[0].getUnit( )
        xMin = evaluateValueAsUnit( unit, xMin )
        xMax = evaluateValueAsUnit( unit, xMax )
        s = pointwiseXY.xSlice( self, xMin = xMin, xMax = xMax, fill = fill, dullEps = dullEps )
        return( self.returnAsClass( self, s ) )

    def __mod__( self, other ) : raise NotImplementedError( 'Currently, mod is not implemented' )

    def __pow__( self, other ) : raise NotImplementedError( 'Currently, pow is not implemented' )

    def __exp__( self, other ) : raise NotImplementedError( 'Currently, __exp__ is not implemented' )

    def convolute( self, func, fpars ) : 
        raise NotImplementedError( 'Currently, convolute is not implemented' )

    def _group( self, groupBoundaries, f2 = None, f3 = None, norm = None, asXYs = False ) :
        """
        This function will be renamed to group when it is completed. Called ``_group`` for now so that it does not conflict with
        Dave Brown's version when using svn. That is, this version will eventually replace Dave's.

        This function groups self onto the boudaries given in ``groupBoundaries``. If ``asXYs`` is ``False``, the returned 
        list, here in called ``groups``, has ``length = len( groupBoundaries ) - 1`` with the :math:`{i-1}^{th}` item begin the 
        integral along the :math:`x` axis, :math:`dx`, from ``x[i] = groupBoundaries[i]`` to ``x[i+1] = groupBoundaries[i+1]`` defined as
        ::

            groups[i] = self(x) / n[i]                       for f2 = f3 = None
            groups[i] = self(x) * f2(x) / n[i]               for f2 = None
            groups[i] = self(x) * f3(x) / n[i]               for f3 = None
            groups[i] = self(x) * f2(x) * f(x) / n[i]        otherwise

        Here norm can be None, ``dx`` or a list of length ``len( groupBoundaries ) - 1`` of numbers.
        
            +---------------+----------------+
            | norm          |  n[i]          |
            +---------------+----------------+
            | None          |  1             |
            +---------------+----------------+
            | 'dx'          |  x[i+1] - x[i] |
            +---------------+----------------+
            | a python list |  norm[i]       |
            +---------------+----------------+

        If ``asXYs`` is ``True``, an instance of class ``XYs`` is returned with the x-values from ``groupBoundaries`` and with 
        interpolation 'linear,flat'.
        """

        accuracy, yUnit = self.getAccuracy( ), PhysicalQuantityWithUncertainty( 1, self.axes[1].getUnit( ) )
        if( f2 is None ) :
            if( f3 is None ) : 
                grouped = self.groupOneFunctions( groupBoundaries, norm = norm )
            else :
                grouped = self.groupTwoFunctions( groupBoundaries, f3, norm = norm )
                yUnit = yUnit * PhysicalQuantityWithUncertainty( 1,  f3.axes[1].getUnit( ) )
                accuracy = max( accuracy, f3.getAccuracy( ) )
        else :
            yUnit = yUnit * PhysicalQuantityWithUncertainty( 1, f2.axes[1].getUnit( ) )
            accuracy = max( accuracy, f2.getAccuracy( ) )
            if( f3 is None ) :
                grouped = self.groupTwoFunctions( groupBoundaries, f2, norm = norm )
            else :
                grouped = self.groupThreeFunctions( groupBoundaries, f2, f3, norm = norm )
                yUnit = yUnit * PhysicalQuantityWithUncertainty( 1, f3.axes[1].getUnit( ) )
                accuracy = max( accuracy, f3.getAccuracy( ) )
        if( norm is None ) :
            yUnit = PhysicalQuantityWithUncertainty( 1, self.axes[0].getUnit( ) ) * yUnit
        elif( norm != 'dx' ) :
            pass                    # Need to add units to norm. That is, norm, grouped and groupBoundaries should be an instance of Ys.
        if( asXYs ) :
            grouped.append( grouped[-1] )
            axes_ = axes.defaultAxes( dependentInterpolation = axes.flatToken, 
                labelsUnits = { 0 : [ self.axes[0].getLabel( ), self.axes[0].getUnit( ) ], 1 : [ "", yUnit.getUnitName( ) ] } )
            grouped = XYs( axes_, [ groupBoundaries, grouped ], accuracy, dataForm = 'xsandys' )
        return( grouped )
        
    def group( self, groupBoundaries = [], groupUnit = '' ):
        '''
        A convience fuction that groups an XYs instance.  
        
        :param self: the ``XYs`` instance to group
        :param groupBoundaries:  a list of group boundaries in ascending order
        :param groupUnits: the units the group boundaries are given in
        
        Returns a copy of self, but in the new group structure, 
        in the group structure's units (at least for the X axis).
            
        Note: We do not do flux weighting!!!

        This version uses Bret's norm='dx' keyword
        '''
        scaled = self.convertAxisToUnit( 0, groupUnit )
         
        # New axes are same as old ones, just with flat interpolation
        newAxes = scaled.axes.copy( )
        newAxes.axes[0].interpolation.dependent = axes.flatToken
    
        # Make the new instance with the grouped data
        grouped = XYs( newAxes, [ ( groupBoundaries[i], y ) for i,y in enumerate( scaled.groupOneFunction( groupBoundaries, norm = 'dx' ) ) ], 0.0001 )
        
        # Last point gets lost, so put it back
        grouped.setValue( scaled[-1][0], grouped[-1][1] )
        
        return grouped
        
    def groupOneFunction( self, groupBoundaries, norm = None ) :
        '''.. note:: Need unit of groupBoundaries.'''
        return( pointwiseXY.groupOneFunction( self, groupBoundaries, norm = norm ) )

    def groupTwoFunctions( self, groupBoundaries, f2, norm = None ) :
        '''.. note:: Need unit of groupBoundaries.'''
        return( pointwiseXY.groupTwoFunctions( self, f2, groupBoundaries, norm = norm ) )

    def groupThreeFunctions( self, groupBoundaries, f2, f3, norm = None ) :
        '''.. note:: Need unit of groupBoundaries.'''
        return( pointwiseXY.groupThreeFunctions( self, f2, f3, groupBoundaries, norm = norm ) )

    def integrate( self, xMin = None, xMax = None ) :
        '''
        Definite integral of current ``XYs`` instance from ``xMin`` to ``xMax``:
        
        .. math::
            \int_{xMin}^{xMax} dx \; XYs(x)

        If ``xMin`` or ``xMax`` are unspecified, they are taken from domain of the ``XYs`` instance.
        '''

        if( xMin is None ) : xMin = self.xMin( )
        if( xMax is None ) : xMax = self.xMax( )
        xMin = max( xMin, self.xMin( ) )
        xMax = min( xMax, self.xMax( ) )
        return( pointwiseXY.integrate( self, xMin = xMin, xMax = xMax ) )

    def indefiniteIntegral( self, xMin = None, xMax = None ) :
        '''
        Indefinite integral of current XYs instance:
        
        .. math::
            \int_0^x dx \; XYs(x)
            
        The new ``XYs`` instance is defined on the range of the old one and the units are wrong.
        '''
        myAxes = self.axes
        myData = [ [ self.xMin( ), 0.0 ] ]
        for i in range( len( self ) - 1 ):
            xMin = self[i][0]
            xMax = self[i+1][0]
            myData.append( [ xMax, myData[-1][1]+pointwiseXY.integrate( self, xMin = xMin, xMax = xMax ) ] )
        return XYs( myAxes, myData, 1e-6 )

    def integrateTwoFunctions( self, f2, xMin = None, xMax = None ) :
        '''.. note:: Need unit of f2.'''

        if( xMin is None ) : xMin = self.xMin( )
        if( xMax is None ) : xMax = self.xMax( )
        xMin = max( xMin, self.xMin( ), f2.xMin( ) )
        xMax = min( xMax, self.xMax( ), f2.xMax( ) )
        return( pointwiseXY.groupTwoFunctions( self, f2, [ xMin, xMax ] )[0] )

    def integrateThreeFunctions( self, f2, f3, xMin = None, xMax = None ) :
        '''.. note:: Need units for f2 and f3.'''
        if( xMin is None ) : xMin = self.xMin( )
        if( xMax is None ) : xMax = self.xMax( )
        xMin = max( xMin, self.xMin( ), f2.xMin( ), f3.xMin( ) )
        xMax = min( xMax, self.xMax( ), f2.xMax( ), f3.xMax( ) )
        return( pointwiseXY.groupThreeFunctions( self, f2, f3, [ xMin, xMax ] )[0] )

    def toPointwiseLinear( self, accuracy = None, lowerEps = 0, upperEps = 0, cls = None ) :

        return( self.changeInterpolation( axes.linearToken, axes.linearToken, accuracy = self.getAccuracy( ), lowerEps = lowerEps, upperEps = upperEps,
            cls = cls ) )

    def toXML( self, tag = 'xData', indent = '', incrementalIndent = '  ', pairsPerLine = 10, xyFormatter = None, xySeparater = ' ' ) :

        return( '\n'.join( self.toXMLList( indent = indent, incrementalIndent = incrementalIndent, pairsPerLine = pairsPerLine, xyFormatter = xyFormatter, 
            xySeparater = xySeparater ) ) )

    def toXMLList( self, tag = None, indent = '', incrementalIndent = '  ', pairsPerLine = 10, xyFormatter = None, xySeparater = ' ', oneLine = False ) :

        from pqu.physicalQuantityWithUncertainty import toShortestString
        if( hasattr( self, 'tag' ) ) :
            tag = self.tag
        elif( tag is None ) :
            tag = self.moniker
        indent2 = indent + incrementalIndent
        indexValueStr = ""
        if( self.value is not None ) : indexValueStr = ' value="%s"' % self.value
        if( self.index is not None ) : indexValueStr += ' index="%s"' % self.index
        xDataString = ''
        if( hasattr( self, 'isPrimaryXData' ) ) :
            if( self.isPrimaryXData ) : xDataString = ' xData="%s"' % self.xData
        s = [ '%s<%s%s%s length="%s" accuracy="%s">' % ( indent, tag, indexValueStr, xDataString, len( self ), self.getAccuracy( ) ) ] 
        s += self.axes.toXMLList( indent = indent2 )
        xyStrs = []
        for x, y in self :
            if( xyFormatter is None ) :
                xyStrs.append( '%s %s' % ( toShortestString(x), toShortestString(y) ) )
            else :
                xyStrs.append( xyFormatter( x, y ) )

        if( oneLine ) :
            s[-1] += ' ' + xySeparater.join( xyStrs )
        else :
            n, xySubStrs, indent2Tag = 0, [], indent2 + '<data> '
            for xyStr in xyStrs :
                xySubStrs.append( xyStr )
                if( len( xySubStrs ) >= pairsPerLine ) :
                    s.append( indent2Tag + xySeparater.join( xySubStrs ) )
                    xySubStrs, indent2Tag = [], indent2 + '       '
            if( len( xySubStrs ) > 0 ) : s.append( indent2Tag + xySeparater.join( xySubStrs ) )
            s[-1] += '</data>'
        s[-1] += '</%s>' % tag
        return( s )

    def plot( self, xylog = 0, xMin = None, xMax = None, yMin = None , yMax = None, title = '' ) :

        from fudge.core import fudgemisc
        from fudge.core.utilities import fudgeFileMisc
        from fudge.vis.gnuplot import plotbase
        import os

        def getUnitlessNumber( value, unit, default ) :

            if( value is None ) : value = default
            if( fudgemath.isNumber( value ) ) : return( float( value ) )
            if( type( '' ) == type( value ) ) :
                try :
                    return( float( value ) )
                except :
                    print 'value = "%s"' % value
                    value = physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty( value )
            if( isinstance( value, physicalQuantityWithUncertainty.PhysicalQuantityWithUncertainty ) ) : return( value.getValueAs( unit ) )
            raise Exception( 'Cannot convert %s to a unitless number' % str( value ) )

        xLabel = self.axes[0].plotLabel( )
        yLabel = self.axes[1].plotLabel( )

        xMin = getUnitlessNumber( xMin, self.axes[0].getUnit( ), self.xMin( ) )
        xMax = getUnitlessNumber( xMax, self.axes[0].getUnit( ), self.xMax( ) )
        yMin = getUnitlessNumber( yMin, self.axes[1].getUnit( ), self.yMin( ) )
        yMax = getUnitlessNumber( yMax, self.axes[1].getUnit( ), self.yMax( ) )

        dt = plotbase.parsePlotOptions( xMin, xMax, yMin, yMax, xLabel, yLabel, title )
        f = fudgeFileMisc.fudgeTempFile( )
        f.write( self.toString( ) )
        f.close( )
        p = os.path.join( os.path.realpath( __file__ ).split( '/fudge/core/' )[0], "fudge", "vis", "gnuplot", "endl2dplot.py" )
        s = [ "python", p, 'xylog', str( xylog ) ] + dt + [ f.getName( ) ]
        subprocessing.spawn( s )

    @classmethod
    def returnAsClass( cls, self, other, parent = None, index = None, value = None, axes_ = None, moniker = xData, isPrimaryXData = None ) :
        """
        Returns other as a class of cls. Other must be an instance derived from the class pointwiseXY.
        cls must be a sub-class of XYs. If index and axes are not specified, they are taken from self. The method 
        changeInterpolationSubFunction is taken from self. The main use of this classmethod is for methods like __add__ where
        the addends may be a class derived from XYs. For example, the crossSection.pointwise class is derived from The XYs class.
        If two instances xSec1 and xSec2 of the crossSection.pointwise class are added together (i.e., xSec1 + xSec2) then, since
        the __add__ method used this classmethod, the returned instance will also be an instance of crossSection.pointwise class.
        """

        if( index is None ) : index = self.index
        if( axes_ is None ) : axes_ = self.axes
        return( cls( axes_, other, other.getAccuracy( ), overflowSize = 10, biSectionMax = other.getBiSectionMax( ),
            infill = other.getInfill( ), safeDivide = other.getSafeDivide( ), index = index, value = value, parent = parent, 
            moniker = moniker, isPrimaryXData = isPrimaryXData, changeInterpolationSubFunction = self.changeInterpolationSubFunction ) )

    @classmethod
    def parseXMLNode( cls, xdataElement, linkData={} ):
        """ translate anything directly inheriting XYs (pointwise multiplicity or xsc for example) back from xml """
        attrs = dict( xdataElement.items() )
        if "xData" in attrs:
            assert attrs.pop("xData") == monikerXYs
            attrs['isPrimaryXData'] = True
        else: attrs['isPrimaryXData'] = False
        if 'length' in attrs: del attrs['length']
        if 'index' in attrs: attrs['index'] = int(attrs['index'])
        attrs['moniker'] = xdataElement.tag
        axes_ = axes.parseXMLNode( xdataElement[0] )
        if axes_[0].interpolation.dependent==axes.chargedParticleToken:
            from fudge.gnd.reactionData import crossSection
            attrs['changeInterpolationSubFunction'] = crossSection.chargeParticle_changeInterpolationSubFunction
        data = map(float, xdataElement[1].text.split())
        data = zip( data[::2], data[1::2] )
        return cls( axes_, data, float(attrs.pop('accuracy')), **attrs )

    @staticmethod
    def createFromFunction( axes_, Xs, func, parameters, accuracy, biSectionMax, checkForRoots = False, infill = 1, safeDivide = 1 ) :
        """Given an ascending list of x-values (Xs) and a function (func), create a linear pointwise representation 
        of the function over the domain of the x-values. There must be at least 2 x-values in Xs. The axes interpolation
        must be linear,linear, else a raise is executed. See pointwiseXY_C.createFromFunction for all other arguments."""

        import math
        if( not( axes_.isLinear( ) ) ) : raise Exception( "Interpolation can only be '%s,%s', not %s" % ( axes.linearToken, axes.linearToken, 
            `axes_.interpolation` ) )
        xys = pointwiseXY_C.createFromFunction( Xs, func, parameters, accuracy, biSectionMax, checkForRoots = checkForRoots, infill = infill, safeDivide = safeDivide )
        biSectionMax = max( 0, biSectionMax - math.log( len( xys ) / len( Xs ) ) / math.log( 2 ) )
        return( XYs( axes_, xys, accuracy, infill = infill, safeDivide = safeDivide ) )

if( __name__ == '__main__' ) :

    vl1 = axes.axes( )
    vl1[0] = axes.axis( 'energy_in', 0, 'eV', frame = axes.labToken, interpolation = axes.interpolationXY( axes.linearToken, axes.linearToken ) )
    vl1[1] = axes.axis( 'crossSection', 1, 'b', frame = axes.labToken )
    pXY1 = XYs( vl1, [ [ 1, 0 ], [ 3, 2 ], [ 4, 1 ] ], 1e-3, safeDivide = True, biSectionMax = 7 )
    pXY2 = XYs( vl1, [ [ 1, 0 ], [ 4, 2 ] ], 1e-3, safeDivide = True, biSectionMax = 7 )


    print pXY1.axes
    print pXY1.toString( )
    print (pXY1 + '320 b').toString( )
    print
    print pXY2.toString( )
    print
    print
    print (pXY2 / pXY1).toString( )

    pXY1 = XYs( vl1, [ [ 1, 1 ], [ 2, 2 ], [ 4, 3 ], [ 7, 1 ] ], 1e-3, safeDivide = True, biSectionMax = 7 )
    print pXY1.toString( )

    print
    print (5. / pXY1).toString( )

    pXY2 = XYs( vl1, [ [ 1, 0 ], [ 2.5, 2 ], [ 3, 3 ], [ 7, 0 ] ], 1e-3, safeDivide = True, biSectionMax = 7 )
    print
    print pXY2.toString( )

    print
    print (3. / pXY2).toString( )

    print
    print (pXY1 * pXY2).toString( )

    mul = pXY1 * pXY2
    mul = XYs( vl1, mul, 1e-3, safeDivide = True, biSectionMax = 0 )
    print
    print mul.toXML( indent = '    ', pairsPerLine = 5 )

    tXY1 = pXY1.thicken( 10, 0, 1 )
    tXY2 = pXY2.thicken( 10, 0, 1 )
    print
    print (tXY1 * tXY2).toString( )

    print
    print (pXY1 * pXY2 - tXY1 * tXY2).toString( )

    print
    r = pXY1 * pXY2 / ( tXY1 * tXY2 )
    print r.toXML( indent = '    ', pairsPerLine = 5 )
    rm1 = pXY1 * pXY2 / ( tXY1 * tXY2 ) - 1
    print rm1.toString( )
