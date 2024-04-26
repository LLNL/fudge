# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import abc
import math

from pqu import PQU as PQUModule

from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from . import enums as enumsModule
from . import XYs1d as XYs1dModule
from . import series1d as series1dModule
from . import base as baseModule

domainEpsilon = 1e-15

class Regions( baseModule.XDataFunctional ) :
    """Abstract base class for Regions."""

    ancestryMembers = baseModule.XDataFunctional.ancestryMembers                # See note in multiD_XYs.py about setting ancestryMembers.

    def __init__( self, axes = None, index = None, valueType = enumsModule.ValueType.float64, outerDomainValue = None, label = None ) :

        baseModule.XDataFunctional.__init__(self, axes, index=index, valueType=valueType, outerDomainValue=outerDomainValue, label=label)
        self.__regions = []

    def __len__( self ) :
        """Returns the number of regions in self."""

        return( len( self.__regions ) )

    @staticmethod
    @abc.abstractmethod
    def allowedSubElements( ):
        pass

    def __getitem__( self, index ) :
        """Returns the (i-1)^th region of self."""

        return( self.__regions[index] )

    def __setitem__( self, index, region ) :
        """
        Set region (the right-hand-side) to the :math:`{i-1}^{th}` region of self. Region must be an instance of 
        base.XDataFunctional with the same dimension as self.
        If :math:`i > 0`, the following must also be met:
        
            - all the prior regions must already exists,
            - region's minimum domain value must be equal to the prior region's maximum domain value.
        """

# BRB need to check axes.
        if( not( isinstance( region, self.allowedSubElements( ) ) ) ) : raise TypeError( 'Invalid class for insertion: %s' % region.__class__ )

        n1 = len( self )
        if( index < 0 ) : index += n1
        if( not( 0 <= index <= n1 ) ) : raise IndexError( 'Index = %s not in range 0 <= index <= %d' % ( index, n1 ) )
        if( index > 0 ) :
            if( not( math.isclose( self.__regions[index-1].domainMax, region.domainMin ) ) ) :
                raise ValueError( "Prior region's domainMax %s != new region's domainMin = %s" % ( self.__regions[index-1].domainMax, region.domainMin ) )
        if( ( n1 > 0 ) and ( index < ( n1 - 1 ) ) ) :
            if( not math.isclose( self.__regions[index+1].domainMin, region.domainMax ) ) :
                raise ValueError(  "Next region's domainMin %s != new region's domainMax = %s" % ( self.__regions[index-1].domainMin, region.domainMax ) )

        if( index == n1 ) :
            self.__regions.append( region )             # Append to the end.
        else :
            self.__regions[index] = region              # Replaces the current contents of index with region.

        region.setAncestor( self )
        region.index = index

    def __add__(self, other):

        self2, other2 = self.copyToCommonRegions(other)
        for i1, region in enumerate(self2):
            self2[i1] = region + other2[i1]

        return self2

    __radd__ = __add__

    def __sub__( self, other ) :

        self2, other2 = self.copyToCommonRegions( other )
        for i1, region in enumerate( self2 ) : self2[i1] = region - other2[i1]
        return( self2 )

    @property
    def domainMin( self ) :

        return( self.__regions[0].domainMin )

    @property
    def domainMax( self ) :

        return( self.__regions[-1].domainMax )

    @property
    def domainGrid( self ) :

        grid = set()
        for region in self: grid.update( region.domainGrid )
        return sorted( grid )

    @property
    def domainUnit( self ) :

        return( self.getAxisUnitSafely( self.dimension ) )

    @property
    def rangeMin( self ) :

        return( min( [ region.rangeMin for region in self ] ) )

    @property
    def rangeMax( self ) :

        return( max( [ region.rangeMax for region in self ] ) )

    @property
    def rangeUnit( self ) :

        return( self.getAxisUnitSafely( 0 ) )

    @property
    def regions( self ) :
        """Returns self's __region."""

        return( self.__regions )

    @property
    def functionNdsName( self ) :
        """Returns the node name for the child "function#ds"."""

        return( "function%dds" % self.dimension )

    def append(self, curve):
        '''
        Appends *curve* to the end of *self*. If *curve* is a **Regions** instance, then **append** is call for 
        each region in curve.
        '''

        if isinstance(curve, Regions):
            for region in curve:
                self.append(region.copy())
        else:
            self[len(self)] = curve
            curve.keyName = 'index'

    def prepend(self, curve):
        '''
        Prepends *curve* to the beginning of *self*. If *curve* is a **Regions** instance, then **append** is call 
        for each region in curve.
        '''

        if isinstance(curve, Regions):
            for index in range(len(curve), 0, -1):
                self.prepend(curve[index].copy())
        else:
            if not isinstance(curve, self.allowedSubElements()):
                raise TypeError('Invalid class for insertion: %s' % region.__class__)

            if len(self) > 0:
                if not math.isclose(curve.domainMax, self[0].domainMin):
                    raise ValueError('''Prepending region's domainMax %s != first region's domainMin = %s''' % (region.domainMax, self[0].domainMin))

            self.__regions.insert(0, curve)
            curve.setAncestor(self)
            curve.index = 0

    def convertUnits( self, unitMap ) :

        factors = self.axes.convertUnits( unitMap )
        for region in self : region.convertUnits( unitMap )
        self.fixValuePerUnitChange( factors )

    def copy( self ) :
        # FIXME some of this should probably move to 'returnAsClass' method

        axes = self.axes
        if( axes is not None ) : axes = axes.copy( )
        newRegions = self.__class__( axes = axes, index = self.index, valueType = self.valueType, outerDomainValue = self.outerDomainValue, label = self.label )
        for child in self : newRegions.append( child.copy( ) )
        return( newRegions )

    def hasData(self):
        '''
        Returns True if one of self's regions hasData returns True and False otherwise.
        '''

        for region in self:
            if region.hasData():
                return True

        return False

    def splitInTwo( self, domainValue, epsilon = domainEpsilon ) :
        """
        Splits the region containing domainValue into two regions.
        """

        for i1, region in enumerate( self ) :
            domainMin, domainMax = region.domainMin, region.domainMax
            if( domainMin < domainValue < domainMax ) :
                r1, r2 = region.splitInTwo( domainValue, epsilon = domainEpsilon )
                self.__regions[i1] = r2
                r1.setAncestor( self )
                self.__regions.insert( i1, r1 )
                r2.setAncestor( self )
                return

    def domainUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.domainUnit ).getValueAs( unitTo ) )

    def domainSlice(self, domainMin=None, domainMax=None, fill=1, dullEps=0.0, **kwargs):
        """
        Returns a new instance with self sliced between ``domainMin`` and ``domainMax``.
        Result may be a regions container, or it may be XYs, multiD_XYs, etc.

        :param domainMin:   [optional] the lower x-value of the slice, default is domain minimum of self,
        :param domainMax:   [optional] the upper x-value of the slice, default is domain maximum of self,
        :param fill:        [optional] if True, points are added at domainMin and domainMax if they are not in self,
                                       else only existing points in the range [domainMin, domainMax] are included.
        :param dullEps:     [optional] (Currently not implemented) the lower and upper points are dulled, default is 0.
        """

        if( domainMin is None ) : domainMin = self.domainMin
        domainMin = max( domainMin, self.domainMin )
        if( domainMax is None ) : domainMax = self.domainMax
        domainMax = min( domainMax, self.domainMax )

        for ridx1, region in enumerate( self ) :
            if( region.domainMax  > domainMin ) : break
        for ridx2, region in enumerate( self ) :
            if( region.domainMax >= domainMax ) : break

        if ridx1 == ridx2:  # only one region left after slicing, return as XYs1d or multiD_XYs
            return self[ridx1].domainSlice( domainMin=domainMin, domainMax=domainMax, fill=fill, dullEps=dullEps )
        else:
            newRegions = self.__class__( axes = self.axes, index = self.index, valueType = self.valueType, outerDomainValue = self.outerDomainValue, label = self.label )
            newRegions.append( self[ridx1].domainSlice( domainMin = domainMin, domainMax = self[ridx1].domainMax, fill = fill, dullEps = dullEps ) )
            for idx in range(ridx1+1,ridx2):
                newRegions.append(self[idx])
            newRegions.append( self[ridx2].domainSlice(
                domainMin=self[ridx2].domainMin, domainMax=domainMax, fill=fill, dullEps=dullEps ) )

            return newRegions

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """
        Sets domain minimum and maximum per the arguments.
        """

        OldDomainMin = self.domainMin
        OldDomainMax = self.domainMax

        if fixToDomain == enumsModule.FixDomain.lower:
            regions1 = self.domainSlice(domainMin = domainMin, fill = True)
        elif fixToDomain == enumsModule.FixDomain.upper:
            regions1 = self.domainSlice(domainMax = domainMax, fill = True)
        else:
            regions1 = self.domainSlice(domainMin = domainMin, domainMax = domainMax, fill = True)

        self.__regions = []
        if isinstance(regions1, Regions):
           for index, region in enumerate(regions1): self[index] = region
        else:
            if len(regions1) == 0: return 1             # Caller beware.
            self[0] = regions1

        if OldDomainMin == self.domainMin and OldDomainMax == self.domainMax: return 0
        return 1

    def rangeUnitConversionFactor( self, unitTo ) :

        if( unitTo is None ) : return( 1. )
        return( PQUModule.PQU( '1 ' + self.rangeUnit ).getValueAs( unitTo ) )

    def evaluate( self, domainValue ) :
        """
        Evaluate function at requested domain point. If at discontinuity, return upper region's value.

        :param domainValue:
        :return interpolated point at domainValue:
        """

        for region in self.__regions :
            if( domainValue < region.domainMax ) : return( region.evaluate( domainValue ) )

        return( self.__regions[-1].evaluate( domainValue ) )  # Domain value is above the last region. Let the last region determine what to do.

    def findInstancesOfClassInChildren( self, cls, level = 9999 ) :
        """
        Finds all instances of class *cls* in self's children, grand-children, etc.
        """

        foundInstances = []
        level -= 1
        if( level < 0 ) : return( foundInstances )
        for region in self :
            if( isinstance( region, cls ) ) : foundInstances.append( region )
            foundInstances += region.findInstancesOfClassInChildren( cls, level = level )

        return( foundInstances )

    def integrate( self, **limits ):
        """
        Integrate a piecewise function. Supports limits for each axis.
        Example::

            regions.integrate( energy_in = ('1e-5 eV', '10 eV'), energy_out = ('1 keV', '10 keV') )

        :param limits: dictionary containing limits for each independent axis (keyed by axis label or index).
        If an independent axis is missing from the dictionary, integrate over the entire domain of that axis.

        :return: float
        """

        from . import multiD_XYs as multiD_XYsModule

        integral = 0
        for region in self:
            if( isinstance( region, ( XYs1dModule.XYs1d, series1dModule.Series1d ) ) ) :
                Min, Max = None, None
                if self.axes[-1].label in limits:
                    Min, Max = limits.pop( self.axes[-1].label )
                elif self.axes[-1].index in limits:
                    Min, Max = limits.pop( self.axes[-1].index )

                integral += region.integrate( domainMin = Min, domainMax = Max )
            elif isinstance( region, multiD_XYs1dModule.XYsnd ):
                integral += region.integrate( **limits )
            else:
                raise TypeError( "Unsupported class for integration: %s" % type( region ) )
        return( integral )

    def toString( self ) :

        s1 = ''
        for i1, region in enumerate( self ) : s1 += '# region %s\n' % i1 + region.toString( )
        return( s1 )

    def toXML_strList(self, indent = '', **kwargs):

        formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        indent3 = indent2 + kwargs.get('incrementalIndent', '  ')
        if formatVersion == GNDS_formatVersionModule.version_1_10: indent3 = indent2

        attributeStr = baseModule.XDataFunctional.attributesToXMLAttributeStr(self)
        XML_strList = [ '%s<%s%s>' % ( indent, self.moniker, attributeStr ) ]
        if self.isPrimaryXData():
            if self.axes is not None: XML_strList += self.axes.toXML_strList(indent2, **kwargs)

        if formatVersion != GNDS_formatVersionModule.version_1_10: XML_strList.append('%s<%s>' % ( indent2, self.functionNdsName ))
        for region in self.__regions: XML_strList += region.toXML_strList(indent3, **kwargs)
        if formatVersion != GNDS_formatVersionModule.version_1_10: XML_strList[-1] += "</%s>" % self.functionNdsName

        if self.uncertainty: XML_strList += self.uncertainty.toXML_strList(indent2, **kwargs)
        XML_strList[-1] += '</%s>' % self.moniker

        return XML_strList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        attributes, extraAttributes = baseModule.XDataFunctional.parseBareNodeCommonAttributes(node, xPath)     # parseBareNodeCommonAttributes adds to xPath.
        if len(extraAttributes) > 0: raise Exception('Invalid attributes: %s.' % ( ', '.join(list(extraAttributes.keys())) ))

        regions1 = cls(**attributes)

        extraNodes = baseModule.XDataFunctional.parseNodeStandardChildren(regions1, node, xPath, linkData, **kwargs)

        if len(extraNodes) > 0:                                                     # The next few line support GNDS 1.10 and 2.0.
            for index, extraNode in enumerate(extraNodes):
                if extraNode.tag == regions1.functionNdsName: break
            if extraNode.tag == regions1.functionNdsName: extraNodes = extraNodes.pop(index)

        if 'axes' not in kwargs: kwargs['axes'] = regions1.axes
        allowedSubElements = cls.allowedSubElements()
        for child in extraNodes:
            childClass = None
            for allowedChildClass in allowedSubElements:
                if allowedChildClass.moniker == child.tag:
                    childClass = allowedChildClass
                    break
            if childClass is None: raise TypeError('Invalid child "%s" for node "%s"' % (child.tag, cls.moniker))
            xdata = childClass.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            regions1.append(xdata)

        xPath.pop()

        return regions1

class Regions1d( Regions ) :

    moniker = 'regions1d'
    dimension = 1

    def __neg__(self):
        """
        This method returns self with all regions negated.
        """

        regions1d = self.__class__(axes=self.axes)
        for region in self:
            regions1d.append(-region)

        return regions1d

    def __mul__( self, other ) :

        _self, _other = self.copyToCommonRegions( other )
        _regions1d = self.__class__( )
        for index, region1 in enumerate( _self ) :
            region2 = _other[index]
            region = region1 * region2
            _regions1d.append( region )

        return( _regions1d )

    def copyToCommonRegions(self, other, epsilon = domainEpsilon):
        """
        Returns two regions instances that are copies of self and other but with regions added as
        needed so that each has the same number of regions and the region boundaries align.
        """

        if self.dimension != other.dimension:
            raise ValueError('self.dimension = %s not equal to other.dimension = %s' % (self.dimension, other.dimension))

        self2 = self.copy()
        other2 = other.copy()

        if len(self2) == 0:
            if len(other2) == 0:
               return self2, other2

            data = [[other2.domainMin, 0.0], [other2.domainMax, 0.0]]
            xys1d = self2.toLinearXYsClass()(data=data, axes=self2.axes)
            if not isinstance(self2, Regions1d):
                self2 = Regions1d(axes=self2.axes)
            self2.append(xys1d)
        elif len(other2) == 0:
            data = [[self2.domainMin, 0.0], [self2.domainMax, 0.0]]
            xys1d = other2.toLinearXYsClass()(data=data, axes=other2.axes)
            if not isinstance(other2, Regions1d):
                other2 = Regions1d(axes=other2.axes)
            other2.append(xys1d)

        if isinstance(other2, XYs1dModule.XYs1d):
            temp = self.__class__()
            temp.append(other2)
            other2 = temp 
        elif not isinstance(other2, Regions):
            raise NotImplementedError('object of instance "%s" not implemented' % other2.__class__)

        if   self2.domainMin < other2.domainMin:
            region1 = other2[0].copy()
            region1.setData([[self2.domainMin, 0], [other2.domainMin, 0]])
            other2.prepend(region1)
        elif self2.domainMin > other2.domainMin:
            region1 = self2[0].copy()
            region1.setData([[other2.domainMin, 0], [self2.domainMin, 0]])
            self2.prepend(region1)

        if   self2.domainMax > other2.domainMax:
            region1 = other2[-1].copy()
            region1.setData([[other2.domainMax, 0], [self2.domainMax, 0]])
            other2.append(region1)
        elif self2.domainMax < other2.domainMax:
            region1 = self2[-1].copy()
            region1.setData([[self2.domainMax, 0], [other2.domainMax, 0]])
            self2.append(region1)

        boundaries = set()
        for region in self2[1:]:
            boundaries.add(region.domainMin)
        for region in other2[1:]:
            boundaries.add(region.domainMin)
        boundaries = sorted(boundaries)

        count = 0
        priorBoundary = None
        boundariesToMove = []
        for boundary in boundaries:
            if priorBoundary is not None:
                if (boundary - priorBoundary) < epsilon * max(abs(boundary), abs(priorBoundary)):
                    boundariesToMove.append([boundary, priorBoundary])
                    boundary = priorBoundary
                    count += 1
                    if count > 1:
                        raise ValueError('more than one boundary within epsilon = %s of %s' % (epsilon, priorBoundary))
                else:
                    count = 0
            priorBoundary = boundary

        for boundary, priorBoundary in boundariesToMove:
            boundaries.remove(boundary)

        def processBoundaries(regions_, boundaries, boundariesToMove):

            for region in regions_:
                domainMin, domainMax = region.domainMin, region.domainMax
                for boundary, priorBoundary in boundariesToMove:
                    if   domainMin == boundary:
                        region.tweakDomain(domainMin = priorBoundary, epsilon = epsilon)
                    elif domainMax == boundary:
                        region.tweakDomain(domainMax = priorBoundary, epsilon = epsilon)

            for boundary in boundaries:
                for region in regions_:
                    if region.domainMin < boundary < region.domainMax:
                        regions_.splitInTwo(boundary, epsilon = epsilon)
                        break

        if len(boundaries) > 0:
            processBoundaries(self2, boundaries, boundariesToMove)
            processBoundaries(other2, boundaries, boundariesToMove)

        return self2, other2

    def copyDataToXsAndYs(self):
        '''Copies data in regions to a list of xs and ys. Only works if each has  "copyDataToXsAndYs" method.'''

        xs = []
        ys = []
        for region in self:
            subXs, subYs = region.copyDataToXsAndYs()
            xs += subXs
            ys += subYs
        return xs, ys

    def normalize( self, insitu = False, dimension = 1 ) :
        """
        The dimension argument is currently ignored, but kept to be compatible with calling from XYsnd.normalize.
        """

        factor = 1.0 / self.integrate()
        if( insitu ) :
            copy = self
        else :
            copy = self.copy()

        for region in copy : region *= factor
        return( copy )

    def plot(self, **kwargs):
        '''
        Calls multiPlot with only *self* as a curve.
        '''

        self.multiPlot([self], **kwargs)

    def toPointwiseLinear( self, **kwargs ) :
        """
        Returns a new instance, converted to lin-lin interpolation with added points to maintain desired accuracy.
        For more details see toPointwise_withLinearXYs.
        """

        return( self.toPointwise_withLinearXYs( **kwargs ) )

    def toPointwise_withLinearXYs( self, **kwargs ) :
        """
        Converts the regions of self into a single ``XYs1d.XYs1d`` instance that has 'lin-lin' interpolation. At the
        boundary between two abutting regions, the x-values are the same, which is not allowed for an ``XYs1d.XYs1d`` instance.

        Optional (key-word) arguments:
        :param accuracy: indicates desired accuracy. Controls how many points are added when switching interpolation
        :param lowerEps:
        :param upperEps: These arguments are used to smear the x-values at a boundary as follows. Let :math:`(x_l, y_l)` and
        :math:`(x_u, y_u)` be the abutting points for two abutting regions. If :math:`y_l = y_u` then the point :math:`(x_u, y_u)` is removed.
        Otherwise, if( lowerEps > 0 ) the point :math:`(x_l, y_l)` is moved to :math:`x = x_l * ( 1 - lowerEps )` (or :math:`x = x_l * ( 1 + lowerEps )`
        if :math:`x_l < 0`) and the :math:`y` value is interpolated at :math:`x`. If :math:`x` is less than the x-value of the point below :math:`(x_l, y_l)`
        and ``removeOverAdjustedPoints`` is True then the point :math:`(x_l, y_l)` is removed; otherwise, a raise is executed. Similarly
        for upperEps and the point :math:`(x_u, y_u)`.
        :param cls: class to return. Defaults to xData.regions.Regions1d
        """

        def getAdjustedX( x, eps ) :

            if( x == 0. ) :
                x_ = eps
            elif( x < 0 ) :
                x_ = x * ( 1. - eps )
            else :
                x_ = x * ( 1. + eps )
            return( x_ )

        pointwiseClass = self.toLinearXYsClass( )
        if( len( self.regions ) == 0 ) : return( pointwiseClass( data = [], axes = self.axes ) )

        arguments = self.getArguments( kwargs, { 'accuracy' : XYs1dModule.defaultAccuracy, 'lowerEps' : 0, 'upperEps' : 0, 
                'removeOverAdjustedPoints' : False, 'axes' : None } )
        accuracy = arguments['accuracy']
        lowerEps = arguments['lowerEps']
        upperEps = arguments['upperEps']
        removeOverAdjustedPoints = arguments['removeOverAdjustedPoints']
        axes = arguments['axes']

        if( lowerEps < 0. ) : raise ValueError( 'lowerEps = %s must >= 0.' % lowerEps )
        if( upperEps < 0. ) : raise ValueError( 'upperEps = %s must >= 0.' % upperEps )
        if( ( lowerEps == 0. ) and ( upperEps == 0. ) ) : raise ValueError( 'lowerEps and upperEps cannot both be 0.' )

        xys = []
        for iRegion, region in enumerate( self.regions ) :
            _region = region.changeInterpolation(enumsModule.Interpolation.linlin, accuracy, lowerEps = 2 * lowerEps, upperEps = 2 * upperEps )
            _region = _region.copyDataToXYs( )
            if len(_region) < 2: continue
            if iRegion > 0 and len(xys) > 0:
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
                                raise ValueError( 'Adjustment at %s makes new x = %s >= prior x = %s; eps = %s' % ( x12, x, x11, lowerEps ) )
                        else :
                            xys[-1] = [x, XYs1dModule.pointwiseXY_C.interpolatePoint(str(enumsModule.Interpolation.linlin), x, x11, y11, x12, y12)]
                    if( upperEps != 0. ) :
                        x22, y22 = _region[1]
                        x = getAdjustedX( x21, upperEps )
                        if( x >= x22 ) :
                            if( removeOverAdjustedPoints ) :
                                del _region[0]
                            else :
                                raise ValueError( 'Adjustment at %s makes new x = %s >= next x = %s; eps = %s' % ( x21, x, x22, upperEps ) )
                        else :
                            _region[0] = [x, XYs1dModule.pointwiseXY_C.interpolatePoint(str(enumsModule.Interpolation.linlin), x, x21, y21, x22, y22)]
            xys += _region
        pointwise = pointwiseClass( data = xys, axes = self.axes, outerDomainValue = self.outerDomainValue ) # FIXME - need more work to insure all parameters are set properly.
        return( pointwise )

    def toLinearXYsClass( self ) :

        return( XYs1dModule.XYs1d )

    @staticmethod
    def allowedSubElements( ) :

        return( ( XYs1dModule.XYs1d, series1dModule.Series1d ) )

    @staticmethod
    def multiPlot(curve1ds, **kwargs):
        '''
        Plots a list of 1d curves on the same plot. Uses each curve's 'plotLabel' as the legend key. Each curve must have a copyDataToXsAndYs method.
        '''

        XYs1dModule.XYs1d.multiPlot(curve1ds, **kwargs)

class RegionsMultiD( Regions ) :

    def getBoundingSubFunctions( self, domainValue ) :

        for region in self.regions :
            if( domainValue < region.domainMax ) : return( region.getBoundingSubFunctions( domainValue ) )

        return( self.regions[-1].getBoundingSubFunctions( domainValue ) )  # Domain value is above the last region. Let the last region determine what to do.

class Regions2d( RegionsMultiD ) :

    moniker = 'regions2d'
    dimension = 2

    def toLinearXYsClass( self ) :

        return( Regions )

    @staticmethod
    def allowedSubElements( ) :

        from . import multiD_XYs as multiD_XYsModule

        return( ( multiD_XYsModule.XYs2d, ) )

class Regions3d( RegionsMultiD ) :

    moniker = 'regions3d'
    dimension = 3

    def toLinearXYsClass( self ) :

        return( Regions )

    @staticmethod
    def allowedSubElements( ) :

        from . import multiD_XYs as multiD_XYsModule

        return( ( multiD_XYsModule.XYs3d, ) )
