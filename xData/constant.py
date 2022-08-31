# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from numericalFunctions import pointwiseXY_C as pointwiseXY_CModule
floatToShortestString = pointwiseXY_CModule.floatToShortestString

from pqu import PQU as PQUModule

from . import enums as enumsModule
from . import base as baseModule
from . import XYs1d as XYs1dModule

class Constant( baseModule.XDataFunctional )  :

    def __init__(self, _value, domainMin, domainMax, axes=None, label=None, index=None, outerDomainValue=None):

        baseModule.XDataFunctional.__init__(self, label=label, axes=axes)

        self.value = _value

        if( isinstance( domainMin, int ) ) : domainMin = float( domainMin )
        if( not( isinstance( domainMin, float ) ) ) : TypeError( 'domainMin not a float instance' )
        self.__domainMin = domainMin

        if( isinstance( domainMax, int ) ) : domainMax = float( domainMax )
        if( not( isinstance( domainMax, float ) ) ) : TypeError( 'domainMax not a float instance' )
        self.__domainMax = domainMax

    def copy( self ) :

        axes = self.axes
        if( axes is not None ) : axes = self.axes.copy( )
        return( self.__class__( self.value, self.domainMin, self.domainMax, axes = axes, label = self.label ) )

    __copy__ = copy

    @property
    def value( self ) :

        return( self.__value )

    @value.setter
    def value( self, _value ) :

        if( isinstance( _value, int ) ) : _value = float( _value )
        if( not( isinstance( _value, float ) ) ) : TypeError( 'value not a float instance' )
        self.__value = _value

    @property
    def domainMin( self ) :

        return( self.__domainMin )

    @domainMin.setter
    def domainMin(self, domainMin):
        '''Sets self's domainMin to *domainMin*. *domainMin* must be less than self.domainMax or a **raise** is executed.'''

        if domainMin >= self.domainMax:
            raise ValueError('domainMin = %.17 >= self.domainMax = %.17e' % (domainMin, self.domainMax))
        self.__domainMin = domainMin

    @property
    def domainMax( self ) :

        return( self.__domainMax )

    @domainMax.setter
    def domainMax(self, domainMax):
        '''Sets self's domainMax to *domainMax*. *domainMax* must be greater than self.domainMin or a **raise** is executed.'''

        if domainMax <= self.domainMin:
            raise ValueError('domainMax = %.17 <= self.domainMin = %.17e' % (domainMax, self.domainMin))
        self.__domainMax = domainMax

    @property
    def domainGrid(self):

        return self.__domainMin, self.__domainMax

    @property
    def domainUnit( self ) :

        return( self.getAxisUnitSafely( self.dimension ) )

    @property
    def rangeMin( self ) :

        return( self.__value )

    rangeMax = rangeMin

    @property
    def rangeUnit( self ) :

        return( self.getAxisUnitSafely( 0 ) )

    def fixDomainPerUnitChange( self, factors ) :

        self.__domainMin *= factors[self.dimension]
        self.__domainMax *= factors[self.dimension]

    def fixDomains(self, domainMin, domainMax, fixToDomain):
        """
        Sets *domainMin* and *domainMax* per the arguments.
        """

        OldDomainMin = self.domainMin
        OldDomainMax = self.domainMax

        domainMin = max(domainMin, self.domainMin)
        domainMax = min(domainMax, self.domainMax)
        if fixToDomain == enumsModule.FixDomain.lower:
            self.__domainMin = domainMin
        elif fixToDomain == enumsModule.FixDomain.upper:
            self.__domainMax = domainMax
        else:
            self.__domainMin = domainMin
            self.__domainMax = domainMax

        if OldDomainMin == self.domainMin and OldDomainMax == self.domainMax: return 0
        return 1

class Constant1d( Constant ) :

    moniker = 'constant1d'
    dimension = 1

    def __truediv__(self, other):

        if isinstance(other, (int, float)):
            return Constant1d(self.value / other, self.domainMin, self.domainMax,
                              self.axes.copy(), self.label, self.index, self.outerDomainValue)
        else:
            raise NotImplementedError(f"Dividing Constant1d by {type(other)}")

    def __rtruediv__(self, other):

        if isinstance(other, (int, float)):
            return Constant1d(other / self.value, self.domainMin, self.domainMax,
                              self.axes.copy(), self.label, self.index, self.outerDomainValue)
        else:
            raise NotImplementedError(f"Dividing Constant1d by {type(other)}")

    def __itruediv__(self, other):

        if isinstance(other, (int, float)):
            self.value /= other
            return self
        else:
            raise NotImplementedError(f"Dividing Constant1d by {type(other)}")

    def __mul__(self, other):

        if isinstance(other, (int, float)):
            return Constant1d(self.value * other, self.domainMin, self.domainMax,
                              self.axes.copy(), self.label, self.index, self.outerDomainValue)
        else:
            raise NotImplementedError(f"Multiplying Constant1d by {type(other)}")

    __rmul__ = __mul__

    def __imul__(self, other):

        if isinstance(other, (int, float)):
            self.value *= other
            return self
        else:
            raise NotImplementedError(f"Multiplying Constant1d by {type(other)}")

    def convertUnits( self, unitMap ) :
        """
        unitMap is a dictionary of the for { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        if len(self.axes) == 0: return

        factors = self.axes.convertUnits( unitMap )
        self.value *= factors[0]
        self.fixDomainPerUnitChange( factors )
        self.fixValuePerUnitChange( factors )

    def convertAxisToUnit( self, indexOrName, newUnit ) :

        if len(self.axes) == 0:
            n = self
        else:
            n = self.copy()
            index = n.getAxisIndexByIndexOrName( indexOrName )
            axis = n.axes[index]
            factor = PQUModule.PQU( '1 ' + axis.unit ).getValueAs( newUnit )
            if index == 0:
                n.fixDomainPerUnitChange({0: factor})
            elif index == 1:
                n.value *= factor
            axis.unit = newUnit
        return n

    def evaluate( self, x ) :

        return( self.value )

    def toPointwise_withLinearXYs(self, **kwargs):

        xys1d = XYs1dModule.XYs1d(data=[[self.domainMin, self.value], [self.domainMax, self.value]],
                      axes=self.axes, outerDomainValue=self.outerDomainValue, label=self.label)
        return xys1d

    def toXML_strList(self, indent = '', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        valueFormatter = kwargs.get('valueFormatter', floatToShortestString)
        significantDigits = kwargs.get('significantDigits', 15)

        attributeStr = baseModule.XDataCoreMembers.attributesToXMLAttributeStr(self)
        startTag = [ '%s<%s%s value="%s" domainMin="%s" domainMax="%s">' % ( indent, self.moniker, attributeStr, 
                valueFormatter(self.value, significantDigits = significantDigits),
                valueFormatter(self.domainMin, significantDigits = significantDigits),
                valueFormatter(self.domainMax, significantDigits = significantDigits) ) ]

        return self.buildXML_strList(indent2, startTag, [], **kwargs)

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        attributes, extraAttributes = baseModule.XDataFunctional.parseBareNodeCommonAttributes(node, xPath)     # parseBareNodeCommonAttributes adds to xPath.

        value = float(extraAttributes.pop('value'))
        domainMin = float(extraAttributes.pop('domainMin'))
        domainMax = float(extraAttributes.pop('domainMax'))

        if len(extraAttributes) > 0: raise Exception('Invalid attributes: %s.' % ( ', '.join(list(extraAttributes.keys())) ))

        instance = cls(value, domainMin, domainMax, **attributes)

        extraNodes = baseModule.XDataFunctional.parseNodeStandardChildren(instance, node, xPath, linkData, **kwargs)
        if len(extraNodes) > 0: raise Exception('Invalid nodes: %s.' % (', '.join([extraNode.tag for extraNode in extraNodes])))

        xPath.pop()                                     # Per comment above, parseBareNodeCommonAttributes adds to xPath.

        return instance
