# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module defines the Ys1d class, used to store union-grid functions such as gridded cross sections.
When multiple functions can be mapped onto the same independent axis grid, the Ys1d class saves space by
linking to the common independent grid rather than repeating it.
Otherwise, Ys1d functions should behave like XYs1d.
"""

from pqu import PQU as PQUModule

from . import enums as enumsModule
from . import base as baseModule
from . import axes as axesModule
from . import values as valuesModule
from . import XYs1d as XYs1dModule


class Ys1d(baseModule.XDataFunctional):

    moniker = 'Ys1d'
    dimension = 1

    ancestryMembers = ('Ys',)

    def __init__(self, Ys=None, interpolation=enumsModule.Interpolation.linlin, axes=None,
                 index=None, valueType=enumsModule.ValueType.float64, outerDomainValue=None, label=None):

        baseModule.XDataFunctional.__init__(self, axes, index=index, valueType=valueType,
                                            outerDomainValue=outerDomainValue, label=label)

        self.interpolation = enumsModule.Interpolation.checkEnumOrString(interpolation)

        if Ys is None:
            Ys = valuesModule.Values([])
        if not isinstance(Ys, valuesModule.Values): raise TypeError('Ys must be an instance of values.values.')
        self.__Ys = Ys
        self.__Ys.setAncestor(self)

    def __len__(self):
        # FIXME should this account for self.__Ys.start?

        return len(self.__Ys)

    def __getitem__(self, index):
        # FIXME should this account for self.__Ys.start?

        return self.__Ys[index]

    def __add__(self, other):
        """
        Add two Ys1d instances, checking to ensure they are on the same x grid.
        @param other: must be Ys1d instance or 0
        @return: new Ys1d instance
        """

        if len(self.__Ys) == 0: return other.copy()

        assert isinstance(other, Ys1d), "Adding Ys1d to %s not supported" % type(other)
        if self.grid != other.grid:
            raise NotImplementedError("Adding Ys1d with different grids")
        if self.Ys.length != other.Ys.length:
            raise Exception('self.Ys.length = %d != other.Ys.length = %d' % (self.Ys.length, other.Ys.length))

        if self.__Ys.start <= other.Ys.start:
            ys1d_1 = self.copy()
            ys1d_2 = other
        else:
            ys1d_1 = other.copy()
            ys1d_2 = self 

        offset = ys1d_2.Ys.start - ys1d_1.Ys.start
        values = [y for y in ys1d_1.Ys.values]
        for i1, y2 in enumerate(ys1d_2.Ys): values[i1+offset] += y2
        ys1d_1.Ys.values = values

        return ys1d_1

    @property
    def grid(self):
        """
        Access the X-grid (following link if necessary). Raises AttributeError if no grid found.
        @return: axesModule.Grid
        """

        if isinstance(self.axes[1], axesModule.Grid):
            return self.axes[1]
        elif hasattr(self.axes[1], 'link') and isinstance(self.axes[1].link, axesModule.Grid):
            return self.axes[1].link
        else:
            raise AttributeError("No Grid found for Ys1d!")

    @property
    def Ys(self):
        """
        Access the Values container storing y-values
        @return: xData.values.Values instance
        """

        return self.__Ys

    def convertUnits(self, unitMap):
        """
        unitMap is a dictionary of the form { 'eV' : 'MeV', 'b' : 'mb' }.
        """

        factors = self.axes.convertUnits(unitMap)
        yFactor = factors[0]
        if yFactor != 1:
            self.__Ys.values = [yFactor * value for value in self.__Ys]
        self.fixValuePerUnitChange(factors)

    def copy(self):

        axes = self.axes
        if axes is not None: axes = axes.copy()
        Ys = self.__class__(Ys=self.__Ys.copy(), interpolation=self.interpolation, axes=axes,
                            index=self.index, outerDomainValue=self.outerDomainValue, label=self.label)

        return Ys

    __copy__ = copy

    def evaluate(self, domainValue, extrapolation=enumsModule.Extrapolation.none, epsilon=0):

        raise NotImplementedError("Still TBD. As a temporary work-around, call toPointwise_withLinearXYs to get XYs1d")

    @property
    def domainMin(self):

        return self.grid.domainMin

    @property
    def domainMax(self):

        return self.grid.domainMax

    @property
    def domainUnit(self):

        return self.grid.domainUnit

    def domainUnitConversionFactor(self, unitTo):

        return self.grid.domainUnitConversionFactor(unitTo)

    @property
    def domainGrid(self):

        return self.grid.domainGrid

    @property
    def rangeMin(self):

        return min(self.__Ys.values)

    @property
    def rangeMax(self):

        return max(self.__Ys.values)

    @property
    def rangeUnit(self):

        return self.getAxisUnitSafely(0)

    def rangeUnitConversionFactor(self, unitTo):

        if unitTo is None: return 1.
        return PQUModule.PQU('1 ' + self.rangeUnit).getValueAs(unitTo)

    def toPointwise_withLinearXYs(self, **kwargs):

        cls = kwargs.pop('cls', XYs1dModule.XYs1d)

        xys = [self.domainGrid[self.__Ys.start:self.__Ys.end], self.__Ys]
        return cls(data=xys, dataForm='xsandys', interpolation=self.interpolation, axes=self.axes, index=self.index,
                   valueType=self.valueType, outerDomainValue=self.outerDomainValue, label=self.label)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        outline = kwargs.get('outline', False)
        if len(self) < 6: outline = False

        attributeStr = baseModule.XDataFunctional.attributesToXMLAttributeStr(self)
        if self.interpolation != enumsModule.Interpolation.linlin:
            attributeStr += ' interpolation="%s"' % self.interpolation

        XML_strList = ['%s<%s%s>' % (indent, self.moniker, attributeStr)]
        if self.isPrimaryXData() and self.axes is not None: XML_strList += self.axes.toXML_strList(indent2)
        XML_strList += self.__Ys.toXML_strList(indent2, **kwargs)
        XML_strList[-1] += '</%s>' % self.moniker

        return XML_strList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Translates XML Ys1d into a Ys1d instance.
        """

        attributes, extraAttributes = baseModule.XDataFunctional.parseBareNodeCommonAttributes(node, xPath, True) # parseBareNodeCommonAttributes adds to xPath.
        if len(extraAttributes) > 0: raise Exception('Invalid attributes: %s.' % (', '.join(list(extraAttributes.keys()))))

        values = None
        for child in node:
            if child.tag == 'values': values = valuesModule.Values.parseNodeUsingClass(child, xPath, linkData, **kwargs)
        if values is None: raise Exception('No "values" child node found.')

        ys1d = cls(Ys=values, **attributes)

        extraNodes = baseModule.XDataFunctional.parseNodeStandardChildren(ys1d, node, xPath, linkData, **kwargs)
        if len(extraNodes) == 1: values = extraNodes.pop()

        if len(extraNodes) > 0: raise Exception('Invalid nodes: %s.' % (', '.join([extraNode.tag for extraNode in extraNodes])))

        xPath.pop()                             # Per comment above, parseBareNodeCommonAttributes adds to xPath.

        return ys1d

    @staticmethod
    def defaultAxes(labelsUnits=None):
        """
        :param labelsUnits: dictionary of form { 0 : ( 'dependent label',   'dependent unit' ),
                                                 1 : ( 'independent label', 'independent unit' ) }
        :return: new axes instance
        """

        return axesModule.Axes(2, labelsUnits=labelsUnits)
