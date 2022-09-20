# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import string

from LUPY import ancestry as ancestryModule
from fudge import GNDS_formatVersion as GNDS_formatVersionModule
from pqu import PQU as PQUModule

from . import enums as enumsModule
from . import link as linkModule
from . import values as valuesModule

# only used for GNDS-1.10 and older:
linkGridToken = 'link'

class Axis(ancestryModule.AncestryIO):

    moniker = 'axis'

    def __init__(self, label, index, unit):
        """
        Constructor for the axis class.
        """

        ancestryModule.AncestryIO.__init__(self)

        if not isinstance(label, str): raise TypeError('label = "%s" is not a string' % label)
        self.__label = label

        self.index = index

        self.unit = unit

    def __str__(self):
        """Returns a simple string representation of *self*."""

        return 'label="%s", index="%s", unit="%s"' % ( self.label, self.__index, self.unit )

    def __eq__( self, other ) :

        return isinstance(other, Axis) and self.label == other.label and self.unit == other.unit

    def __ne__(self, other):

        return not self.__eq__(other)

    @property
    def keyName(self):

        return('index')

    @property
    def keyValue(self):

        return(self.index)

    @property
    def index(self):

        return(self.__index)

    @index.setter
    def index(self, value):

        self.__index = int(value)

    @property
    def label(self):

        return(self.__label)

    @label.setter
    def label(self, value):
        """Set the label to *value*."""

        self.__label = value

    @property
    def unit(self):

        return self.__unit

    @unit.setter
    def unit(self, value):
        """Sets self's unit. Only checks that unit is a string. If unit is None, it is set to an empty string (i.e., '')."""

        if value is None: value = ''
        if not isinstance(value, str): raise TypeError('unit type "%s" is not a string' % type(value))
        self.__unit = value.strip()

    def convertUnits(self, unitMap):

        unit, factor = PQUModule.convertUnits(self.unit, unitMap)
        self.unit = unit
        return factor

    def copy(self):
        """Returns a new instance that is a copy of self."""

        return Axis(self.label, self.index, self.unit)

    __copy__ = copy

    def divideUnit(self, other):
        """
        Returns the unit obtained by the division of self.unit by other.unit. Other must be an axis based instance.
        """

        pqu = PQUModule.PQU(1, self.unit) / PQUModule.PQU(1, other.unit)

        return str(pqu.unit)

    def multiplyUnit(self, other):
        """
        Returns the unit obtained by the product of self.unit times other.unit. Other must be an axis based instance.
        """

        pqu = PQUModule.PQU(1, self.unit) * PQUModule.PQU(1, other.unit)

        return str(pqu.unit)

    def plotLabel(self):

        label = self.label
        if label == '': label = 'unknown'
        if self.unit != '': label += ' [%s]' % self.unit

        return label

    def toXML_strList( self, indent = '', **kwargs ) :

        return [ '%s<%s index="%d" label="%s" unit="%s"/>' % ( indent, self.moniker, self.index, self.label, self.unit ) ]

    def unitConversionFactor(self, newUnit):
        """Returns as a float the factor needed to convert self's unit to newUnit. If units are not compatible, a raise is executed."""

        return PQUModule.PQU(1., self.unit).getValueAs(newUnit)

    def parseNode(self, node, xPath, linkData, **kwargs):

        xPath.append( '%s[@index="%s"]' % ( Axis.moniker, node.get( 'index' ) ) )

        self.label = node.get('label')
        self.index = node.get('index')
        self.unit = node.get('unit')

        xPath.pop()

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        axis1 = cls('', 0, '')
        axis1.parseNode(node, xPath, linkData, **kwargs)

        return axis1

class Grid(Axis):

    moniker = 'grid'
    ancestryMembers = ( 'values', )

    def __init__(self, label, index, unit, style, values, interpolation=enumsModule.Interpolation.linlin):
        """
        Returns a new instance of grid.
        """

        Axis.__init__(self, label, index, unit)

        self.__style = enumsModule.GridStyle.checkEnumOrString(style)
        if self.__style == enumsModule.GridStyle.none:                 # Required for GNDS-1.10 support.
            if not isinstance(values, linkModule.Link):
                raise ValueError('style = %s not supported' % style)

        if not isinstance(values, (valuesModule.Values, linkModule.Link)):
            raise TypeError('Unsupported grid values: %s' % type(values))
        self.__values = values
        self.values.setAncestor( self )

        self.interpolation = enumsModule.Interpolation.checkEnumOrString(interpolation)

    @property
    def style(self):

        if self.__style is None and self.isLink():
            # follow link to determine actual style:
            self.__style = self.values.link.ancestor.style
        return self.__style

    @property
    def values(self):

        return self.__values

    def isLink(self):

        return isinstance(self.__values, linkModule.Link)

    @property
    def domainMin(self):

        return self.values[0]

    @property
    def domainMax(self):

        return self.values[-1]

    @property
    def domainUnit(self) :

        return self.unit

    def domainUnitConversionFactor(self, unitTo):

        return self.unitConversionFactor(unitTo)

    @property
    def domainGrid(self):

        return [ value for value in self.values ]

    def convertToUnit(self, unit):

        factor = self.unitConversionFactor(unit)
        self.unit = unit
        if not self.isLink():
            self.__values = valuesModule.Values([ factor * value for value in self.values ])

    def convertUnits(self, unitMap):

        factor = Axis.convertUnits(self, unitMap)
        if factor != 1:
            if not self.isLink():
                self.__values.offsetScaleValues(0, factor)
        return factor

    def copy(self):
        """Returns a new grid instance that is a copy of self."""

        grid1 = Grid(self.label, self.index, self.unit, self.style, self.values.copy(), interpolation=self.interpolation)
        return grid1

    __copy__ = copy

    def getIndexOfValue(self, v):
        """
        Get the index of the value in values where x would fit
        :param v:

        :return:
        """

        for ival, val in enumerate(self.values[:-1]):
            if v >= val and v <= self.values[ival+1]: return ival

        return None

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        attributeStr = ' style="%s"' % self.style
        formatVersion = kwargs.get('formatVersion', GNDS_formatVersionModule.default)
        if self.isLink() and formatVersion == GNDS_formatVersionModule.version_1_10:
            attributeStr = ' style="link"'
        if self.interpolation is not enumsModule.Interpolation.linlin: attributeStr += ' interpolation="%s"' % self.interpolation
        XML_strList = [ '%s<%s index="%d" label="%s" unit="%s"%s>' % ( indent, self.moniker, self.index, self.label, self.unit, attributeStr ) ]
        XML_strList += self.values.toXML_strList(indent2, **kwargs)
        XML_strList[-1] += '</%s>' % self.moniker

        return XML_strList

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        xPath.append('%s[@index="%s"]' % ( Grid.moniker, node.get('index') ))

        for key in list( node.keys( ) ) :
            if( 'href' == key[-4:] ) :
                xPath.pop( )
                return linkModule.Link2.parseNodeUsingClass(node, xPath, linkData, **kwargs)

        style = node.get('style')
        if style == linkGridToken and linkData['formatVersion'] in (
                    GNDS_formatVersionModule.version_1_10,
                    GNDS_formatVersionModule.version_2_0_LLNL_3,
                    GNDS_formatVersionModule.version_2_0_LLNL_4):
            # style needs to be same as the linked-to grid.
            # Set to None for now, look it up later once links are fixed
            style = enumsModule.GridStyle.none

        gridClass = valuesModule.Values
        if node[0].tag == linkModule.Link.moniker:
            gridClass = linkModule.Link

        gridData = gridClass.parseNodeUsingClass(node[0], xPath, linkData, **kwargs)
        interpolation=node.get('interpolation', enumsModule.Interpolation.linlin)
        grid1 = cls(node.get('label'), int(node.get('index')), node.get('unit'), style, gridData, interpolation)

        xPath.pop()

        return grid1

class Axes(ancestryModule.AncestryIO):

    moniker = 'axes'
    ancestryMembers = ( 'axes', )

    def __init__(self, size=0, labelsUnits={}):
        """
        Constructor for ``axes`` class. For example::

            _axes = Axes(labelsUnits = { 0: ( 'crossSection' , 'b' ), 1: ( 'energy_in', 'eV' ) })
        """

        ancestryModule.AncestryIO.__init__(self)

        size = max(size, len(labelsUnits))
        self.axes = []                                  # FIXME2, self.axes needs to be a suite instance.
        if size <= 0: return

        if not size < 26: raise Exception('Size = %d must be less than 26.' % size)

        abcsOffset = string.ascii_lowercase.index('y')
        for index in range(size):
            label, unit = string.ascii_lowercase[abcsOffset-index], ''
            if index in labelsUnits: label, unit = labelsUnits[index]
            self.axes.append(Axis(label, index, unit))

    def __eq__(self, other):

        if not isinstance(other, Axes): raise ValueError('Other not an Axes instance')
        if len(self) != len(other): return False
        for index, _axis in enumerate(self.axes) :
            if _axis != other[index]: return False
        return True

    def __ne__(self, other):

        return not self.__eq__(other)

    def __len__(self):

        return len(self.axes)

    def __getitem__(self, index):

        return self.axes[index]

    def __setitem__(self, index, axisOrGrid):

        if not isinstance(axisOrGrid, ( Axis, Grid, linkModule.Link2 )): raise TypeError('axisOrGrid is not an instance of Axis or Grid')

        size = len(self.axes)
        if index < 0: index += size
        if index == size:
            self.axes.append(axisOrGrid)
        else:
            if not 0 <= index < size: raise IndexError("index = %s out of range for self of size %s" % ( index, size ))
            self.axes[index] = axisOrGrid

        axisOrGrid.index = index

        axisOrGrid.setAncestor(self)

    def __str__(self):
        """Returned a simple string representation of each **Axes** of *self*."""

        return '\n'.join([ str(axis) for axis in self ])

    def convertUnits(self, unitMap):
        """
        Converts each axis units.
        unitMap is a dictionary of mapping old units to new units (e.g., { 'eV' : 'MeV', 'b' : 'mb' }).
        """

        factors = []
        for axis in self :
            if isinstance(axis, linkModule.Link2): continue
            factors.append(axis.convertUnits(unitMap))
        return factors

    def copy(self):

        newAxes = Axes(-1)
        for index, axis in enumerate(self.axes): newAxes[index] = axis.copy()
        for index, child in enumerate(newAxes.axes):
            if isinstance(child, Grid):
                if isinstance(child.values, linkModule.Link): child.values.link = child.values.follow(child.values)

        return newAxes

    __copy__ = copy

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XML_strList = [ '%s<%s>' % ( indent, self.moniker ) ]
        for index in reversed(range(len(self))):
            XML_strList += self[index].toXML_strList(indent=indent2, **kwargs)
        XML_strList [-1] += '</%s>' % self.moniker
        return XML_strList

    def parseNode(self, node, xPath, linkData, **kwargs):       # FIXME2, needed until self.axes is a suite instance.

        xPath.append( node.tag )

        self.axes = len(node) * [None]
        for child in node:
            childClass = { Axis.moniker: Axis, Grid.moniker: Grid }.get(child.tag)
            if childClass is None: raise TypeError("Unexpected child node '%s' encountered in axes" % child.tag)
            childAxis = childClass.parseNodeUsingClass(child, xPath, linkData, **kwargs)
            self[childAxis.index] = childAxis

        xPath.pop()

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):

        axes = cls(1)
        axes.parseNode(node, xPath, linkData, **kwargs)

        return axes
