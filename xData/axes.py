# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module containes all the classes for handling GNDS axes and its child nodes.

This module contains the following classes:

    +-----------------------------------+-----------------------------------------------------------------------+
    | Class                             | Description                                                           |
    +===================================+=======================================================================+
    | Axis                              | This class represents a GNDS axis node.                               |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Grid                              | This class represents a GNDS grid node.                               |
    +-----------------------------------+-----------------------------------------------------------------------+
    | Axes                              | This class represents a GNDS axes node.                               |
    +-----------------------------------+-----------------------------------------------------------------------+
"""

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
    """
    This class represents a GNDS axis node.

    The following table list the primary members of this class:

    +-----------+---------------------------------------------------------------+
    | Member    | Description                                                   |
    +===========+===============================================================+
    | label     | The label for the axis.                                       |
    +-----------+---------------------------------------------------------------+
    | index     | The index of the axis in the parent :py:class: Axes instance. |
    +-----------+---------------------------------------------------------------+
    | unit      | The unit for the data for the axis.                           |
    +-----------+---------------------------------------------------------------+
    """

    moniker = 'axis'

    def __init__(self, label, index, unit):
        """
        :param label:       The label for the axis.
        :param index:       The index of the axis in the parent :py:class: Axes instance.
        :param unit:        The unit for the data for the axis. 
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
        """
        This method returns True if *self* and *other* are equal and False otherwise.
        For two :py:class:`Axis` instances to be equal, they must have the same *label* and *unit*.

        :param other:   Another :py:class:`Axis` instance to compare with *self*.

        :returns:       A boolean instance.
        """

        return isinstance(other, Axis) and self.label == other.label and self.unit == other.unit

    def __ne__(self, other):
        """
        This method returns True if *self* and *other* are not equal and False otherwise.
        For two :py:class:`Axis` instances to be equal, they must have the same *label* and *unit*.

        :param other:   Another :py:class:`Axis` instance to compare with *self*.

        :returns:       A boolean instance.
        """

        return not self.__eq__(other)

    @property
    def keyName(self):
        """ 
        This method returns the key name for *self*.

        :returns:       A python str instance. 
        """

        return('index')

    @property
    def keyValue(self):
        """
        This method returns the key value for *self*.

        :returns:       Whatever the type of the keyValue is.
        """

        return(self.index)

    @property
    def index(self):
        """
        This method returns the index for *self*.

        :returns:       A python int.
        """

        return(self.__index)

    @index.setter
    def index(self, value):
        """
        This method sets the index for *self* to *value*.

        :param value:       The new index.
        """

        self.__index = int(value)

    @property
    def label(self):
        """
        This method returns the label for *self*.

        :returns:       A python str.
        """

        return(self.__label)

    @label.setter
    def label(self, value):
        """
        This method sets the label for *self* to *value*.

        :param value:       The new label.
        """

        self.__label = value

    @property
    def unit(self):
        """
        This method returns the unit for *self*.

        :returns:       A python str.
        """

        return self.__unit

    @unit.setter
    def unit(self, value):
        """
        This method sets the unit of *self* to *value. This method only checks that unit is a string or None. 
        If None, the unit is set to an empty string (i.e., '').

        :param value:   The new unit.
        """

        if value is None: value = ''
        if not isinstance(value, str): raise TypeError('unit type "%s" is not a string' % type(value))
        self.__unit = value.strip()

    def convertUnits(self, unitMap):
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        unit, factor = PQUModule.convertUnits(self.unit, unitMap)
        self.unit = unit
        return factor

    def copy(self):
        """
        This method returns a new instance that is a copy of *self*.

        :returns:           An instance of :py:class:`Axis`.
        """

        return Axis(self.label, self.index, self.unit)

    __copy__ = copy

    def divideUnit(self, other):
        """
        This method returns the unit obtained by the division of self.unit by other.unit. Other must be an instance of :py:class:`Axis`.

        :param other:   An :py:class:`Axis` instance.

        :returns:       A python str.
        """

        pqu = PQUModule.PQU(1, self.unit) / PQUModule.PQU(1, other.unit)

        return str(pqu.unit)

    def multiplyUnit(self, other):
        """
        This method returns the unit obtained by the product of self.unit times other.unit. Other must be an instance of :py:class:`Axis`.

        :param other:   An :py:class:`Axis` instance.

        :returns:       A python str.
        """

        pqu = PQUModule.PQU(1, self.unit) * PQUModule.PQU(1, other.unit)

        return str(pqu.unit)

    def plotLabel(self):
        """
        This method returns the string composed of self's label and unit that can be used as an axis label.

        :returns:       A python str.
        """

        label = self.label
        if label == '': label = 'unknown'
        if self.unit != '': label += ' [%s]' % self.unit

        return label

    def toXML_strList( self, indent = '', **kwargs ) :
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        return [ '%s<%s index="%d" label="%s" unit="%s"/>' % ( indent, self.moniker, self.index, self.label, self.unit ) ]

    def unitConversionFactor(self, newUnit):
        """
        This method returns a scale factor as a float that is needed to convert self's unit to *newUnit*. If units are not compatible, a raise is executed.

        :param newUnit:     A unit.

        :returns:           A float.
        """

        return PQUModule.PQU(1., self.unit).getValueAs(newUnit)

    def parseNode(self, node, xPath, linkData, **kwargs):
        """
        This method sets data in *self* using the contents of *node*.

        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        """

        xPath.append( '%s[@index="%s"]' % ( Axis.moniker, node.get( 'index' ) ) )

        self.label = node.get('label')
        self.index = node.get('index')
        self.unit = node.get('unit')

        xPath.pop()

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return: an instance of *cls* representing *node*.
        """

        axis1 = cls('', 0, '')
        axis1.parseNode(node, xPath, linkData, **kwargs)

        return axis1

class Grid(Axis):
    """
    This class represents a GNDS grid node.

    The following table list the primary members of this class:

    +-------------------+---------------------------------------------------------------+
    | Member            | Description                                                   |
    +===================+===============================================================+
    | label             | The label for the axis.                                       |
    +-------------------+---------------------------------------------------------------+
    | index             | The index of the axis in the parent :py:class: Axes instance. |
    +-------------------+---------------------------------------------------------------+
    | unit              | The unit for the data for the axis.                           |
    +-------------------+---------------------------------------------------------------+
    | style             | The style of the grid.                                        |
    +-------------------+---------------------------------------------------------------+
    | values            | The list of numbers for the grid.                             |
    +-------------------+---------------------------------------------------------------+
    | interpolaction    | The interpolation rule for the values in the grid.            |
    +-------------------+---------------------------------------------------------------+
    """

    moniker = 'grid'
    ancestryMembers = ( 'values', )

    def __init__(self, label, index, unit, style, values, interpolation=enumsModule.Interpolation.linlin):
        """
        :param label:           The label for the axis.
        :param index:           The index of the axis in the parent :py:class: Axes instance.
        :param unit:            The unit for the data for the axis. 
        :param style:           The style of the grid.
        :param values:          The list of grid values.
        :param interpolation:   The interpolation rule for the values in the grid.
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
        """
        Thid method returns *self*'s style.

        :returns:               A python str.
        """

        if self.__style is None and self.isLink():
            # follow link to determine actual style:
            self.__style = self.values.link.ancestor.style
        return self.__style

    @property
    def values(self):
        """
        Thid method returns *self*'s value instance.

        :returns:               A :py:class:`valuesModule.Values` or :py:class:`linkModule.Link` instance.
        """

        return self.__values

    def isLink(self):
        """
        This method returns True if *self* is a :py:class:`linkModule.Link` instance and False otherwise.

        :returns:               A boolean.
        """

        return isinstance(self.__values, linkModule.Link)

    @property
    def domainMin(self):
        """
        This method returns the minimum domain value for *self*'s grid.

        :returns:       A number.
        """

        return self.values[0]

    @property
    def domainMax(self):
        """
        This method returns the maximum domain value for *self*'s grid.

        :returns:       A number.
        """

        return self.values[-1]

    @property
    def domainUnit(self) :
        """
        This method returns the domain unit for *self*'s grid with is the same as self's unit.

        :returns:       A python str.
        """

        return self.unit

    def domainUnitConversionFactor(self, unitTo):
        """
        This method returns the factor needed to convert self's domain to unit *unitTo*.

        :param unitTo:      The unit for converting self's domain.

        :returns:           A float.
        """

        return self.unitConversionFactor(unitTo)

    @property
    def domainGrid(self):
        """
        This method returns all domain values for *self* as a python list.

        :returns:           A python list.
        """

        return [ value for value in self.values ]

    def convertToUnit(self, unit):
        """
        This method changes *self*'s unit to *unit* and scales the grid to the new unit.

        :param unit:    The new unit.
        """

        factor = self.unitConversionFactor(unit)
        self.unit = unit
        if not self.isLink():
            self.__values = valuesModule.Values([ factor * value for value in self.values ])

    def convertUnits(self, unitMap):
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        factor = Axis.convertUnits(self, unitMap)
        if factor != 1:
            if not self.isLink():
                self.__values.offsetScaleValues(0, factor)
        return factor

    def copy(self):
        """
        Returns a new grid instance that is a copy of self.

        :returns:       An instance of :py:class:`Grid`.
        """

        grid1 = Grid(self.label, self.index, self.unit, self.style, self.values.copy(), interpolation=self.interpolation)
        return grid1

    __copy__ = copy

    def getIndexOfValue(self, v):
        """
        Thie method returns the lower index of the two elements of the grid where *v* is between. If *v* is not in the domain
        of the grid, None is returned.

        :param v:       The value whose index is returned.

        :returns:       A python int or None.
        """

        for ival, val in enumerate(self.values[:-1]):
            if v >= val and v <= self.values[ival+1]: return ival

        return None

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

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
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:            An instance of *cls* representing *node*.
        """

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
    """
    This class represents a GNDS Axes node. Basically, this class stores a list :py:class:`Axis` and/or py:class:`Grid` children.

    The following table list the primary members of this class:

    +-------------------+---------------------------------------------------------------+
    | Member            | Description                                                   |
    +===================+===============================================================+
    | axes              | The list of axis and/or grid children.                        |
    +-------------------+---------------------------------------------------------------+
    """


    moniker = 'axes'
    ancestryMembers = ( 'axes', )

    def __init__(self, size=0, labelsUnits={}):
        """
        Constructor for ``axes`` class. For example::

            axes = Axes(labelsUnits = { 0: ( 'crossSection' , 'b' ), 1: ( 'energy_in', 'eV' ) })

        The *labelsUnits* argument is a dictionary where the keys are the index for each child :py:class:`Axis`
        and the corresponding value is a list of (label, unit). Child :py:class:`Axis` not specified by *labelsUnits*
        are given default labels and an empty unit.

        :param size:            The number of :py:class:`Axis` children to create.
        :param labelsUnits:     A dictionary of initial labels and units for each :py:class:`Axis` children.
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
        """
        This method returns True if *self* and *other* have the same child axis nodes.

        :param other:   Another :py:class:`Axes` instance to compare with *self*.

        :returns:       A boolean instance.
        """

        if not isinstance(other, Axes): raise ValueError('Other not an Axes instance')
        if len(self) != len(other): return False
        for index, _axis in enumerate(self.axes) :
            if _axis != other[index]: return False
        return True

    def __ne__(self, other):
        """
        This method returns True if *self* and *other* have the same child axis nodes.

        :param other:   Another :py:class:`Axes` instance to compare with *self*.

        :returns:       A boolean instance.
        """

        return not self.__eq__(other)

    def __len__(self):
        """
        This method returns the number of child axis nodes.

        :returns:       A python int.
        """

        return len(self.axes)

    def __getitem__(self, index):
        """
        This method returns the child node at *index*.

        :param index:       The index of the child node to return.

        :returns:           An :py:class:`Axes` or :py:class:`Grid` instance.
        """

        return self.axes[index]

    def __setitem__(self, index, axisOrGrid):
        """
        This method sets the child node at *index*.

        :param index:       The index to where *axisOrGrid* is put.
        :param axisOrGrid:  An :py:class:`Axes` or :py:class:`Grid` instance.
        """

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
        """
        Returned a simple string representation of each **Axes** of *self*.

        :returns:   A python str.
        """

        return '\n'.join([ str(axis) for axis in self ])

    def convertUnits(self, unitMap):
        """
        Converts all data in *self* per *unitMap*.

        :param unitMap:     A dictionary in which each key is a unit that will be replaced by its value which must be an equivalent unit.
        """

        factors = []
        for axis in self :
            if isinstance(axis, linkModule.Link2): continue
            factors.append(axis.convertUnits(unitMap))
        return factors

    def copy(self):
        """
        This method returns a new instance that is a copy of *self*.

        :returns:           An instance of :py:class:`Axes`.
        """

        newAxes = Axes(-1)
        for index, axis in enumerate(self.axes): newAxes[index] = axis.copy()
        for index, child in enumerate(newAxes.axes):
            if isinstance(child, Grid):
                if isinstance(child.values, linkModule.Link): child.values.link = child.values.follow(child.values)

        return newAxes

    __copy__ = copy

    def toXML_strList(self, indent='', **kwargs):
        """
        Returns a list of str instances representing the XML lines of *self*.

        :param indent:          The minimum amount of indentation.
        :param kwargs:          A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:                List of str instances representing the XML lines of self.
        """

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        XML_strList = [ '%s<%s>' % ( indent, self.moniker ) ]
        for index in reversed(range(len(self))):
            XML_strList += self[index].toXML_strList(indent=indent2, **kwargs)
        XML_strList [-1] += '</%s>' % self.moniker
        return XML_strList

    def parseNode(self, node, xPath, linkData, **kwargs):       # FIXME2, needed until self.axes is a suite instance.
        """
        This method sets data in *self* using the contents of *node*.

        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.
        """

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
        """
        Parse *node* into an instance of *cls*.

        :param cls:         Form class to return.
        :param node:        Node to parse.
        :param xPath:       List containing xPath to current node, useful mostly for debugging.
        :param linkData:    dict that collects unresolved links.
        :param kwargs:      A dictionary of extra arguments that controls how *self* is converted to a list of XML strings.

        :return:            An instance of *cls* representing *node*.
        """

        axes = cls(1)
        axes.parseNode(node, xPath, linkData, **kwargs)

        return axes
