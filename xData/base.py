# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the xData class that must be the base class for all xData classes.
"""

import abc

from pqu import PQU as PQUModule

from LUPY import ancestry as ancestryModule

from . import enums as enumsModule

def getArguments(kwargs, arguments):

    for argument in arguments: arguments[argument] = kwargs.get(argument, arguments[argument])
    return arguments

class XDataBase:

    pass

class XData(XDataBase, ancestryModule.AncestryIO):

    def __init__(self):

        ancestryModule.AncestryIO.__init__(self)

class XDataCoreMembers(XData):

    def __init__(self, index=None, label=None):

        XData.__init__(self)

        if index is not None: index = int(index)
        self.index = index
        self.label = label

    def attributesToXMLAttributeStr(self):

        attributeStr = ''
        if self.index is not None: attributeStr += ' index="%s"' % self.index
        if self.label is not None: attributeStr += ' label="%s"' % self.label
        return attributeStr

    @property
    def index(self):

        return self.__index

    @index.setter
    def index(self, value):

        if value is not None: value = int(value)
        self.__index = value

    @property
    def label(self):

        return self.__label

    @label.setter
    def label(self, value):

        if value is not None:
            if not isinstance(value, str): raise TypeError('label must be a string, got %s' % type(value))
        self.__label = value

    @staticmethod
    def getArguments(kwargs, arguments):

        return getArguments(kwargs, arguments)

class XDataFunctional(XDataCoreMembers):

    ancestryMembers = ( 'axes', 'uncertainty' )

    def __init__(self, axes, index=None, valueType=enumsModule.ValueType.float64, outerDomainValue=None, label=None):

        from . import axes as axesModule
        from . import uncertainties as uncertaintiesModule

        XDataCoreMembers.__init__(self, index=index, label=label)

        self.valueType = valueType

        if outerDomainValue is not None: outerDomainValue = float(outerDomainValue)
        self.__outerDomainValue = outerDomainValue

        if axes is None: axes = axesModule.Axes(0)
        self.axes = axes
        self.uncertainty = uncertaintiesModule.Uncertainty()                # FIXME2, uncertainty should be a suite like object.
        self.uncertainty.setAncestor(self)

    @property
    @abc.abstractmethod
    def dimension(self):

        pass

    def attributesToXMLAttributeStr(self):

        attributeStr = XDataCoreMembers.attributesToXMLAttributeStr(self)
        if self.outerDomainValue is not None: attributeStr += ' outerDomainValue="%s"' % PQUModule.floatToShortestString(self.outerDomainValue, 12)
        return attributeStr

    @property
    def axes(self):
        """Returns self's axes."""

        return self.__axes

    @axes.setter
    def axes(self, axes1):

        from . import axes as axesModule

        if not isinstance(axes1, axesModule.Axes): raise TypeError('axes in not an axes instance')
        if len(axes1) > 0:
            if len(axes1) <= self.dimension:
                raise Exception('len( axes ) = %d != ( self.dimension + 1 ) = %d' % (len(axes1), self.dimension + 1))
        if axes1.ancestor is not None:
            axes1 = axes1.copy()
        axes1.setAncestor(self)
        self.__axes = axes1

    @property
    def outerDomainValue(self):

        return self.__outerDomainValue

    @outerDomainValue.setter
    def outerDomainValue(self, outerDomainValue):

        if outerDomainValue is not None: outerDomainValue = float(outerDomainValue)
        self.__outerDomainValue = outerDomainValue

    @property
    def valueType(self):

        return self.__valueType

    @valueType.setter
    def valueType(self, valueType):

        valueType = enumsModule.ValueType.checkEnumOrString(valueType)
        if valueType != enumsModule.ValueType.float64:
            raise TypeError('valueType = "%s" not supported' % valueType)
        self.__valueType = valueType

    def fixValuePerUnitChange(self, factors):

        if self.__outerDomainValue is not None: self.__outerDomainValue *= factors[self.dimension+1]

    def getAxisIndexByIndexOrName(self, indexOrName):

        if isinstance(indexOrName, int): return indexOrName
        if isinstance(indexOrName, str):
            for index, axis in enumerate(self.axes):
                if axis.label == indexOrName: return index
        raise TypeError('argument must be an integer or a string')

    def getAxisUnitSafely(self, index):

        if len(self.axes) == 0: return ''
        return self.axes[index].unit

    def getPrimaryXData(self):

        ancestor = self
        while not ancestor.isPrimaryXData(): ancestor = ancestor.ancestor
        return ancestor

    def isPrimaryXData(self):
        """Returns False if self is contained in a higher dimension XDataFunctional and False otherwise."""

        ancestry = self.ancestor
        if ancestry is None: return True
        return not isinstance(ancestry, XDataFunctional)

    def buildXML_strList(self, indent, startTag, extraData, **kwargs):
        """
        This methods builds the list of XML strings by adding (if present) the axes, uncertainty data and the end tag name.
        It is mainly for internal use.
        """

        XML_strList = startTag
        if self.isPrimaryXData() and len(self.axes) > 0: XML_strList += self.axes.toXML_strList(indent=indent, **kwargs)
        XML_strList += extraData
        XML_strList += self.uncertainty.toXML_strList(indent=indent, **kwargs)
        XML_strList[-1] += '</%s>' % self.moniker

        return XML_strList

    def parseNodeStandardChildren(self, node, xPath, linkData, **kwargs):
        """
        Adds 'axes' and 'uncertainty' parsing to *self* and returns the child nodes not parsed as a dictionary with keys of the child.tag values.
        """

        extraNodes = []
        for child in node:
            if child.tag == 'axes':
                self.axes.parseNode(child, xPath, linkData, **kwargs)
            elif child.tag == 'uncertainty':
                self.uncertainty.parseNode(child, xPath, linkData, **kwargs)
            else:
                extraNodes.append(child)

        if len(self.axes) == 0 and 'axes' in kwargs: self.axes = kwargs['axes']

        return extraNodes

    @staticmethod
    def parseBareNodeCommonAttributes(node, xPath, allowInterapolation=False):

        keyExist = False
        for attrName in ( 'label', 'outerDomainValue' ):
            if node.get(attrName) is not None:
                keyExist = True
                xPath.append('%s[@%s="%s"]' % (node.tag, attrName, node.get(attrName) ))
                break
        if not keyExist: xPath.append(node.tag)

        extraAttributes = {}
        attributes =     { 'label': None, 'index': None, 'outerDomainValue': None  }
        attributeTypes = { 'label': str,  'index': int,  'outerDomainValue': float }
        if allowInterapolation:
            attributes['interpolation'] = enumsModule.Interpolation.linlin
            attributeTypes['interpolation'] = str

        for key, item in list(node.items()):
            if key in attributes:
                attributes[key] = attributeTypes[key](item)
            else:
                extraAttributes[key] = item

        return attributes, extraAttributes

def isXDataFunctional(object):
    """Returns True if object is an instance of XDataBase and False otherwise."""

    return isinstance(object, XDataBase)

def getDomainValue(value, unit, default):

    if value is None: return default
    if isinstance(value, ( str, PQUModule.PQU )): return PQUModule.PQU(value).getValueAs(unit)
    return value

def getDomainLimits(self, domainMin, domainMax, unit):

    defaultMin, defaultMax = self.domainMin, self.domainMax
    return getDomainValue(domainMin, unit, defaultMin), getDomainValue(domainMax, unit, defaultMax)

def processUnits(unit1, unit2, operator):

    if operator not in [ '*', '/' ]: raise ArithmeticError('unsupported unit operation "%s"' % operator)
    result = eval('PQUModule.PQU(1, unit1) %s PQUModule.PQU(1, unit2)' % operator)
    return result.getUnitSymbol()

def getDomainValue2(domainValue):

    if isinstance(domainValue, PQUModule.PQU): return domainValue
    if isinstance(domainValue, str):
        try:
            return float(domainValue)
        except:
            return PQUModule.PQU(domainValue)
    return float(domainValue)
