# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Defines classes used in both resolved and unresolved regions
"""

import fractions

from LUPY import ancestry as ancestryModule

from fudge import suites as suitesModule
from fudge.outputChannelData import Q as QModule
from pqu import PQU as PQUModule
from .scatteringRadius import ScatteringRadius, HardSphereRadius

from fudge import GNDS_formatVersion as GNDS_formatVersionModule

from xData import table as tableModule
from xData import link as linkModule


class Spin(fractions.Fraction):
    """
    Store spins for a collection of resonances. Check denominator (must be integer or half-integer)
    """

    def __new__(cls, *args):
        self = fractions.Fraction.__new__(cls, *args)
        if self.denominator not in (1, 2):
            raise ValueError("Illegal spin '%s': must be integer or half-integer" % self)
        return self

    @property
    def value(self):
        return float(self)


class Parity:
    """
    Store parity for a collection of resonances. Allowed values for the parity are +1 or -1.
    """

    def __init__(self, parity):
        self.value = int(parity)
        if self.value not in (1, -1):
            raise ValueError("%d is not a legal value for parity!" % self.value)

    def __str__(self):
        return str(self.value)

    def __int__(self):
        return self.value

    def __copy__(self):
        return Parity(self.value)


class ResonanceReactions(suitesModule.Suite):
    """
    Stores a list of resonanceReaction
    """

    moniker = 'resonanceReactions'

    def __init__(self):

        suitesModule.Suite.__init__(self, [ResonanceReaction])

    def check(self, info):
        from fudge import warning
        warnings = []
        for c in self:
            warningList = c.check(info)
            if warningList:
                warnings.append(warning.Context(str(c.moniker) + ' ' + str(c.label), warningList))
        return warnings


class ResonanceReaction(ancestryModule.AncestryIO):
    """
    Describes one reaction channel that opens up in the resonance region. In an R-Matrix section,
    all open reaction channels should be described in the list of resonanceReaction elements
    """

    moniker = 'resonanceReaction'

    fission = 'fission'  # special tokens to support fission reactions
    fissionProduct = 'fissionProduct'

    def __init__(self, label, link, ejectile, Q=None, scatteringRadius=None, hardSphereRadius=None,
                 boundaryConditionValue=None, eliminated=False):
        """
        :param label: unique label for this reaction
        :param link: href pointing to corresponding reaction in reactionSuite/reactions
        :param ejectile: id for the light particle emitted by this reaction
        :param Q: optional, overloads the linked reaction Q-value
        :param scatteringRadius: optional, overrides resonances.scatteringRadius
        :param hardSphereRadius: optional, overrides resonances.scatteringRadius
        :param boundaryConditionValue: optional, set numeric value of boundary condition
        :param eliminated: boolean, default=False
        """

        ancestryModule.AncestryIO.__init__(self)
        self.label = label
        self.__link = link
        self.__link.setAncestor(self)
        self.__ejectile = ejectile
        self.__residual = None
        self.Q = Q
        self.scatteringRadius = scatteringRadius
        self.hardSphereRadius = hardSphereRadius
        self.boundaryConditionValue = boundaryConditionValue
        self.eliminated = eliminated

    @property
    def ejectile(self):
        return self.__ejectile

    @property
    def residual(self):
        if self.__residual is None:
            if self.ejectile == ResonanceReaction.fission: return ResonanceReaction.fissionProduct

            products = set([p.pid for p in self.link.link.outputChannel.products])
            products.remove(self.ejectile)
            if len(products) != 1: raise ValueError("Cannot compute resonanceReaction residual!")
            self.__residual = products.pop()
        return self.__residual

    @property
    def link(self):
        return self.__link

    @property
    def Q(self):
        return self.__Q

    def getQ(self):
        """ Return Q-value, referring to linked reaction if no Q defined locally. """
        if self.__Q is not None:
            return self.__Q
        else:
            reaction = self.link.link
            return reaction.outputChannel.Q

    @Q.setter
    def Q(self, value):
        """Can be set to None or to a Q.Component instance."""
        if value is not None:
            if not isinstance(value, QModule.Component):
                raise TypeError("ResonanceReaction Q value can't be set to type '%s'" % type(value))
            value.setAncestor(self)
        self.__Q = value

    @property
    def scatteringRadius(self):

        return self.__scatteringRadius

    def getScatteringRadius(self):
        """Return ScatteringRadius, looking up ancestry if necessary"""
        if self.__scatteringRadius is not None:
            return self.__scatteringRadius
        else:
            from .resonances import Resonances
            return self.findClassInAncestry(Resonances).getScatteringRadius()

    @scatteringRadius.setter
    def scatteringRadius(self, value):
        """Can be set to None or to a ScatteringRadius instance."""

        if value is not None:
            if not isinstance(value, ScatteringRadius):
                raise TypeError("Scattering radius can't be set to type '%s'" % type(value))
            value.setAncestor(self)
        self.__scatteringRadius = value

    @property
    def hardSphereRadius(self):

        return self.__hardSphereRadius

    def getHardSphereRadius(self):
        """
        Return HardSphereRadius, looking up ancestry if necessary.
        If no HardSphereRadius is defined, return ScatteringRadius instead
        """
        if self.hardSphereRadius is not None:
            return self.hardSphereRadius
        else:
            from .resonances import Resonances
            resonancesHSR = self.findClassInAncestry(Resonances).hardSphereRadius
            if resonancesHSR is not None:
                return resonancesHSR
            else:
                return self.getScatteringRadius()

    @hardSphereRadius.setter
    def hardSphereRadius(self, value):
        """Can be set to None or to a hardSphereRadius instance."""

        if value is not None:
            if not isinstance(value, HardSphereRadius):
                raise TypeError("Hard sphere radius can't be set to type '%s'" % type(value))
            value.setAncestor(self)
        self.__hardSphereRadius = value

    def check(self, info):
        from fudge import warning
        warnings = []

        # check the reaction link
        theLinkTarget = None
        try:
            theLinkTarget = self.link.follow(self.rootAncestor)
            if theLinkTarget is None: warnings.append(warning.UnresolvedLink(self.link))
        except:
            warnings.append(warning.UnresolvedLink(self.link))

        # check the radii
        for thing in [self.scatteringRadius, self.hardSphereRadius]:
            if thing is None: continue
            warningList = thing.check(info)
            if warningList:
                warnings.append(warning.Context(thing.moniker, warningList))
        return warnings

    def convertUnits(self, unitMap):
        for child in ('Q', 'scatteringRadius', 'hardSphereRadius'):
            if getattr(self, child) is not None:
                getattr(self, child).convertUnits(unitMap)

    def isFission(self):

        return self.link.link.isFission()

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + '  '
        attrstring = ''
        if self.ejectile is not None: attrstring += ' ejectile="%s"' % self.ejectile
        if kwargs.get('formatVersion') == GNDS_formatVersionModule.version_1_10:
            if self.eliminated or self.isFission():
                attrstring += ' calculatePenetrability="false"'
        if self.boundaryConditionValue is not None:
            attrstring += ' boundaryConditionValue="%r"' % self.boundaryConditionValue
        if self.eliminated: attrstring += ' eliminated="true"'
        xmlString = ['%s<%s label="%s"%s>' % (indent, self.moniker, self.label, attrstring)]
        xmlString += self.link.toXML_strList(indent=indent2, **kwargs)
        if self.Q is not None: xmlString += self.Q.toXML_strList(indent=indent2, **kwargs)
        if self.__scatteringRadius is not None:
            xmlString += self.__scatteringRadius.toXML_strList(indent=indent2, **kwargs)
        if self.hardSphereRadius is not None:
            xmlString += self.hardSphereRadius.toXML_strList(indent=indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        xPath.append(element.tag)
        reactionLink = linkModule.Link.parseNodeUsingClass(element.find('link'), xPath, linkData, **kwargs)
        Qval, scatRad, hsRad = None, None, None
        if element.find(QModule.Component.moniker):
            Qval = QModule.Component()
            Qval.parseNode(element.find(QModule.Component.moniker), xPath, linkData, **kwargs)
        if element.find(ScatteringRadius.moniker):
            scatRad = ScatteringRadius.parseNodeUsingClass(
                element.find(ScatteringRadius.moniker), xPath, linkData, **kwargs)
        if element.find(HardSphereRadius.moniker):
            hsRad = HardSphereRadius.parseNodeUsingClass(
                element.find(HardSphereRadius.moniker), xPath, linkData, **kwargs)

        boundaryConditionVal = element.get('boundaryConditionValue')
        if boundaryConditionVal is not None:
            boundaryConditionVal = float(boundaryConditionVal)

        tmp = cls(element.get('label'), link=reactionLink, ejectile=element.get('ejectile'),
                  Q=Qval, scatteringRadius=scatRad, hardSphereRadius=hsRad, boundaryConditionValue=boundaryConditionVal,
                  eliminated=getBool(element.get('eliminated', 'false')))
        xPath.pop()
        return tmp


class ResonanceParameters(ancestryModule.AncestryIO):
    """
    Light-weight wrapper around a table.
    """

    moniker = 'resonanceParameters'

    def __init__(self, table):
        ancestryModule.AncestryIO.__init__(self)
        self.table = table
        self.table.setAncestor(self)

    def check(self, info):
        warnings = []
        return warnings

    def convertUnits(self, unitMap):
        self.table.convertUnits(unitMap)

    def toXML_strList(self, indent='', **kwargs):
        indent2 = indent + '  '
        xmlList = ['%s<%s>' % (indent, self.moniker)]
        xmlList += self.table.toXML_strList(indent2, **kwargs)
        xmlList[-1] += '</%s>' % self.moniker
        return xmlList

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):
        xPath.append(element.tag)

        linkData['conversionTable'] = {'L': int}
        rps = cls(tableModule.Table.parseNodeUsingClass(
                  element.find(tableModule.Table.moniker), xPath, linkData, **kwargs))
        del linkData['conversionTable']

        xPath.pop()
        return rps


class EnergyIntervals(ancestryModule.AncestryIO):
    """ Resonance region may be broken up into multiple energy intervals (deprecated) """

    moniker = 'energyIntervals'
    keyName = 'label'

    def __init__(self, label):

        ancestryModule.AncestryIO.__init__(self)
        self.label = label
        self.__intervals = []

    def __len__(self):
        return len(self.__intervals)

    def __getitem__(self, item):
        return self.__intervals[item]

    def append(self, item):
        self.__intervals.append(item)
        item.setAncestor(self)

    def check(self, info):
        from fudge import warning
        warnings = []
        for idx, interval in enumerate(self):
            info['energyIntervalIndex'] = idx
            warningList = interval.check(info)
            if warningList:
                warnings.append(warning.Context('%s[@index="%d"' % (interval.moniker, interval.index), warningList))
        del info['energyIntervalIndex']
        return warnings

    def convertUnits(self, unitMap):
        for interval in self:
            interval.convertUnits(unitMap)

    @property
    def useForSelfShieldingOnly(self):
        for interval in self:
            if interval.evaluated.useForSelfShieldingOnly: return True
        return False

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        xmlString = ['%s<%s label="%s">' % (indent, self.moniker, self.label)]
        for interval in self:
            xmlString += interval.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker

        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        EIs = cls(element.get('label'))
        for child in element:
            EIs.append(EnergyInterval.parseNodeUsingClass(child, xPath, linkData, **kwargs))
        xPath.pop()
        return EIs


class EnergyInterval(ancestryModule.AncestryIO):
    """ single energy interval, for use inside energyIntervals """

    moniker = 'energyInterval'
    keyName = 'index'

    _requiredAttributes = (
        ('index', int),
        ('domainMin', float),
        ('domainMax', float),
        ('domainUnit', str)
    )

    def __init__(self, index, data, domainMin, domainMax, domainUnit):
        ancestryModule.AncestryIO.__init__(self)
        self.index = index
        self.domainMin = domainMin
        self.domainMax = domainMax
        self.domainUnit = domainUnit
        self.evaluated = data
        self.evaluated.setAncestor(self)

    def check(self, info):
        from fudge import warning
        warnings = []
        warningList = self.evaluated.check(info)
        if warningList:
            warnings.append(warning.Context(self.evaluated.moniker, warningList))
        return warnings

    def convertUnits(self, unitMap):

        if self.domainUnit in unitMap:
            newUnit = unitMap[self.domainUnit]
            factor = PQUModule.PQU(1, self.domainUnit).getValueAs(newUnit)
            self.domainMin *= factor
            self.domainMax *= factor
            self.domainUnit = newUnit
        self.evaluated.convertUnits(unitMap)

    def toString(self, simpleString=False):
        return ("%s resonances, %s to %s. Contains %d resonances" %
                (self.evaluated, self.domainMin, self.domainMax, len(self.evaluated)))

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlString = [
            '%s<%s index="%s" domainMin="%s" domainMax="%s" domainUnit="%s">' %
            (indent, self.moniker, self.index,
             PQUModule.floatToShortestString(self.domainMin, 12), PQUModule.floatToShortestString(self.domainMax, 12),
             self.domainUnit)]
        xmlString += self.evaluated.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @property
    def BreitWigner(self):
        from .resolved import BreitWigner
        if isinstance(self.evaluated, BreitWigner): return self.evaluated
        return None

    @property
    def RMatrix(self):
        from .resolved import RMatrix
        if isinstance(self.evaluated, RMatrix): return self.evaluated
        return None

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        from .resolved import BreitWigner, RMatrix
        from .unresolved import TabulatedWidths

        xPath.append('%s[@label="%s"]' % (cls.moniker, element.get('label')))

        formClass = {
            BreitWigner.moniker: BreitWigner,
            RMatrix.moniker: RMatrix,
            TabulatedWidths.moniker: TabulatedWidths,
        }.get(element[0].tag)
        if formClass is None: raise Exception("unknown unresolved resonance form '%s'!" % element[0].tag)
        data = formClass.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        EI = cls(data=data, **getAttrs(element, required=cls._requiredAttributes))

        xPath.pop()
        return EI


# helper functions for reading in from xml:
def getBool(value):
    return {'true': True, '1': True, 'false': False, '0': False}[value]


def floatOrint(value):
    if float(value).is_integer(): return int(value)
    return float(value)


def getAttrs(element, required=(), optional=(), attributeRenames=()):
    """
    Parse all attributes, convert to type specified in required / optional lists.
    Warns if extra attributes encountered.

    :param element: Node instance with attributes to be parsed
    :param required: list of tuples: (attributeName, Type)
    :param optional: list of tuples: (attributeName, Type, default)
    :param attributeRenames: optional list of tuples for GNDS backwards-compatibility: (attributeName, attributeValDict)
    :return: dictionary of attributes converted to specified type
    """
    attrs = dict(element.items())
    for attrName, attrVals in attributeRenames:
        if attrName in attrs:
            oldVal = attrs[attrName]
            if oldVal in attrVals:
                attrs[attrName] = attrVals[oldVal]
    typed = {}
    for key, Type in required:
        assert key in attrs, 'Missing required attribute "%s"!' % key
        val = attrs.pop(key)
        if Type is bool:
            val = getBool(val)
        else:
            val = Type(val)
        typed[key] = val
    for key, Type, default in optional:
        if key in attrs:
            val = attrs.pop(key)
            if Type is bool:
                val = getBool(val)
            else:
                val = Type(val)
            typed[key] = val
    if attrs:
        print("WARNING: encountered unexpected attributes in node %s: %s" % (element.tag, attrs))
    return typed
