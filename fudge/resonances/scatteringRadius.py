# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU
from LUPY import ancestry as ancestryModule

from xData import constant as constantModule
from xData import regions as regionsModule
from xData import XYs1d as XYs1dModule


class BaseRadius(ancestryModule.AncestryIO):
    """
    Base class for ScatteringRadius and HardSphereRadius. Contains a radius that may be constant or energy-dependent.

    ScatteringRadius and HardSphereRadius forms require a label attribute when they appear inside a 'resonances' node,
    but not when they appear inside the resolved or unresolved region. This class handles both cases.
    """

    def __init__(self, form):

        ancestryModule.AncestryIO.__init__(self)
        self.evaluated = form

    def __eq__(self, other):
        if not isinstance(other, ScatteringRadius): return False
        return self.evaluated == other.evaluated

    def __str__(self):
        return self.evaluated.moniker

    def __bool__(self):

        return bool(self.evaluated)

    def __getitem__(self, key):
        if key == self.evaluated.label:
            return self.evaluated
        raise KeyError(key)

    __nonzero__ = __bool__

    @property
    def evaluated(self):
        return self.__evaluated

    @evaluated.setter
    def evaluated(self, form):
        if not isinstance(form, (constantModule.Constant1d, XYs1dModule.XYs1d, regionsModule.Regions1d)):
            raise TypeError(f"{self.moniker} form must be constant1d, XYs1d or Regions1d, but received {type(form)}")
        form.setAncestor(self)
        self.__evaluated = form

    @property
    def form(self):
        print("WARNING: scatteringRadius.form is deprecated. Use scatteringRadius.evaluated instead")
        return self.evaluated

    def toString(self, simpleString=False):
        """Returns a string representation of self. simpleString option included for compatibility."""
        return str(self)

    def copy(self):

        return self.__class__(self.evaluated.copy())

    def check(self, info):
        """
        Checks that the scattering radius is within a factor of 3 of the expected scattering radius
        of AP = 0.123 * self.target.getMass('amu')**(1./3.) + 0.08

        This default AP is roughly the nuclear radius, so getting within a factor of 3 of the default
        shouldn't be a problem.

        A similar test is included in PSYCHE, but PSYCHE's cannot handle LRF=7
        :param info:
        :return:
        """
        from fudge import warning
        warnings = []
        target = info['reactionSuite'].PoPs[info['reactionSuite'].target]
        if target.id in info['reactionSuite'].PoPs.aliases:
            target = info['reactionSuite'].PoPs[target.pid]
        expectedAP = 10.0 * (0.123 * target.getMass('amu') ** (1. / 3.) + 0.08)  # expected radius in fm
        factor = 3.0
        if self.isEnergyDependent():
            egrid = self.evaluated.domainGrid
            APs = self.getValueAs('fm', energy_grid=egrid)
            for iE, AP in enumerate(APs):
                if AP / expectedAP > factor or AP / expectedAP < 1. / factor:
                    warning.BadScatteringRadius(factor=factor, gotAP=AP, expectedAP=expectedAP, E=egrid[iE])
        else:
            AP = self.evaluated.value
            if self.evaluated.rangeUnit != 'fm':
                AP *= PQU.PQU(1, self.evaluated.rangeUnit).getValueAs('fm')
            if AP / expectedAP > factor or AP / expectedAP < 1. / factor:
                warning.BadScatteringRadius(factor=factor, gotAP=AP, expectedAP=expectedAP)
        return warnings

    def convertUnits(self, unitMap):

        self.evaluated.convertUnits(unitMap)

    def isEnergyDependent(self):
        return isinstance(self.evaluated, XYs1dModule.XYs1d)

    def getValueAs(self, unit, energy_grid=None):
        if self.isEnergyDependent():
            if energy_grid is None:
                raise NameError("Missing: energy_grid to evaluate E-dependent scattering radius")
            energy_unit = self.evaluated.axes[-1].unit
            xScale = self.evaluated.domainUnitConversionFactor(energy_unit)
            yScale = self.evaluated.rangeUnitConversionFactor(unit)
            return [yScale * self.evaluated.evaluate(xScale * e) for e in energy_grid]
        else:
            oldUnit = self.evaluated.rangeUnit
            factor = PQU.PQU(1, oldUnit).getValueAs(unit)
            return self.evaluated.value * factor

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xml = ['%s<%s>' % (indent, self.moniker)]
        xml += self.evaluated.toXML_strList(indent2, **kwargs)
        xml[-1] += '</%s>' % self.moniker
        return xml

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)

        childClass = {
            constantModule.Constant1d.moniker: constantModule.Constant1d,
            XYs1dModule.XYs1d.moniker: XYs1dModule.XYs1d,
            regionsModule.Regions1d.moniker: regionsModule.Regions1d,
        }[element[0].tag]

        form = childClass.parseNodeUsingClass(element[0], xPath, linkData, **kwargs)
        SR = cls(form)

        xPath.pop()

        return SR


class ScatteringRadius(BaseRadius):
    """
    The ScatteringRadius determines scattering length and is used for computing penetrability, shift function and phase shift.
    The Resonances class must define a scatteringRadius. It may also appear in a BreitWigner, unresolved TabulatedWidths,
    or may be overridden for a specific ResonanceReaction or Channel.
    """
    moniker = 'scatteringRadius'


class HardSphereRadius(BaseRadius):
    """
    The HardSphereRadius may be used to override the ScatteringRadius. This value is only used for computing the phase shift.
    May appear everywhere that a ScatteringRadius may be defined.
    """
    moniker = 'hardSphereRadius'
