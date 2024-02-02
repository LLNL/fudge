# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import abc
import numpy

from LUPY import ancestry as ancestryModule
from PoPs.quantities import quantity as quantityModule

"""
Defines incident-energy-dependent functions representing the contribution to the resolved region cross section
from resonances external to the evaluated set.
The external R-Matrix is added to the R-Matrix diagonal during resonance reconstruction.
"""


class ExternalRMatrix(ancestryModule.AncestryIO, metaclass=abc.ABCMeta):
    """
    Abstract base class inherited by the Froehner and SAMMY classes.
    """
    moniker = 'externalRMatrix'

    def __init__(self, **kwargs):

        super().__init__()
        provided_terms = set(kwargs.keys())
        required_terms = {'singularityEnergyBelow', 'singularityEnergyAbove'}
        if not required_terms.issubset(provided_terms):
            missing = required_terms.difference(provided_terms)
            raise AttributeError("%s external R-Matrix is missing required terms: %s" % (self.type, ", ".join(missing)))
        extra = provided_terms.difference(self.ancestryMembers)
        if extra:
            raise AttributeError("%s external R-Matrix received unexpected terms: %s" % (self.type, ", ".join(extra)))
        self._terms = kwargs

    @property
    @abc.abstractmethod
    def type(self): pass

    @property
    @abc.abstractmethod
    def terms(self): pass

    @abc.abstractmethod
    def evaluate(self, energies): pass

    def getTerm(self, key, unit):
        result = self.terms.get(key)
        if result is None:
            return 0
        return result.float(unit)

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + '  '
        xmlString = ['%s<%s type="%s">' % (indent, self.moniker, self.type)]
        for key in self.ancestryMembers:
            term = self.terms.get(key)
            if term is not None:
                xmlString += term.toXML_strList(indent=indent2, **kwargs)
        xmlString[-1] += ('</%s>' % self.moniker)

        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)

        terms = {term.get("label"): quantityModule.Double.parseNodeUsingClass(term, xPath, linkData, **kwargs)
                 for term in element.findall("double")}

        class_ = {
            'Froehner': Froehner,
            'SAMMY': SAMMY
        }[element.get("type")]

        result = class_(**terms)            # FIXME2, this is not using cls.

        xPath.pop()

        return result


class Froehner(ExternalRMatrix):
    """
    Froehner's external R-Matrix parametrization.
    """
    ancestryMembers = ('averageRadiationWidth', 'constantExternalR', 'poleStrength', 'singularityEnergyBelow',
                       'singularityEnergyAbove')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @property
    def type(self):
        return "Froehner"

    @property
    def terms(self):
        return self._terms

    def evaluate(self, energies):
        """
        Evaluate Fröhner's external R-Matrix parametrization at the given energy or energies.
        @param energies: single energy or numpy array of energies
        @return: tuple(real part, imaginary part)
        """
        R0 = self.getTerm('constantExternalR', '')
        sc = self.getTerm('poleStrength', '')
        Gamma = self.getTerm('averageRadiationWidth', 'eV')
        Edown = self.getTerm('singularityEnergyBelow', 'eV')
        Eup = self.getTerm('singularityEnergyAbove', 'eV')
        Ebar = (Eup + Edown) / 2
        I = Eup - Edown
        realTerm = R0 + 2 * sc * numpy.arctan2(energies - Ebar, I/2)
        imaginaryTerm = (Gamma * I / 4) / (I**2 / 4 - (energies - Ebar)**2)
        return realTerm, imaginaryTerm

    def convertUnits(self, unitMap):
        energyUnit = str(self.terms['singularityEnergyBelow'].unit)
        if energyUnit not in unitMap:
            return

        newUnit = unitMap[energyUnit]
        for key in ('singularityEnergyBelow', 'singularityEnergyAbove', 'averageRadiationWidth'):
            value = self.terms[key]
            self.terms[key] = quantityModule.Double(key, value.float(newUnit), newUnit)


class SAMMY(ExternalRMatrix):
    """
    External R-Matrix parametrization from SAMMY
    """
    ancestryMembers = ('constantExternalR', 'linearExternalR', 'quadraticExternalR',
                       'constantLogarithmicCoefficient', 'linearLogarithmicCoefficient',
                       'singularityEnergyBelow', 'singularityEnergyAbove')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @property
    def type(self):
        return "SAMMY"

    @property
    def terms(self):
        return self._terms

    def evaluate(self, energies):
        """
        Evaluate the SAMMY external R-Matrix parametrization at the given energy or energies.
        @param energies: single energy or numpy array of energies
        @return: tuple(real part, imaginary part)
        """
        Rcon = self.getTerm('constantExternalR', '')
        Rlin = self.getTerm('linearExternalR', '1/eV')
        Rquad = self.getTerm('quadraticExternalR', '1/eV**2')
        scon = self.getTerm('constantLogarithmicCoefficient', '')
        slin = self.getTerm('linearLogarithmicCoefficient', '1/eV')
        Edown = self.getTerm('singularityEnergyBelow', 'eV')
        Eup = self.getTerm('singularityEnergyAbove', 'eV')
        logTerm = numpy.log((Eup - energies) / (energies - Edown))
        realTerm = Rcon + Rlin * energies + Rquad * energies**2 - slin * (Eup - Edown) - (scon + slin * energies) * logTerm
        return realTerm, numpy.zeros_like(energies)

    def convertUnits(self, unitMap):
        energyUnit = str(self.terms['singularityEnergyBelow'].unit)
        if energyUnit not in unitMap:
            return

        newUnit = unitMap[energyUnit]
        inverse = f"1/{newUnit}"
        inverseSq = f"1/{newUnit}**2"
        for key, unit in zip(self.ancestryMembers,
                             ('', inverse, inverseSq, '', inverse, newUnit, newUnit)):
            value = self.terms[key]
            self.terms[key] = quantityModule.Double(key, value.float(unit), unit)
