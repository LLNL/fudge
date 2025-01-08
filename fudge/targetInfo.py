# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""This module defines the TargetInfo class, used to store isotopic abundances for TNSL targets."""

from LUPY import ancestry as ancestryModule
from fudge import suites as suitesModule


class TargetInfo(ancestryModule.Ancestry):
    moniker = "targetInfo"

    def __init__(self):
        super().__init__()
        self.__isotopicAbundances = IsotopicAbundances()

    @property
    def isotopicAbundances(self):
        return self.__isotopicAbundances

    def toXML_strList(self, indent='', **kwargs):
        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        xml = [f'{indent}<{self.moniker}>']
        xml += self.isotopicAbundances.toXML_strList(indent2, **kwargs)
        xml[-1] += f'</{self.moniker}>'
        return xml

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        xPath.append(node.tag)

        instance = cls()
        instance.isotopicAbundances.parseNode(
            node.find(IsotopicAbundances.moniker), xPath, linkData, **kwargs)
        xPath.pop()
        return instance


class IsotopicAbundances(ancestryModule.Ancestry):
    moniker = 'isotopicAbundances'

    def __init__(self):
        super().__init__()
        self.__chemicalElements = ChemicalElements()

    @property
    def chemicalElements(self):
        return self.__chemicalElements

    def toXML_strList(self, indent='', **kwargs):
        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        xml = [f'{indent}<{self.moniker}>']
        xml += self.chemicalElements.toXML_strList(indent2, **kwargs)
        xml[-1] += f'</{self.moniker}>'
        return xml

    def parseNode(self, node, xPath, linkData, **kwargs):
        xPath.append(node.tag)
        self.chemicalElements.parseNode(node.find(ChemicalElements.moniker), xPath, linkData, **kwargs)
        xPath.pop()


class ChemicalElement(ancestryModule.Ancestry):
    moniker = "chemicalElement"

    def __init__(self, symbol):
        super().__init__()
        self.symbol = symbol
        self.__nuclides = Nuclides()

    @property
    def label(self):
        # FIXME: hack to support making Elements inherit from fudge.suites.Suite
        return self.symbol

    @property
    def nuclides(self):
        return self.__nuclides

    def toXML_strList(self, indent='', **kwargs):
        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        attrs = ''
        xml = [f'{indent}<{self.moniker} symbol="{self.symbol}"{attrs}>']
        xml += self.nuclides.toXML_strList(indent2, **kwargs)
        xml[-1] += f'</{self.moniker}>'
        return xml

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        xPath.append(node.tag)

        instance = cls(symbol=node.get("symbol"))
        instance.nuclides.parseNode(node.find(Nuclides.moniker), xPath, linkData, **kwargs)
        xPath.pop()
        return instance


class ChemicalElements(suitesModule.ExclusiveSuite):
    moniker = 'chemicalElements'

    def __init__(self):

        suitesModule.ExclusiveSuite.__init__(self, [ChemicalElement])


class Nuclide(ancestryModule.Ancestry):
    moniker = "nuclide"

    def __init__(self, pid, atomFraction):
        super().__init__()
        self.pid = pid
        self.atomFraction = atomFraction

    @property
    def label(self):
        # FIXME: hack to support making Nuclides inherit from fudge.suites.Suite
        return self.pid

    @property
    def atomFraction(self):
        return self.__atomFraction

    @atomFraction.setter
    def atomFraction(self, value):
        value = float(value)
        assert 0 < value <= 1
        self.__atomFraction = value

    def toXML_strList(self, indent='', **kwargs):
        indent2 = indent + kwargs.get('incrementalIndent', '  ')
        xml = [f'{indent}<{self.moniker} pid="{self.pid}" atomFraction="{self.atomFraction}"/>']
        return xml

    @classmethod
    def parseNodeUsingClass(cls, node, xPath, linkData, **kwargs):
        xPath.append(node.tag)

        instance = cls(pid=node.get("pid"), atomFraction=node.get("atomFraction"))
        xPath.pop()
        return instance


class Nuclides(suitesModule.ExclusiveSuite):
    moniker = 'nuclides'

    def __init__(self):

        suitesModule.ExclusiveSuite.__init__(self, [Nuclide])
