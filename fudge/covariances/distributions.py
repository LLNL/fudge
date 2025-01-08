# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from fudge import suites as suitesModule
from . import mixed, covarianceMatrix

"""Special classes for storing covariances for product distributions."""


class LegendreOrderCovarianceForm(ancestryModule.AncestryIO):
    """ 
    Stores covariance between energy-dependent Legendre coefficients for a reaction.
    This class contains one or more LegendreLValue sections, each section containing the matrix
    between a pair of L-values 
    """

    moniker = 'LegendreOrderCovariance'

    def __init__(self, label=None, lvalues=None):

        ancestryModule.AncestryIO.__init__(self)

        self.__label = label
        self.lvalues = lvalues or []  #: the l values of course

    def __getitem__(self, index):
        return self.lvalues[index]

    @property
    def label(self):

        return (self.__label)

    def addLegendreOrder(self, LValue):
        LValue.setAncestor(self)
        self.lvalues.append(LValue)

    def check(self, info):

        from fudge import warning

        warnings = []
        for Lvalue in self:
            Lvalue_warnings = Lvalue.check(info)
            if Lvalue_warnings:
                warnings.append(
                    warning.Context('%s L=%i vs %i' % (Lvalue.moniker, Lvalue.L1, Lvalue.L2), Lvalue_warnings))
        return warnings

    def convertUnits(self, unitMap):

        for lvalue in self: lvalue.convertUnits(unitMap)

    def fix(self, **kw):
        return []

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlString = ['%s<%s label="%s">' % (indent, self.moniker, self.label)]
        for lvalue in self.lvalues: xmlString += lvalue.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)
        form = LegendreOrderCovarianceForm(label=element.get("label"))
        # add L's to each component:
        for lValue in element:
            form.addLegendreOrder(LegendreLValue.parseNodeUsingClass(lValue, xPath, linkData, **kwargs))
        xPath.pop()
        return form


class LegendreLValue(suitesModule.Suite):
    """ 
    Represents one subsection of the Legendre coefficient covariance matrix:
    covariance between coefficients for two Legendre orders at various energies 
    """

    moniker = 'LegendreLValue'

    def __init__(self, L1, L2, frame):
        suitesModule.Suite.__init__(self, [covarianceMatrix, mixed.MixedForm])
        self.L1 = L1  #:
        self.L2 = L2  #:
        self.frame = frame  #:

    def check(self, info):

        warnings = []
        if self.L1 == self.L2:
            for form in self:
                warnings += form.check(info)
        # FIXME what about cross-terms between L-values?
        return warnings

    def fix(self, **kw):
        return []

    def toXML_strList(self, indent='', **kwargs):

        indent2 = indent + kwargs.get('incrementalIndent', '  ')

        xmlString = ['%s<%s L1="%i" L2="%i" frame="%s">' %
                     (indent, self.moniker, self.L1, self.L2, self.frame)]
        for form in self:
            xmlString += form.toXML_strList(indent2, **kwargs)
        xmlString[-1] += '</%s>' % self.moniker
        return xmlString

    @classmethod
    def parseNodeUsingClass(cls, element, xPath, linkData, **kwargs):

        xPath.append(element.tag)

        component = cls(int(element.get("L1")), int(element.get("L2")), element.get("frame"))
        for form in element:
            formCls = {
                covarianceMatrix.moniker: covarianceMatrix,
                mixed.MixedForm.moniker: mixed.MixedForm
            }.get(form.tag)
            if formCls is None:
                raise TypeError("Encountered unknown covariance type '%s'" % form.tag)
            component.add(formCls.parseNodeUsingClass(form, xPath, linkData, **kwargs))
        xPath.pop()
        return component
