# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from xData import enums as xDataEnumsModule
from fudge.productData.distributions import LLNL_angularEnergy as LLNL_angularEnergyModule
from fudge.productData.distributions import angularEnergy as angularEnergyModule

from ... import gndsToENDF6 as gndsToENDF6Module
from ... import endfFormats as endfFormatsModule

#
# LLNLAngularEnergyForm
#
def toENDF6(self, MT, endfMFList, flags, targetInfo):

    angularForm = self.angularSubform.data
    angularEnergyForm = self.angularEnergySubform.data
    if len(angularForm) != len(angularEnergyForm):
        raise Exception('len(angularForm) = %s != len(angularEnergyForm) = %s.', (len(angularForm), len(angularEnergyForm)))

    axes = angularEnergyModule.defaultAxes(angularForm.domainUnit)

    xys3d = angularEnergyModule.XYs3d(axes=axes, interpolation=xDataEnumsModule.Interpolation.linlin, 
            interpolationQualifier=xDataEnumsModule.InterpolationQualifier.unitBase)

    energyConversionFactor = PQUModule.PQU(1, angularForm.axes[-1].unit ).getValueAs('eV')
    for energyInIndex, P_ofMus in enumerate(angularForm):
        P_ofEpGiveEMus = angularEnergyForm[energyInIndex]
        if len(P_ofMus) != len(P_ofEpGiveEMus):
            raise Exception('For incident energy %.7e, len(P_ofMus) %s != len(P_ofEpGiveEMus) = %s' % 
                    (P_ofMus.outerDomainValue, len(P_ofMus), len(P_ofEpGiveEMus)))

        xys2d = angularEnergyModule.XYs2d(axes=axes, interpolation=xDataEnumsModule.Interpolation.linlin, outerDomainValue=P_ofMus.outerDomainValue)
        for muIndex, P_ofMu in enumerate(P_ofMus):
            mu, muProbability = P_ofMu
            P_ofEpGiveEMu = muProbability * P_ofEpGiveEMus[muIndex]
            xys1d = angularEnergyModule.XYs1d(data=P_ofEpGiveEMu, axes=axes, outerDomainValue=mu)
            xys2d.append(xys1d)
        xys3d.append(xys2d)

    angularEnergy = angularEnergyModule.Form('temp', xDataEnumsModule.Frame.lab, xys3d)
    angularEnergy.toENDF6(MT, endfMFList, flags, targetInfo)

LLNL_angularEnergyModule.LLNLAngularEnergyForm.toENDF6 = toENDF6
