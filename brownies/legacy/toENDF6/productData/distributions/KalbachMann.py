# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from xData import enums as xDataEnumsModule
from xData import XYs1d as XYs1dModule
from xData import regions as regionsModule

from fudge.productData.distributions import KalbachMann as KalbachMannModule

from ... import endfFormats as endfFormatsModule
from ... import gndsToENDF6 as gndsToENDF6Module


#
# KalbachMann
#
def toENDF6(self, MT, endfMFList, flags, targetInfo):

    fSubform = self.fSubform.data
    rSubform = self.rSubform.data
    aSubform = self.aSubform
    if not aSubform.isEmptyASubform(): raise 'hell - FIXME'

    def getFlatListOfInterpolations(functional2d):
        if isinstance(functional2d, regionsModule.Regions2d):
            return [xys1d.interpolation for xys2d in functional2d for xys1d in xys2d]
        else:
            return [xys1d.interpolation for xys1d in functional2d]

    outgoingInterpolation = set(getFlatListOfInterpolations(fSubform) + getFlatListOfInterpolations(rSubform))
    if len(outgoingInterpolation) != 1:
        raise NotImplementedError("Only one outgoing interpolation supported when writing Kalbach-Mann to ENDF-6")

    EInUnit = fSubform.axes[2].unit
    EpUnit = fSubform.axes[1].unit
    fUnit = fSubform.axes[0].unit
    EInFactor = PQUModule.PQU(1, EInUnit).getValueAs('eV')
    EpFactor = PQUModule.PQU(1, EpUnit).getValueAs('eV')
    fFactor = PQUModule.PQU(1, fUnit).getValueAs('1/eV')

    LEP = {xDataEnumsModule.Interpolation.flat: 1,
           xDataEnumsModule.Interpolation.linlin: 2}[outgoingInterpolation.pop()]
    if isinstance(fSubform, regionsModule.Regions2d):
        N1 = len(fSubform)
        NE = sum([len(xys2d) for xys2d in fSubform])
        interp = fSubform[0].interpolation
        qualifier = fSubform[0].interpolationQualifier
    else:
        N1 = 1
        NE = len(fSubform)
        interp = fSubform.interpolation
        qualifier = fSubform.interpolationQualifier
        fSubform = [fSubform]
        rSubform = [rSubform]
    ENDFDataList = [endfFormatsModule.endfContLine(0, 0, 2, LEP, N1, NE)]
    interpolation = gndsToENDF6Module.gndsToENDF2PlusDInterpolationFlag(interp, qualifier)
    ENDFDataList += endfFormatsModule.endfInterpolationList([NE, interpolation])
    for f_region, r_region in zip(fSubform, rSubform):
        for i1, fAtEnergy in enumerate(f_region):
            rAtEnergy = r_region[i1]
            if not isinstance(fAtEnergy, XYs1dModule.XYs1d): raise 'hell - FIXME'
            if not isinstance(rAtEnergy, XYs1dModule.XYs1d): raise 'hell - FIXME'
            value = fAtEnergy.outerDomainValue
            if value != rAtEnergy.outerDomainValue: raise 'hell - FIXME'
            value *= EInFactor
            coefficients = []
            for i1, (Ep1, f1) in enumerate(fAtEnergy):
                Ep2, r1 = rAtEnergy[i1]
                if Ep1 != Ep2: raise 'hell - FIXME'
                coefficients += [EpFactor * Ep1, fFactor * f1, r1]
            length = len(coefficients)
            ENDFDataList += [endfFormatsModule.endfContLine(0, value, 0, 1, length, length / 3)]
            ENDFDataList += endfFormatsModule.endfDataList(coefficients)
    gndsToENDF6Module.toENDF6_MF6(MT, endfMFList, flags, targetInfo, 1, self.productFrame, ENDFDataList)

KalbachMannModule.Form.toENDF6 = toENDF6
