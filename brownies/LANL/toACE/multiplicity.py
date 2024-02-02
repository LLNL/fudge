# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module adds the method toACE to the classes in the fudge.productData.multiplicity module.
"""

from xData import enums as xDataEnumsModule

from fudge.productData import multiplicity as multiplicityModule


#
#   Polynomial1d multiplicity.
#
def toACE(self):

    factor = 1
    coefficients = []
    for coefficient in self:
        coefficients.append(coefficient / factor)
    return [1, len(self)] + coefficients

multiplicityModule.Polynomial1d.toACE = toACE

#
#   XYs1d multiplicity.
#
def toACE(self):

    interpolation = 0
    if self.interpolation == xDataEnumsModule.Interpolation.flat:
        interpolation = 1
    elif self.interpolation == xDataEnumsModule.Interpolation.linlin:
        interpolation = 2
    if interpolation == 0: raise Exception('Interpolation "%s" not supported' % self.interpolation)

    return [interpolation, 0, len(self)] + [E1 for E1, m1 in self] + [m1 for E1, m1 in self]

multiplicityModule.XYs1d.toACE = toACE

#
#   Regions1d multiplicity.
#
def toACE(self):

    interpolations = set([region.interpolation for region in self])
    if len(interpolations) != 1:
        raise Exception('Multiple different interpolations "%s" not supported' % interpolations)

    interpolation = interpolations.pop()
    if interpolation == xDataEnumsModule.Interpolation.flat:
        interpolation = 1
    elif interpolation == xDataEnumsModule.Interpolation.linlin:
        interpolation = 2
    else:
        raise Exception('Interpolation "%s" not supported' % interpolation)

    nPoints = sum([len(region) for region in self])
    xys = [xy for region in self for xy in region]
    xs, ys = zip(*xys)
    return [interpolation, 0, nPoints] + list(xs) + list(ys)

multiplicityModule.Regions1d.toACE = toACE
