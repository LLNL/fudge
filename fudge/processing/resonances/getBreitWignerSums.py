# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import numpy
from . import _getBreitWignerSums

# #### Functions & Classes ###########################################

def getBreitWignerSums(E, Eres, captureWidth, elasticWidth, fissionWidth, shiftFactorAtRes,
                       penetrationFactor, shiftFactor, phi, approximation):
    """

    :param E: incident energies (numpy array)
    :param Eres: resonance energies
    :param captureWidth: resonance capture widths
    :param elasticWidth: resonance neutron widths
    :param fissionWidth: resonance fission widths (0 for non-fissile)
    :param shiftFactorAtRes: shift factor computed at Eres
    :param penetrationFactor: computed at E
    :param shiftFactor: computed at E
    :param phi: computed at E
    :param approximation: 'SingleLevel' or 'MultiLevel'
    :return:
    """
    # .... Check arguments:
    for dat in (E,Eres,captureWidth,elasticWidth,fissionWidth,shiftFactorAtRes,penetrationFactor,shiftFactor,phi):
        if type(dat) != numpy.ndarray: raise TypeError("All inputs must be numpy arrays!")
        if dat.dtype != numpy.float64: raise TypeError("All inputs must be float arrays")

    # check dimensions. Row or column is fine for incident energy:
    NE = len(E)
    nRes = len(Eres)
    for arr,name in zip((E, penetrationFactor, shiftFactor, phi),('E','penetrationFactor','shiftFactor','phi')):
        if arr.shape not in ((NE,),(NE,1)):
            raise TypeError("%s has wrong shape: %s, should be %s" % (name, arr.shape, (NE,)))
    for arr,name in zip((Eres, captureWidth, elasticWidth, fissionWidth, shiftFactorAtRes),
            ('Eres','captureWidth','elasticWidth','fissionWidth','shiftFactorAtRes')):
        if arr.shape != (nRes,):
            raise TypeError("%s has wrong shape: %s, should be %s" % (name, arr.shape, (nRes,)))

    # resonance energy & width may not be contiguous arrays:
    Eres = numpy.ascontiguousarray(Eres)
    captureWidth = numpy.ascontiguousarray(captureWidth)
    elasticWidth = numpy.ascontiguousarray(elasticWidth)
    fissionWidth = numpy.ascontiguousarray(fissionWidth)
    shiftFactorAtRes = numpy.ascontiguousarray(shiftFactorAtRes)

    approx = {'SingleLevel': 1, 'MultiLevel': 2}[ approximation ]
    # .... Call C extension function
    capture,elastic,fission = _getBreitWignerSums.getBreitWignerSums(E,Eres[:],captureWidth[:],elasticWidth[:],
            fissionWidth[:],shiftFactorAtRes[:],penetrationFactor,shiftFactor,phi,approx)
    return (capture,elastic,fission)
