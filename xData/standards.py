# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains standard xData types.
"""

__metaclass__ = type

import sys

class types :

    integer32Token = 'Integer32'
    float32Token = 'Float32'
    float64Token = 'Float64'
    validTypes = ( integer32Token, float32Token, float64Token )

noExtrapolationToken = 'noExtrapolation'
flatExtrapolationToken = 'flatExtrapolation'

validExtrapolations = ( noExtrapolationToken, flatExtrapolationToken )

class floats :

    epsilonMultiplier = 4
    epsilon = sys.float_info.epsilon * epsilonMultiplier

class frames :

    labToken = 'lab'
    centerOfMassToken = 'centerOfMass'
    productToken = 'product'
    allowedFrames = ( labToken, centerOfMassToken )

class interpolation :

    linlinToken = 'lin-lin'
    linlogToken = 'lin-log'
    loglinToken = 'log-lin'
    loglogToken = 'log-log'
    flatToken = 'flat'
    chargedParticleToken = 'charged-particle'
    defaultAllowedTokens = [ linlinToken, linlogToken, loglinToken, loglogToken, flatToken ]

# Should add a interpolationQualifier class with the following and allowed.

    noneQualifierToken = ''
    unitBaseToken = 'unitbase'
    unitBaseUnscaledToken = 'unitbase-unscaled'
    correspondingPointsToken = 'correspondingPoints'
    correspondingPointsUnscaledToken = 'correspondingPoints-unscaled'
# We also need direct, correspondingEnergies and correspondingEnergies-unscaled
