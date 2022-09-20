# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module adds the method toACE to the reaction class.
"""

from fudge import outputChannel as outputChannelModule


def toACE( self, cdf_style, MTData, MT, verbose ) :

    for product in self : product.toACE( cdf_style, MTData, MT, verbose )

outputChannelModule.OutputChannel.toACE = toACE
