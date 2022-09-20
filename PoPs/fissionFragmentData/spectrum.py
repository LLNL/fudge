# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import XYs1d as XYs1dModule

from .. import suite as suiteModule

class XYs1d( XYs1dModule.XYs1d ) :

    pass

class Suite( suiteModule.SortedSuite ) :

    moniker = 'spectrum'

    def __init__( self ) :

        suiteModule.SortedSuite.__init__( self, allowedClasses = ( XYs1d, ) )
