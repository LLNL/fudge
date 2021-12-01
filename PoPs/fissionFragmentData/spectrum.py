# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import XYs as XYsModule

from .. import suite as suiteModule

class XYs1d( XYsModule.XYs1d ) :

    pass

class suite( suiteModule.sortedSuite ) :

    moniker = 'spectrum'

    def __init__( self ) :

        suiteModule.sortedSuite.__init__( self, allowedClasses = ( XYs1d, ) )
