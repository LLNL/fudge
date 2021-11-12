# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import ancestry as ancestryModule
from xData import xDataArray as arrayModule

from .. import suite as suiteModule

class full( arrayModule.full ) :

    pass

class QMatrix( ancestryModule.ancestry ) :

    moniker = 'QMatrix'

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

class suite( suiteModule.sortedSuite ) :

    moniker = 'QMatrix'

    def __init__( self ) :

        suiteModule.sortedSuite.__init__( self, allowedClasses = ( QMatrix, ) )
