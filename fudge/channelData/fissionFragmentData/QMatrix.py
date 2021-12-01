# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import ancestry as ancestryModule
from xData import xDataArray as arrayModule

from fudge import suites as suitesModule

class full( arrayModule.full ) :

    pass

class QMatrix( ancestryModule.ancestry ) :

    moniker = 'QMatrix'

    def __init__( self ) :

        ancestryModule.ancestry.__init__( self )

class suite( suitesModule.suite ) :

    moniker = 'QMatrix'

    def __init__( self ) :

        suitesModule.suite.__init__( self, allowedClasses = ( QMatrix, ) )
