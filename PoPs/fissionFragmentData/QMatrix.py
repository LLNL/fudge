# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from LUPY import ancestry as ancestryModule

from xData import xDataArray as arrayModule

from .. import suite as suiteModule

class Full( arrayModule.Full ) :

    pass

class QMatrix(ancestryModule.AncestryIO):

    moniker = 'QMatrix'

    def __init__( self ) :

        ancestryModule.AncestryIO.__init__( self )

class Suite( suiteModule.SortedSuite ) :

    moniker = 'QMatrix'

    def __init__( self ) :

        suiteModule.SortedSuite.__init__( self, allowedClasses = ( QMatrix, ) )
