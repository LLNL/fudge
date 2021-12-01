# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from fudge.reactionData.doubleDifferentialCrossSection import base as baseModule

from ... import endfFormats as endfFormatsModule
from ... import gndsToENDF6 as gndsToENDF6Module

#
# XYs1d.
#
def toENDF6( self, MT, endfMFList, targetInfo ) :

    if( MT == 502 ) : MT = self.ancestor.ENDFMT
    Z = targetInfo['ZA'] / 1000

    endfInterpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag( self.interpolation )
    data = []
    for xy in self.copyDataToXYs( ) : data += xy
    endfMFList[27][MT] = [ endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'],  0, 0, 0, 0 ),
                           endfFormatsModule.endfContLine( 0, Z,  0, 0, 1, len( data ) / 2 ) ] + \
                           endfFormatsModule.endfInterpolationList( ( len( data ) / 2, endfInterpolation ) ) + \
                           endfFormatsModule.endfDataList( data ) + [ endfFormatsModule.endfSENDLineNumber( ) ]

baseModule.XYs1d.toENDF6 = toENDF6
#
# regions1d
#
def toENDF6( self, MT, endfMFList, targetInfo ) :

    if( MT == 502 ) : MT = self.ancestor.ENDFMT
    Z = targetInfo['ZA'] / 1000

    endfInterpolation, data = [], []
    counter, lastX, lastY = 0, None, None
    for region in self :
        ENDFInterpolation = gndsToENDF6Module.gndsToENDFInterpolationFlag( region.interpolation )
        regionData = region.copyDataToXYs( )
        if( lastX is not None ) :
            if( lastY == regionData[0][1] ) : del regionData[0]
        counter += len( regionData )
        endfInterpolation.append( counter )
        endfInterpolation.append( ENDFInterpolation )
        for xy in regionData : data += xy
        lastX, lastY = regionData[-1]

    endfMFList[27][MT] = [ endfFormatsModule.endfHeadLine( targetInfo['ZA'], targetInfo['mass'],  0, 0, 0, 0 ),
                           endfFormatsModule.endfContLine( 0, Z,  0, 0, len( endfInterpolation ) / 2, len( data ) / 2 ) ] + \
                           endfFormatsModule.endfInterpolationList( endfInterpolation ) + \
                           endfFormatsModule.endfDataList( data ) + [ endfFormatsModule.endfSENDLineNumber( ) ]

baseModule.regions1d.toENDF6 = toENDF6
